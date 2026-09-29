using Pkg
Pkg.activate("./")
using rosetta_distribution_opf
import PowerModelsDistribution
import InfrastructureModels
using Ipopt
using JuMP
using Statistics

const PMD = PowerModelsDistribution
const IM  = InfrastructureModels
PMD.silence!()

ipopt_solver = JuMP.optimizer_with_attributes(
    Ipopt.Optimizer,
    "print_level"          => 0,
    "sb"                   => "yes",
    "warm_start_init_point"=> "yes",
    "max_iter"             => 2000
)

data_path = "./data/ENWL_4w_Network1_Feeder1/Master.dss"

# ─────────────────────────────────────────────────────────
# NETWORK LOADER
# ─────────────────────────────────────────────────────────
function load_base_network(data_path; sbase=1.0, load_multiplier=1.0)
    data_eng = PMD.parse_file(data_path, transformations=[PMD.transform_loops!])
    data_eng["settings"]["sbase_default"] = sbase

    data_math = PMD.transform_data_model(
        data_eng, multinetwork=false, kron_reduce=false, phase_project=false
    )

    # Physical voltage bounds
    for (i, bus) in data_math["bus"]
        bus["vmin"] = [0.5, 0.5, 0.5, 0]
        bus["vmax"] = [1.5, 1.5, 1.5, 1.5]
    end

    # Unconstrain slack generator
    for (i, gen) in data_math["gen"]
        gen["pmax"] =  [1e4, 1e4, 1e4]
        gen["pmin"] = -[1e4, 1e4, 1e4]
        gen["qmax"] =  [1e4, 1e4, 1e4]
        gen["qmin"] = -[1e4, 1e4, 1e4]
    end

    # Scale loads
    for (i, load) in data_math["load"]
        load["pd"] *= load_multiplier
        load["qd"] *= load_multiplier
    end

    return data_math
end

# ─────────────────────────────────────────────────────────
# PV PLACEMENT
# ─────────────────────────────────────────────────────────
function add_pv!(data_math;
        p_solar_kw     = 4.0,    # positive = generation in PMD convention
        s_inverter_kva = 5.0,
        spacing        = 20,
        curtailable    = true)

    @assert p_solar_kw <= s_inverter_kva "Inverter kVA must be >= solar kW"

    sbase   = data_math["settings"]["sbase_default"]
    p_pu    = p_solar_kw    / sbase
    s_pu    = s_inverter_kva / sbase
    q_limit = sqrt(max(0.0, s_pu^2 - p_pu^2))

    # Correct source bus detection — handles "sourcebus" string name
    source_buses = Set([i for (i, bus) in data_math["bus"] if bus["bus_type"] == 3])

    bus_ids = sort(
        collect(keys(data_math["bus"])),
        by = x -> (x ∈ source_buses ? -1 : tryparse(Int, x) === nothing ? 0 : parse(Int, x))
    )

    pv_ids = String[]
    for i in 1:spacing:length(bus_ids)
        target_bus = bus_ids[i]
        target_bus ∈ source_buses && continue

        gen_id = string(length(data_math["gen"]) + 1)
        data_math["gen"][gen_id] = deepcopy(data_math["gen"]["1"])
        gen = data_math["gen"][gen_id]

        gen["gen_bus"] = parse(Int, target_bus)
        gen["type"]    = "PV"
        gen["name"]    = "pv_bus_$(target_bus)"

        if curtailable
            # P free between 0 (curtailed) and p_pu (full output)
            gen["pmax"] =  p_pu * ones(3)
            gen["pmin"] =  zeros(3)
            # Q envelope uses full inverter rating
            gen["qmax"] =  s_pu * ones(3)
            gen["qmin"] = -s_pu * ones(3)
        else
            gen["pmax"] =  p_pu * ones(3)
            gen["pmin"] =  p_pu * ones(3)   # fixed — use with caution
            gen["qmax"] =  q_limit * ones(3)
            gen["qmin"] = -q_limit * ones(3)
        end

        gen["cost"] = [-1.0 0.0]   # incentivise generation
        push!(pv_ids, gen_id)
    end

    println("  Added $(length(pv_ids)) PV units at buses: $([bus_ids[i] for i in 1:spacing:length(bus_ids) if bus_ids[i] ∉ source_buses])")
    return pv_ids
end

# ─────────────────────────────────────────────────────────
# STATCOM PLACEMENT
# ─────────────────────────────────────────────────────────
function add_statcoms!(data_math;
        q_statcom_kvar = 15.0,
        spacing        = 10)

    sbase = data_math["settings"]["sbase_default"]
    q_pu  = q_statcom_kvar / sbase

    load_ids = sort(collect(keys(data_math["load"])), by=x->parse(Int,x))
    statcom_ids = String[]

    for i in 1:spacing:length(load_ids)
        load   = data_math["load"][load_ids[i]]
        gen_id = string(length(data_math["gen"]) + 1)

        data_math["gen"][gen_id] = deepcopy(data_math["gen"]["1"])
        gen = data_math["gen"][gen_id]

        gen["gen_bus"] = copy(load["load_bus"])
        gen["type"]    = "STATCOM"
        gen["name"]    = "statcom_load_$(load_ids[i])"

        # Pure reactive — zero real power
        gen["pmax"] =  zeros(3)
        gen["pmin"] =  zeros(3)
        gen["qmax"] =  q_pu * ones(3)
        gen["qmin"] = -q_pu * ones(3)

        gen["cost"] = [0.01 0.0]   # small penalty to avoid unnecessary Q dispatch
        push!(statcom_ids, gen_id)
    end

    println("  Added $(length(statcom_ids)) STATCOMs")
    return statcom_ids
end

# ─────────────────────────────────────────────────────────
# SOLVE + REPORT
# ─────────────────────────────────────────────────────────
function solve_and_report(data_math, label)
    PMD.add_start_vrvi!(data_math)
    model = PMD.instantiate_mc_model(
        data_math,
        PMD.IVRENPowerModel,
        PMD.build_mc_opf
    )
    result = PMD.optimize_model!(model, optimizer=ipopt_solver)
    status = result["termination_status"]

    println("\n  [$label] Status: $status")

    if string(status) in ["LOCALLY_SOLVED", "OPTIMAL"]
        sol = result["solution"]

        vm_all = Float64[]
        for (b, bus) in data_math["bus"]
            bus["bus_type"] == 3 && continue          # skip slack
            !haskey(sol["bus"], b) && continue
            for p in 1:3
                vr_val = get(sol["bus"][b], "vr", zeros(4))[p]
                vi_val = get(sol["bus"][b], "vi", zeros(4))[p]
                push!(vm_all, abs(vr_val + im * vi_val))
            end
        end

        isempty(vm_all) && return result

        println("    Min  voltage : $(round(minimum(vm_all), digits=4)) pu")
        println("    Mean voltage : $(round(mean(vm_all),    digits=4)) pu")
        println("    Max  voltage : $(round(maximum(vm_all), digits=4)) pu")
        println("    Buses > 1.1 : $(count(v -> v > 1.1, vm_all))")
        println("    Buses < 0.90 : $(count(v -> v < 0.90, vm_all))")
    else
        println("    WARNING: Did not solve — check model feasibility")
        # Print Ipopt's raw output for debugging
        println("    Objective: $(get(result, "objective", "N/A"))")
    end

    return result
end

# ─────────────────────────────────────────────────────────
# FOUR CASES
# ─────────────────────────────────────────────────────────
println("="^55)
println(" CASE 1: Baseline")
println("="^55)
dm1 = load_base_network(data_path, load_multiplier=1.0)
r1  = solve_and_report(dm1, "Baseline")

println("\n" * "="^55)
println(" CASE 2: PV only")
println("="^55)
dm2 = load_base_network(data_path, load_multiplier=1.0)
add_pv!(dm2; p_solar_kw=4.0, s_inverter_kva=5.0, spacing=20)
r2  = solve_and_report(dm2, "PV only")

println("\n" * "="^55)
println(" CASE 3: STATCOM only")
println("="^55)
dm3 = load_base_network(data_path, load_multiplier=1.0)
add_statcoms!(dm3; q_statcom_kvar=20.0, spacing=2)
r3  = solve_and_report(dm3, "STATCOM only")

println("\n" * "="^55)
println(" CASE 4: PV + STATCOM")
println("="^55)
dm4 = load_base_network(data_path, load_multiplier=1.0)
add_pv!(dm4;       p_solar_kw=4.0, s_inverter_kva=5.0, spacing=20)
add_statcoms!(dm4; q_statcom_kvar=20.0, spacing=2)
r4  = solve_and_report(dm4, "PV + STATCOM")