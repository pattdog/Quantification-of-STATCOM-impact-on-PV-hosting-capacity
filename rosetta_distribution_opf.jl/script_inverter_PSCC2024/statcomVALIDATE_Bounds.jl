#=using Pkg
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
    "print_level"           => 0,
    "sb"                    => "yes",
    "warm_start_init_point" => "yes",
    "max_iter"              => 50000
)

data_path = "./data/ENWL_4w_Network1_Feeder1/Master.dss"

# ═══════════════════════════════════════════════════════════════
# NETWORK LOADER
# ═══════════════════════════════════════════════════════════════
function load_base_network(data_path;
        load_multiplier = 1.0,
        enforce_bounds  = false)

    data_eng  = PMD.parse_file(data_path, transformations=[PMD.transform_loops!])
    data_math = PMD.transform_data_model(
        data_eng, multinetwork=false, kron_reduce=false, phase_project=false
    )

    for (i, bus) in data_math["bus"]
        if bus["bus_type"] == 3
            bus["vmin"] = zeros(4)
            bus["vmax"] = 10.0 * ones(4)
        elseif enforce_bounds
            bus["vmin"] = [0.90, 0.90, 0.90, 0.00]
            bus["vmax"] = [1.10, 1.10, 1.10, 1.10]
        else
            bus["vmin"] = [0.50, 0.50, 0.50, 0.00]
            bus["vmax"] = [2.0, 2.0, 2.0, 2.0]
        end
    end

    for (i, gen) in data_math["gen"]
        gen["pmax"] =  [1e4, 1e4, 1e4]
        gen["pmin"] = -[1e4, 1e4, 1e4]
        gen["qmax"] =  [1e4, 1e4, 1e4]
        gen["qmin"] = -[1e4, 1e4, 1e4]
    end

    for (i, load) in data_math["load"]
        load["pd"] *= load_multiplier
        load["qd"] *= load_multiplier
    end

    return data_math
end

# ═══════════════════════════════════════════════════════════════
# NETWORK SUMMARY
# ═══════════════════════════════════════════════════════════════
function summarise_network(data_math)
    pd_vals = [load["pd"][1] for (i, load) in data_math["load"]]
    qd_vals = [load["qd"][1] for (i, load) in data_math["load"]]
    println("  Network summary (pu, self-consistent):")
    println("    Buses        : $(length(data_math["bus"]))")
    println("    Loads        : $(length(data_math["load"]))")
    println("    pd per load  : $(round(minimum(pd_vals),sigdigits=3)) – $(round(maximum(pd_vals),sigdigits=3)) pu")
    println("    qd per load  : $(round(minimum(qd_vals),sigdigits=3)) – $(round(maximum(qd_vals),sigdigits=3)) pu")
    println("    Total pd     : $(round(sum(pd_vals),sigdigits=4)) pu")
    println("    Total qd     : $(round(sum(qd_vals),sigdigits=4)) pu")
end

# ═══════════════════════════════════════════════════════════════
# PV PLACEMENT
# curtailable=true  → pmin=0,    solver can reduce output (hides stress)
# curtailable=false → pmin=pmax, fixed injection (forces STATCOM to act)
# ═══════════════════════════════════════════════════════════════
function add_pv!(data_math;
        pv_scale    = 3.0,
        q_scale     = 1.25,
        spacing     = 1,
        curtailable = true)

    source_buses = Set([i for (i, bus) in data_math["bus"] if bus["bus_type"] == 3])
    load_ids = sort(
        collect(keys(data_math["load"])),
        by = x -> tryparse(Int, x) === nothing ? 0 : parse(Int, x)
    )

    pv_ids = String[]
    for i in 1:spacing:length(load_ids)
        load       = data_math["load"][load_ids[i]]
        target_bus = string(load["load_bus"])
        target_bus ∈ source_buses && continue

        pd   = load["pd"][1]
        pmax = pv_scale * pd
        smax = q_scale  * pmax
        qlim = sqrt(max(0.0, smax^2 - pmax^2))

        gen_id = string(length(data_math["gen"]) + 1)
        data_math["gen"][gen_id] = deepcopy(data_math["gen"]["1"])
        gen = data_math["gen"][gen_id]

        gen["gen_bus"] = load["load_bus"]
        gen["type"]    = "PV"
        gen["name"]    = "pv_load_$(load_ids[i])"

        if curtailable
            # Solver free to reduce output — curtailment competes with STATCOM
            gen["pmax"] =  pmax * ones(3)
            gen["pmin"] =  zeros(3)
            gen["qmax"] =  smax * ones(3)
            gen["qmin"] = -smax * ones(3)
            gen["cost"] = [-1.0 0.0]   # incentivise generation
        else
            # Fixed injection — STATCOM must handle voltage, not curtailment
            gen["pmax"] =  pmax * ones(3)
            gen["pmin"] =  pmax * ones(3)   # pmin=pmax: no curtailment allowed
            gen["qmax"] =  qlim * ones(3)
            gen["qmin"] = -qlim * ones(3)
            gen["cost"] = [0.0 0.0]    # no cost — output is mandatory
        end

        push!(pv_ids, gen_id)
    end

    mode = curtailable ? "curtailable" : "FIXED (non-curtailable)"
    println("  PV: $(length(pv_ids)) units  pv_scale=$(pv_scale)×pd  q_scale=$(q_scale)  spacing=$(spacing)  mode=$(mode)")
    return pv_ids
end

# ═══════════════════════════════════════════════════════════════
# STATCOM PLACEMENT
# ═══════════════════════════════════════════════════════════════
function add_statcoms!(data_math;
        q_scale = 1.0,
        spacing = 5)

    load_ids = sort(
        collect(keys(data_math["load"])),
        by = x -> tryparse(Int, x) === nothing ? 0 : parse(Int, x)
    )

    statcom_ids = String[]
    for i in 1:spacing:length(load_ids)
        load = data_math["load"][load_ids[i]]
        qd   = load["qd"][1]
        qlim = q_scale * qd

        gen_id = string(length(data_math["gen"]) + 1)
        data_math["gen"][gen_id] = deepcopy(data_math["gen"]["1"])
        gen = data_math["gen"][gen_id]

        gen["gen_bus"] = load["load_bus"]
        gen["type"]    = "STATCOM"
        gen["name"]    = "statcom_load_$(load_ids[i])"

        gen["pmax"] =  zeros(3)
        gen["pmin"] =  zeros(3)
        gen["qmax"] =  qlim * ones(3)
        gen["qmin"] = -qlim * ones(3)
        gen["cost"] = [0.01 0.0]

        push!(statcom_ids, gen_id)
    end

    println("  STATCOM: $(length(statcom_ids)) units  q_scale=$(q_scale)×qd  spacing=$(spacing)")
    return statcom_ids
end

# ═══════════════════════════════════════════════════════════════
# SOLVE AND REPORT
# ═══════════════════════════════════════════════════════════════
function solve_and_report(data_math, label)
    PMD.add_start_vrvi!(data_math)
    model = PMD.instantiate_mc_model(
        data_math,
        PMD.IVRENPowerModel,
        PMD.build_mc_opf
    )
    result = PMD.optimize_model!(model, optimizer=ipopt_solver)
    status = result["termination_status"]
    println("\n  [$label]  →  $(status)")

    if string(status) in ["LOCALLY_SOLVED", "OPTIMAL"]
        sol = result["solution"]

        # Voltage profile
        vm_all     = Float64[]
        bus_labels = String[]
        for (b, bus) in data_math["bus"]
            bus["bus_type"] == 3   && continue
            !haskey(sol["bus"], b) && continue
            for p in 1:3
                vr = get(sol["bus"][b], "vr", zeros(4))[p]
                vi = get(sol["bus"][b], "vi", zeros(4))[p]
                push!(vm_all,     abs(vr + im * vi))
                push!(bus_labels, b)
            end
        end

        if !isempty(vm_all)
            v_min   = minimum(vm_all)
            v_mean  = mean(vm_all)
            v_max   = maximum(vm_all)
            n_over  = count(v -> v > 1.10, vm_all)
            n_under = count(v -> v < 0.90, vm_all)

            println("    Voltage min/mean/max  : $(round(v_min,digits=4)) / $(round(v_mean,digits=4)) / $(round(v_max,digits=4)) pu")
            println("    Violations  > 1.10 pu : $n_over phases")
            println("    Violations  < 0.90 pu : $n_under phases")

            if n_over > 0
                worst = sortperm(vm_all, rev=true)[1:min(5, n_over)]
                println("    Worst overvoltage:")
                for idx in worst
                    vm_all[idx] > 1.10 &&
                    println("      bus $(bus_labels[idx]) → $(round(vm_all[idx],digits=4)) pu")
                end
            end
            if n_under > 0
                worst = sortperm(vm_all)[1:min(5, n_under)]
                println("    Worst undervoltage:")
                for idx in worst
                    vm_all[idx] < 0.90 &&
                    println("      bus $(bus_labels[idx]) → $(round(vm_all[idx],digits=4)) pu")
                end
            end
        end

        # PV dispatch
        pv_gens = [(i, g) for (i, g) in data_math["gen"] if get(g, "type", "") == "PV"]
        if !isempty(pv_gens) && haskey(sol, "gen")
            pg_vals  = [sum(get(sol["gen"][i], "pg", zeros(3))) for (i, g) in pv_gens if haskey(sol["gen"], i)]
            pg_rated = [sum(g["pmax"])                          for (i, g) in pv_gens]
            util     = sum(pg_vals) / max(1e-9, sum(pg_rated)) * 100
            println("    PV utilisation        : $(round(util,digits=1))%  ($(round(sum(pg_vals),sigdigits=4)) / $(round(sum(pg_rated),sigdigits=4)) pu)")
        end

        # STATCOM dispatch
        st_gens = [(i, g) for (i, g) in data_math["gen"] if get(g, "type", "") == "STATCOM"]
        if !isempty(st_gens) && haskey(sol, "gen")
            qg_vals  = [sum(get(sol["gen"][i], "qg", zeros(3))) for (i, g) in st_gens if haskey(sol["gen"], i)]
            qg_rated = [sum(g["qmax"])                          for (i, g) in st_gens]
            util     = abs(sum(qg_vals)) / max(1e-9, sum(qg_rated)) * 100
            direction = sum(qg_vals) >= 0 ? "injecting ↑V" : "absorbing ↓V"
            println("    STATCOM utilisation   : $(round(util,digits=1))%  $(direction)  ($(round(sum(qg_vals),sigdigits=4)) / $(round(sum(qg_rated),sigdigits=4)) pu)")
        end

    else
        println("    WARNING: $(status)")
        println("    Objective : $(get(result, "objective", "N/A"))")
        println("    Tip: check if PV injection is too large for network to absorb")
    end

    return result
end

# ═══════════════════════════════════════════════════════════════
# MAIN — flip booleans to run individual cases
# ═══════════════════════════════════════════════════════════════
run_case1 = true
run_case2 = true
run_case3 = true

# Stress scale identified from Case 2 sweep:
# pv_scale=3.0 → first violation (1716 phases > 1.10 pu, max 1.127 pu)
const STRESS_SCALE = 3.0

dm_ref = load_base_network(data_path)
summarise_network(dm_ref)

# ───────────────────────────────────────────────────────────────
## Case 1: Natural baseline — no bounds, no DER
# ───────────────────────────────────────────────────────────────
if run_case1
    println("\n" * "="^55)
    println(" CASE 1: Natural Baseline (no bounds, no DER)")
    println("="^55)
    dm1 = load_base_network(data_path, enforce_bounds=false)
    r1  = solve_and_report(dm1, "Baseline")
end

# ───────────────────────────────────────────────────────────────
## Case 2: PV stress — no bounds, curtailable
# Shows natural voltage rise as PV penetration increases
# Key result: pv_scale=3.0 is the first violation point
# ───────────────────────────────────────────────────────────────
if run_case2
    println("\n" * "="^55)
    println(" CASE 2: PV Penetration Sweep (no bounds, curtailable)")
    println(" Identifies stress threshold")
    println("="^55)

    for scale in [1.0, 2.0, 3.0]   # stop at known stress point
        println("\n  ── pv_scale = $(scale)×pd ──")
        dm = load_base_network(data_path, enforce_bounds=false)
        add_pv!(dm; pv_scale=scale, q_scale=1.25, spacing=1, curtailable=true)
        solve_and_report(dm, "PV scale=$(scale)×pd")
    end
end

# ───────────────────────────────────────────────────────────────
## Case 3: STATCOM correction — bounds enforced, PV NON-CURTAILABLE
#
# Critical difference from previous runs:
#   curtailable=false → pmin=pmax → PV cannot reduce output
#   The solver MUST use STATCOM reactive absorption to stay within bounds
#   Previously curtailment was doing the work, making STATCOMs idle
#
# Structure:
#   3a — PV fixed, no STATCOM, bounds enforced → expect infeasible/violations
#   3b — PV fixed + STATCOM sweep → find minimum viable STATCOM config
# ───────────────────────────────────────────────────────────────
if run_case3
    println("\n" * "="^55)
    println(" CASE 3: STATCOM Correction")
    println(" PV: FIXED at $(STRESS_SCALE)×pd (non-curtailable)")
    println(" Bounds: enforced (0.90–1.10 pu)")
    println("="^55)

    # 3a: Fixed PV, no STATCOM — establishes whether problem is infeasible
    println("\n  ── 3a: Fixed PV only, bounds enforced ──")
    println("        Expect: infeasible or violations — STATCOM needed")
    dm3a = load_base_network(data_path, enforce_bounds=true)
    add_pv!(dm3a; pv_scale=STRESS_SCALE, q_scale=1.25, spacing=1, curtailable=false)
    r3a = solve_and_report(dm3a, "Fixed PV only, no STATCOM")

    # 3b: Fixed PV + STATCOM sweep
    println("\n  ── 3b: Fixed PV + STATCOM sweep ──")
    configs = [
        (1.0, 10, "baseline      1×qd  sp=10  →  6 STATCOMs"),
        (2.0, 10, "double Q      2×qd  sp=10  →  6 STATCOMs"),
        (4.0, 10, "quad Q        4×qd  sp=10  →  6 STATCOMs"),
        (1.0,  5, "denser        1×qd  sp=5   → 11 STATCOMs"),
        (2.0,  5, "double both   2×qd  sp=5   → 11 STATCOMs"),
        (4.0,  5, "quad+dense    4×qd  sp=5   → 11 STATCOMs"),
        (2.0,  2, "max dense     2×qd  sp=2   → 26 STATCOMs"),
        (1.0,  1, "full cover    1×qd  sp=1   → 52 STATCOMs"),
        (2.0,  1, "full+double   2×qd  sp=1   → 52 STATCOMs"),
    ]

    for (q_sc, sp, lbl) in configs
        println("\n  ── $lbl ──")
        local dm3 = load_base_network(data_path, enforce_bounds=true)
        add_pv!(dm3;       pv_scale=STRESS_SCALE, q_scale=1.25, spacing=1, curtailable=false)
        add_statcoms!(dm3; q_scale=q_sc, spacing=sp)
        solve_and_report(dm3, "PV+STATCOM q=$(q_sc)×qd sp=$(sp)")
    end
end=#

#=
# ==============================================================================
# SCRIPT: STATCOM & PV Hosting Capacity Assessment
# PROJECT: Undergraduate Thesis - Power Systems Engineering
# AUTHOR: Pat
# DATE: May 2026
# ==============================================================================
#
# DESCRIPTION:
# This script performs a multi-stage Optimal Power Flow (OPF) analysis on the 
# ENWL (Electricity North West Limited) 4-wire distribution network. It models 
# the interaction between high-penetration Solar PV and distributed STATCOMs.
#
# CORE FUNCTIONALITY:
# 1. BASELINE ANALYSIS: Establishes natural voltage profiles under stress 
#    without DER or regulatory constraints (Case 1 & 2).
# 2. HOSTING CAPACITY: Quantifies PV curtailment required to maintain 
#    statutory voltage limits (0.90 - 1.10 pu) under high-injection scenarios.
# 3. MITIGATION STUDY: Evaluates the effectiveness of distributed STATCOMs 
#    in "relaxing" voltage constraints through reactive power (Q) absorption.
# 4. SENSITIVITY SWEEPS: Analyzes the impact of STATCOM rating (kVAR) and 
#    spatial density (placement spacing) on total network hosting capacity.
#
# TECHNICAL STACK:
# - Framework: PowerModelsDistribution.jl (PMD)
# - Model: IVREN (Current-Voltage Rectangular Form)
# - Solver: Ipopt (Interior Point Optimizer)
# - Data Source: ENWL 4w_Network1_Feeder1 (OpenDSS format)
# ==============================================================================
=#
using Pkg
Pkg.activate("./")
using rosetta_distribution_opf
import PowerModelsDistribution
import InfrastructureModels
using Ipopt
using JuMP
using Statistics

const PMD  = PowerModelsDistribution
const RPMD = rosetta_distribution_opf
const IM   = InfrastructureModels
PMD.silence!()

ipopt_solver = JuMP.optimizer_with_attributes(
    Ipopt.Optimizer,
    "print_level"           => 0,
    "sb"                    => "yes",
    "warm_start_init_point" => "yes",
    "max_iter"              => 50000
)

data_path = "./data/ENWL_4w_Network1_Feeder1/Master.dss"

# ═══════════════════════════════════════════════════════════════
# NETWORK LOADER
# ═══════════════════════════════════════════════════════════════
function load_base_network(data_path;
        load_multiplier = 1.0,
        enforce_bounds  = false)

    data_eng  = PMD.parse_file(data_path, transformations=[PMD.transform_loops!])
    data_math = PMD.transform_data_model(
        data_eng, multinetwork=false, kron_reduce=false, phase_project=false
    )

    for (i, bus) in data_math["bus"]
        if bus["bus_type"] == 3
            bus["vmin"] = zeros(4)
            bus["vmax"] = 10.0 * ones(4)
        elseif enforce_bounds
            bus["vmin"] = [0.90, 0.90, 0.90, 0.00]
            bus["vmax"] = [1.10, 1.10, 1.10, 0.20]
        else
            bus["vmin"] = [0.50, 0.50, 0.50, 0.00]
            bus["vmax"] = [1.50, 1.50, 1.50, 1.50]
        end
    end

    # Unconstrain slack generator
    for (i, gen) in data_math["gen"]
        gen["pmax"] =  [1e4, 1e4, 1e4]
        gen["pmin"] = -[1e4, 1e4, 1e4]
        gen["qmax"] =  [1e4, 1e4, 1e4]
        gen["qmin"] = -[1e4, 1e4, 1e4]
    end

    for (i, load) in data_math["load"]
        load["pd"] *= load_multiplier
        load["qd"] *= load_multiplier
    end

    return data_math
end

# ═══════════════════════════════════════════════════════════════
# NETWORK SUMMARY
# ═══════════════════════════════════════════════════════════════
function summarise_network(data_math)
    pd_vals = [load["pd"][1] for (i, load) in data_math["load"]]
    qd_vals = [load["qd"][1] for (i, load) in data_math["load"]]
    println("  Network summary (pu, self-consistent):")
    println("    Buses        : $(length(data_math["bus"]))")
    println("    Loads        : $(length(data_math["load"]))")
    println("    pd per load  : $(round(minimum(pd_vals),sigdigits=3)) – $(round(maximum(pd_vals),sigdigits=3)) pu")
    println("    qd per load  : $(round(minimum(qd_vals),sigdigits=3)) – $(round(maximum(qd_vals),sigdigits=3)) pu")
    println("    Total pd     : $(round(sum(pd_vals),sigdigits=4)) pu")
    println("    Total qd     : $(round(sum(qd_vals),sigdigits=4)) pu")
end

# ═══════════════════════════════════════════════════════════════
# PV PLACEMENT
# pv_cost < 0 → maximise generation (hosting capacity objective)
# pv_cost > 0 → penalise generation (wrong for this study)
# ═══════════════════════════════════════════════════════════════
function add_pv!(data_math;
        pv_scale = 3.0,
        q_scale  = 1.25,
        spacing  = 1,
        pv_cost  = -1000.0)

    source_buses = Set([i for (i, bus) in data_math["bus"] if bus["bus_type"] == 3])
    load_ids = sort(
        collect(keys(data_math["load"])),
        by = x -> tryparse(Int, x) === nothing ? 0 : parse(Int, x)
    )

    pv_ids = String[]
    for i in 1:spacing:length(load_ids)
        load       = data_math["load"][load_ids[i]]
        target_bus = string(load["load_bus"])
        target_bus ∈ source_buses && continue

        pd   = load["pd"][1]
        pmax = pv_scale * pd
        smax = q_scale  * pmax
        qlim = sqrt(max(0.0, smax^2 - pmax^2))

        gen_id = string(length(data_math["gen"]) + 1)
        data_math["gen"][gen_id] = deepcopy(data_math["gen"]["1"])
        gen = data_math["gen"][gen_id]

        gen["gen_bus"] = load["load_bus"]
        gen["type"]    = "PV"
        gen["name"]    = "pv_load_$(load_ids[i])"

        gen["pmax"] =  pmax * ones(3)
        gen["pmin"] =  zeros(3)
        gen["qmax"] =  smax * ones(3)
        gen["qmin"] = -smax * ones(3)
        gen["cost"] = [pv_cost 0.0]

        push!(pv_ids, gen_id)
    end

    println("  PV: $(length(pv_ids)) units  pv_scale=$(pv_scale)×pd  q_scale=$(q_scale)  spacing=$(spacing)  cost=$(pv_cost)")
    return pv_ids
end

# ═══════════════════════════════════════════════════════════════
# STATCOM PLACEMENT
# statcom_cost > 0, small → deployed only when needed by solver
# ═══════════════════════════════════════════════════════════════
function add_statcoms!(data_math;
        q_scale      = 1.0,
        spacing      = 1,
        statcom_cost = 1.0)

    load_ids = sort(
        collect(keys(data_math["load"])),
        by = x -> tryparse(Int, x) === nothing ? 0 : parse(Int, x)
    )

    statcom_ids = String[]
    for i in 1:spacing:length(load_ids)
        load = data_math["load"][load_ids[i]]
        qd   = load["qd"][1]
        qlim = q_scale * qd

        gen_id = string(length(data_math["gen"]) + 1)
        data_math["gen"][gen_id] = deepcopy(data_math["gen"]["1"])
        gen = data_math["gen"][gen_id]

        gen["gen_bus"] = load["load_bus"]
        gen["type"]    = "STATCOM"
        gen["name"]    = "statcom_load_$(load_ids[i])"

        gen["pmax"] =  zeros(3)
        gen["pmin"] =  zeros(3)
        gen["qmax"] =  qlim * ones(3)
        gen["qmin"] = -qlim * ones(3)
        gen["cost"] = [statcom_cost 0.0]

        push!(statcom_ids, gen_id)
    end

    println("  STATCOM: $(length(statcom_ids)) units  q_scale=$(q_scale)×qd  spacing=$(spacing)  cost=$(statcom_cost)")
    return statcom_ids
end

# ═══════════════════════════════════════════════════════════════
# SOLVE AND REPORT
# Uses PMD.build_mc_opf with gen["cost"] coefficients
# This correctly wires gen["cost"] coefficients through the
# rosetta OPF builder rather than PMD's default build_mc_opf
# ═══════════════════════════════════════════════════════════════
function solve_and_report(data_math, label)
    PMD.add_start_vrvi!(data_math)

    # PMD standard builder — correctly uses gen["cost"] coefficients
    # RPMD.build_mc_opf_mx requires add_inverter_losses! pre-processing
    # which is not applicable for generic PV/STATCOM devices
    model = PMD.instantiate_mc_model(
        data_math,
        PMD.IVRENPowerModel,
        PMD.build_mc_opf
    )

    result = PMD.optimize_model!(model, optimizer=ipopt_solver)
    status = result["termination_status"]
    println("\n  [$label]  →  $(status)")

    if string(status) in ["LOCALLY_SOLVED", "OPTIMAL"]
        sol = result["solution"]

        # Voltage profile: phases 1–3, skip slack
        vm_all     = Float64[]
        bus_labels = String[]
        for (b, bus) in data_math["bus"]
            bus["bus_type"] == 3   && continue
            !haskey(sol["bus"], b) && continue
            for p in 1:3
                vr = get(sol["bus"][b], "vr", zeros(4))[p]
                vi = get(sol["bus"][b], "vi", zeros(4))[p]
                push!(vm_all,     abs(vr + im * vi))
                push!(bus_labels, b)
            end
        end

        if !isempty(vm_all)
            v_min   = minimum(vm_all)
            v_mean  = mean(vm_all)
            v_max   = maximum(vm_all)
            n_over  = count(v -> v > 1.10, vm_all)
            n_under = count(v -> v < 0.90, vm_all)

            println("    Voltage min/mean/max  : $(round(v_min,digits=4)) / $(round(v_mean,digits=4)) / $(round(v_max,digits=4)) pu")
            println("    Violations  > 1.10 pu : $n_over phases")
            println("    Violations  < 0.90 pu : $n_under phases")

            if n_over > 0
                worst = sortperm(vm_all, rev=true)[1:min(5, n_over)]
                println("    Worst overvoltage:")
                for idx in worst
                    vm_all[idx] > 1.10 &&
                    println("      bus $(bus_labels[idx]) → $(round(vm_all[idx],digits=4)) pu")
                end
            end
            if n_under > 0
                worst = sortperm(vm_all)[1:min(5, n_under)]
                println("    Worst undervoltage:")
                for idx in worst
                    vm_all[idx] < 0.90 &&
                    println("      bus $(bus_labels[idx]) → $(round(vm_all[idx],digits=4)) pu")
                end
            end
        end

        # PV dispatch
        pv_gens = [(i, g) for (i, g) in data_math["gen"] if get(g, "type", "") == "PV"]
        if !isempty(pv_gens) && haskey(sol, "gen")
            pg_vals  = [sum(get(sol["gen"][i], "pg", zeros(3))) for (i, g) in pv_gens if haskey(sol["gen"], i)]
            pg_rated = [sum(g["pmax"])                          for (i, g) in pv_gens]
            pg_total = sum(pg_vals)
            pg_cap   = sum(pg_rated)
            util     = pg_total / max(1e-9, pg_cap) * 100
            curtail  = 100.0 - util
            println("    PV output / capacity  : $(round(pg_total,sigdigits=4)) / $(round(pg_cap,sigdigits=4)) pu")
            println("    PV utilisation        : $(round(util,digits=1))%  →  curtailment: $(round(max(0,curtail),digits=1))%")
        end

        # STATCOM dispatch
        st_gens = [(i, g) for (i, g) in data_math["gen"] if get(g, "type", "") == "STATCOM"]
        if !isempty(st_gens) && haskey(sol, "gen")
            qg_vals  = [sum(get(sol["gen"][i], "qg", zeros(3))) for (i, g) in st_gens if haskey(sol["gen"], i)]
            qg_rated = [sum(g["qmax"])                          for (i, g) in st_gens]
            qg_total = sum(qg_vals)
            qg_cap   = sum(qg_rated)
            util     = abs(qg_total) / max(1e-9, qg_cap) * 100
            direction = qg_total >= 0 ? "injecting ↑V" : "absorbing ↓V"
            println("    STATCOM Q / capacity  : $(round(qg_total,sigdigits=4)) / $(round(qg_cap,sigdigits=4)) pu")
            println("    STATCOM utilisation   : $(round(util,digits=1))%  $(direction)")
        end

    else
        println("    WARNING: $(status)")
        println("    Objective : $(get(result, "objective", "N/A"))")
    end

    return result
end

# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════
do_case1 = true
do_case2 = true
do_case3 = true

# From Case 2: pv_scale=3.0 → first violations at 1.1266 pu (no bounds)
# For Case 3 we need a scale where without STATCOM the solver must curtail
# pv_scale=3.0 already solves at 100% with bounds (voltage just hits 1.10 pu)
# so we push to 4.0 and 5.0 to force genuine curtailment
const STRESS_SCALE   = 5.0    # genuinely requires STATCOM to achieve full output
const PV_COST        = -1000.0
const STATCOM_COST   =  1.0

dm_ref = load_base_network(data_path)
summarise_network(dm_ref)

# ───────────────────────────────────────────────────────────────
## Case 1: Natural baseline — no bounds, no DER
# ───────────────────────────────────────────────────────────────
if do_case1
    println("\n" * "="^55)
    println(" CASE 1: Natural Baseline (no bounds, no DER)")
    println("="^55)
    dm1 = load_base_network(data_path, enforce_bounds=false)
    r1  = solve_and_report(dm1, "Baseline")
end

# ───────────────────────────────────────────────────────────────
## Case 2: PV penetration sweep — no bounds
# Establishes natural voltage rise without enforcement
# ───────────────────────────────────────────────────────────────
if do_case2
    println("\n" * "="^55)
    println(" CASE 2: PV Penetration Sweep (no bounds)")
    println(" Natural voltage profile as PV penetration increases")
    println("="^55)

    for scale in [1.0, 2.0, 3.0, 4.0, 5.0]
        println("\n  ── pv_scale = $(scale)×pd ──")
        dm = load_base_network(data_path, enforce_bounds=false)
        add_pv!(dm; pv_scale=scale, q_scale=1.25, spacing=1, pv_cost=-1000.0)
        solve_and_report(dm, "PV scale=$(scale)×pd  [no bounds]")
    end
end

# ───────────────────────────────────────────────────────────────
## Case 3: STATCOM hosting capacity study — bounds enforced
#
# Objective: MAXIMISE PV output (cost = -1000 per unit P)
# STATCOM: small positive cost (1.0) — deployed only when needed
#
# At pv_scale=5.0:
#   Without STATCOM: solver must curtail PV to satisfy 1.10 pu bound
#   With STATCOM:    Q absorption relaxes the voltage constraint
#                    → solver can host more PV before hitting the bound
#
# Key metric: PV utilisation %
#   Rises as STATCOM capacity increases = more PV hosted
#   Difference vs 3a baseline = STATCOM's hosting capacity contribution
# ───────────────────────────────────────────────────────────────
if do_case3
    println("\n" * "="^55)
    println(" CASE 3: STATCOM Hosting Capacity Study")
    println(" PV stress: $(STRESS_SCALE)×pd  |  bounds: 0.90–1.10 pu")
    println(" PV cost: $(PV_COST)  |  STATCOM cost: $(STATCOM_COST)")
    println(" Builder: PMD.build_mc_opf  |  Objective: cost (gen[\"cost\"])")
    println("="^55)

    # 3a: PV only — establishes hosting capacity WITHOUT STATCOM
    println("\n  ── 3a: PV only — hosting capacity baseline ──")
    dm3a = load_base_network(data_path, enforce_bounds=true)
    add_pv!(dm3a; pv_scale=STRESS_SCALE, q_scale=1.25, spacing=1, pv_cost=PV_COST)
    r3a = solve_and_report(dm3a, "PV only  [no STATCOM]")

    # 3b: STATCOM rating sweep at full density
    # All 52 load buses get a STATCOM, rating increases
    # Shows: at what rating does PV utilisation improve vs 3a?
    println("\n  ── 3b: STATCOM rating sweep (52 STATCOMs, spacing=1) ──")

    for q_sc in [1, 2, 5, 10, 20, 50, 100, 200, 500]
        println("\n  ── q_scale = $(q_sc)×qd  (52 STATCOMs) ──")
        local dm3b = load_base_network(data_path, enforce_bounds=true)
        add_pv!(dm3b;       pv_scale=STRESS_SCALE, q_scale=1.25, spacing=1, pv_cost=PV_COST)
        add_statcoms!(dm3b; q_scale=q_sc,          spacing=1,    statcom_cost=STATCOM_COST)
        solve_and_report(dm3b, "PV+STATCOM q=$(q_sc)×qd  52 units")
    end

    # 3c: Density sweep at rating identified in 3b
    # Fix rating, reduce number of STATCOMs — find minimum viable deployment
    # Update min_viable_q after seeing 3b results
    min_viable_q = 100

    println("\n  ── 3c: Density sweep at q_scale=$(min_viable_q)×qd ──")
    println("        How few STATCOMs achieve the same hosting benefit?")

    for sp in [1, 2, 4, 6, 8, 10]
        n_st = ceil(Int, 52 / sp)
        println("\n  ── spacing=$(sp) → ~$(n_st) STATCOMs ──")
        local dm3c = load_base_network(data_path, enforce_bounds=true)
        add_pv!(dm3c;       pv_scale=STRESS_SCALE, q_scale=1.25,   spacing=1, pv_cost=PV_COST)
        add_statcoms!(dm3c; q_scale=min_viable_q,  spacing=sp,     statcom_cost=STATCOM_COST)
        solve_and_report(dm3c, "PV+STATCOM q=$(min_viable_q)×qd  $(n_st) units")
    end
end