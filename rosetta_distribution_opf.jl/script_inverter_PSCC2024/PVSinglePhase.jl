#=
# ==============================================================================
# SCRIPT: STATCOM & PV Hosting Capacity Assessment — Phase-Realistic PV
# PROJECT: Undergraduate Thesis - Power Systems Engineering
# AUTHOR: Pat
# DATE: Aug 2026
# ==============================================================================
#
# DESCRIPTION:
# This version fixes the "PV as balanced 3-phase injection" issue identified
# in seminarScript.jl. Two structural changes from that script:
#
#  1. PV IS NOW SINGLE-PHASE, MATCHED TO ITS COLOCATED LOAD'S PHASE.
#     Residential rooftop PV is a customer-side asset wired to whichever one
#     phase the household's supply is on — it has no physical access to the
#     other two conductors. add_pv! now reads data_math["load"][id]["connections"]
#     for the colocated load and copies those connections onto the PV gen, so a
#     load on phase B gets a phase-B-only PV unit, not a balanced 3-phase one.
#
#  2. STATCOM STAYS 3-PHASE (unchanged in principle from before), but is now
#     sized directly in kVAr instead of a q_scale×qd multiplier, and that
#     kVAr figure is treated as the TOTAL NAMEPLATE rating of the unit, split
#     evenly across its 3 legs (qmax_per_leg = rating_kvar/3). This matches
#     how real 4-leg D-STATCOM converters are rated (whole-unit apparent
#     power, not independently-rated legs) and removes the "is this number
#     per-unit or total" ambiguity that q_scale×qd created. STATCOM stays
#     3-phase because it is a network-side asset that sees all 3 phases +
#     neutral at every bus regardless of which phase the local load taps,
#     and per-phase-independent Q dispatch is the entire mechanism by which
#     it corrects imbalance (see Case 3 per-phase Q breakdown).
#
# ALSO CHANGED:
#  - All PV/STATCOM sizing is now in real kW / kVAr, not pu multipliers of
#    load demand. (Loads_v2.txt is flat 1kW everywhere, so pv_scale×pd was
#    never really "scaling to demand" — it was scaling to a placeholder.)
#  - STATCOM per-phase Q is now also printed in kVAr, not just pu.
#  - STATCOM utilisation is now reported TWO ways:
#      "net"   = |sum of signed Q across all units & phases| / total capacity
#                (can be misleadingly small near a sign-flip in dispatch)
#      "gross" = sum of |Q| across all units & phases / total capacity
#                (reflects how hard the fleet is actually working)
#  - Cases 4-6 (HC curve, targeted placement, decomposition) are removed for
#    this iteration — this script is deliberately scoped to Cases 1-3 so the
#    effect of phase-matched PV can be seen in isolation before revisiting
#    the more elaborate studies.
#
# TECHNICAL STACK: unchanged
#   Framework : PowerModelsDistribution.jl (PMD)
#   Model     : IVREN (Current-Voltage Rectangular Form)
#   Solver    : Ipopt
#   Data      : ENWL 4w_Network1_Feeder1 (OpenDSS format)
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
using CairoMakie
using Printf

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

data_path = "./rosetta_distribution_opf.jl/data/ENWL_4w_Network1_Feeder1/Master.dss"

# ═══════════════════════════════════════════════════════════════
# NETWORK LOADER  — unchanged from seminarScript.jl
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
# NETWORK SUMMARY  — unchanged
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
# ADD: load phase distribution — shows how unbalanced the feeder's
# customer allocation already is, before any PV is added. This is
# the distribution that phase-matched PV will inherit.
# ═══════════════════════════════════════════════════════════════
function report_load_phase_distribution(data_math)
    tally = Dict(1 => 0, 2 => 0, 3 => 0)
    threephase = 0

    for (id, load) in data_math["load"]
        phase_conns = filter(c -> c != 4, load["connections"])
        if length(phase_conns) == 1
            tally[phase_conns[1]] += 1
        else
            threephase += 1
        end
    end

    println("  Load phase distribution:")
    println("    Phase A: $(tally[1])   Phase B: $(tally[2])   Phase C: $(tally[3])   Three-phase: $threephase")
    return tally
end

# ═══════════════════════════════════════════════════════════════
# ADD: base-value extraction — unchanged from seminarScript.jl
# ═══════════════════════════════════════════════════════════════
function get_base_values(data_eng, data_math)
    sbase_kva = data_math["settings"]["sbase"] * data_math["settings"]["power_scale_factor"] / 1000
    # ^ verify against your PMD version: power_scale_factor is usually 1000 (W<->kW)

    vbases, _ = PMD.calc_voltage_bases(data_eng, data_eng["settings"]["vbases_default"])

    busid2ebus = Dict(
        string(bus["bus_i"]) => split(bus["source_id"], ".")[end]
        for (i, bus) in data_math["bus"]
    )

    vbase_kv = Dict{String,Float64}()
    for (mbus, ebus) in busid2ebus
        if haskey(vbases, ebus)
            vbase_kv[mbus] = vbases[ebus]
        end
    end

    return sbase_kva, vbase_kv
end
pu_to_kw(p_pu, sbase_kva)   = p_pu * sbase_kva
pu_to_kvar(q_pu, sbase_kva) = q_pu * sbase_kva
pu_to_kv(v_pu, vbase_kv)    = v_pu * vbase_kv
kw_to_pu(p_kw, sbase_kva)   = p_kw / sbase_kva
kvar_to_pu(q_kvar, sbase_kva) = q_kvar / sbase_kva

# ═══════════════════════════════════════════════════════════════
# PV PLACEMENT — REWRITTEN
#   • Sized directly in kW (real units), not pv_scale×pd.
#   • Single-phase: connections copied from the colocated load, so
#     a load on phase B gets a phase-B-only PV unit.
#   • Inverter Q headroom still expressed as q_scale over the kW
#     rating (s_max = q_scale × pv_kw), consistent with smart-
#     inverter sizing convention (e.g. a 5kW/6.25kVA inverter).
# ═══════════════════════════════════════════════════════════════
function add_pv!(data_math;
        pv_kw     = 5.0,      # real active rating per installed unit, kW
        q_scale   = 1.25,     # inverter kVA headroom multiplier (s_max = q_scale*pv_kw)
        spacing   = 1,
        pv_cost   = -1000.0,
        sbase_kva = SBASE_KVA)

    pv_pu = kw_to_pu(pv_kw, sbase_kva)
    s_pu  = q_scale * pv_pu
    q_lim = sqrt(max(0.0, s_pu^2 - pv_pu^2))

    source_buses = Set([i for (i, bus) in data_math["bus"] if bus["bus_type"] == 3])
    load_ids = sort(
        collect(keys(data_math["load"])),
        by = x -> tryparse(Int, x) === nothing ? 0 : parse(Int, x)
    )

    pv_ids  = String[]
    n_ph1   = 0   # single-phase units placed
    n_ph3   = 0   # three-phase units placed (rare/none for this feeder)
    n_skip  = 0   # loads with no usable phase connection

    for i in 1:spacing:length(load_ids)
        load       = data_math["load"][load_ids[i]]
        target_bus = string(load["load_bus"])
        target_bus ∈ source_buses && continue

        conns       = load["connections"]
        phase_conns = filter(c -> c != 4, conns)   # actual phase(s), excluding neutral
        n_ph        = length(phase_conns)

        if n_ph == 0
            n_skip += 1
            continue
        end

        gen_id = string(length(data_math["gen"]) + 1)
        data_math["gen"][gen_id] = deepcopy(data_math["gen"]["1"])
        gen = data_math["gen"][gen_id]

        gen["gen_bus"]     = load["load_bus"]
        gen["type"]        = "PV"
        gen["name"]        = "pv_load_$(load_ids[i])"
        gen["connections"] = deepcopy(conns)       # KEY CHANGE: match the load's phase(s)

        gen["pmax"] =  pv_pu * ones(n_ph)
        gen["pmin"] =  zeros(n_ph)
        gen["qmax"] =  q_lim * ones(n_ph)
        gen["qmin"] = -q_lim * ones(n_ph)
        gen["cost"] = [pv_cost 0.0]

        # Defensive: trim any other length-3 template field (e.g. vg/pg/qg
        # starting points) down to n_ph so the gen record stays internally
        # consistent. PMD versions vary in which fields exist here.
        #=
        for (k, v) in gen
            if v isa Vector{<:Real} && length(v) == 3 && n_ph != 3
                gen[k] = v[1:n_ph]
            end
        end
        =#
        for (k, v) in gen
            if v isa Vector{<:Real} && length(v) == 3 && n_ph != 3
                gen[k] = v[phase_conns]
            end
        end

        push!(pv_ids, gen_id)
        n_ph == 1 ? (n_ph1 += 1) : (n_ph3 += 1)
    end

    println("  PV: $(length(pv_ids)) units  pv_kw=$(pv_kw) kW/unit  q_scale=$(q_scale)  spacing=$(spacing)  cost=$(pv_cost)")
    println("      single-phase: $n_ph1   three-phase: $n_ph3   skipped (no phase match): $n_skip")
    return pv_ids
end

# ═══════════════════════════════════════════════════════════════
# STATCOM PLACEMENT — REWRITTEN
#   • Sized directly in kVAr (real units) as a TOTAL nameplate
#     rating per unit, not q_scale×qd.
#   • Stays 3-phase (4-leg converter): total rating split evenly
#     across the 3 legs, qmax_per_leg = rating_kvar/3, so that
#     summing the 3 legs' capacity back up gives exactly the
#     nameplate figure — not 3× it.
#   • target_buses path kept for future targeted-placement studies;
#     unused in Cases 1-3 but left in for continuity.
# ═══════════════════════════════════════════════════════════════
function add_statcoms!(data_math;
        rating_kvar  = 20.0,   # total nameplate rating per unit, kVAr
        spacing      = 1,
        statcom_cost = 1.0,
        target_buses = nothing,
        sbase_kva    = SBASE_KVA)

    q_total_pu = kvar_to_pu(rating_kvar, sbase_kva)
    q_leg_pu   = q_total_pu / 3   # even split across the 3 converter legs

    load_ids = sort(
        collect(keys(data_math["load"])),
        by = x -> tryparse(Int, x) === nothing ? 0 : parse(Int, x)
    )

    statcom_ids = String[]

    if !isnothing(target_buses)
        for bus_id in target_buses
            gen_id = string(length(data_math["gen"]) + 1)
            data_math["gen"][gen_id] = deepcopy(data_math["gen"]["1"])
            gen = data_math["gen"][gen_id]

            gen["gen_bus"] = bus_id
            gen["type"]    = "STATCOM"
            gen["name"]    = "statcom_targeted_bus$(bus_id)"
            gen["pmax"]    =  zeros(3)
            gen["pmin"]    =  zeros(3)
            gen["qmax"]    =  q_leg_pu * ones(3)
            gen["qmin"]    = -q_leg_pu * ones(3)
            gen["cost"]    = [statcom_cost 0.0]

            push!(statcom_ids, gen_id)
        end
        println("  STATCOM: $(length(statcom_ids)) units  TARGETED buses=$(target_buses)  " *
                "rating=$(rating_kvar) kVAr/unit total ($(round(q_leg_pu*sbase_kva,digits=2)) kVAr/leg)  cost=$(statcom_cost)")
        return statcom_ids
    end

    for i in 1:spacing:length(load_ids)
        load = data_math["load"][load_ids[i]]

        gen_id = string(length(data_math["gen"]) + 1)
        data_math["gen"][gen_id] = deepcopy(data_math["gen"]["1"])
        gen = data_math["gen"][gen_id]

        gen["gen_bus"] = load["load_bus"]
        gen["type"]    = "STATCOM"
        gen["name"]    = "statcom_load_$(load_ids[i])"
        gen["pmax"]    =  zeros(3)
        gen["pmin"]    =  zeros(3)
        gen["qmax"]    =  q_leg_pu * ones(3)
        gen["qmin"]    = -q_leg_pu * ones(3)
        gen["cost"]    = [statcom_cost 0.0]

        push!(statcom_ids, gen_id)
    end

    println("  STATCOM: $(length(statcom_ids)) units  rating=$(rating_kvar) kVAr/unit total " *
            "($(round(q_leg_pu*sbase_kva,digits=2)) kVAr/leg)  spacing=$(spacing)  cost=$(statcom_cost)")
    return statcom_ids
end

# ═══════════════════════════════════════════════════════════════
# SOLVE AND REPORT — updated
#   • PV extraction no longer assumes length-3 pg/pmax vectors
#     (single-phase units are length 1).
#   • STATCOM per-phase Q printed in kVAr as well as pu.
#   • STATCOM utilisation reported as both NET and GROSS (see
#     header note for why these can diverge).
# ═══════════════════════════════════════════════════════════════
function solve_and_report(data_math, label; sbase_kva=SBASE_KVA, vbase_kv=VBASE_KV)
    PMD.add_start_vrvi!(data_math)

    model  = PMD.instantiate_mc_model(data_math, PMD.IVRENPowerModel, PMD.build_mc_opf)
    result = PMD.optimize_model!(model, optimizer=ipopt_solver)
    status = result["termination_status"]
    println("\n  [$label]  →  $(status)")

    v_min = NaN;  v_mean = NaN;  v_max = NaN
    n_over = 0;   n_under = 0
    pv_util = NaN;  pv_curtail = NaN
    pv_output = NaN;  pv_capacity = NaN
    st_q_total = NaN;  st_q_cap = NaN
    st_util_net = NaN; st_util_gross = NaN
    vm_per_bus = Dict{String, Float64}()

    if string(status) in ["LOCALLY_SOLVED", "OPTIMAL"]
        sol = result["solution"]

        # ── Voltage profile ─────────────────────────────────────
        vm_all     = Float64[]
        bus_labels = String[]

        for (b, bus) in data_math["bus"]
            bus["bus_type"] == 3   && continue
            !haskey(sol["bus"], b) && continue
            bus_max = 0.0
            for p in 1:3
                vr = get(sol["bus"][b], "vr", zeros(4))[p]
                vi = get(sol["bus"][b], "vi", zeros(4))[p]
                vm = abs(vr + im * vi)
                push!(vm_all,     vm)
                push!(bus_labels, b)
                bus_max = max(bus_max, vm)
            end
            vm_per_bus[b] = bus_max
        end

        if !isempty(vm_all)
            v_min   = minimum(vm_all)
            v_mean  = mean(vm_all)
            v_max   = maximum(vm_all)
            n_over  = count(v -> v > 1.10, vm_all)
            n_under = count(v -> v < 0.90, vm_all)

            println("    Voltage min/mean/max  : $(round(v_min,digits=4)) / $(round(v_mean,digits=4)) / $(round(v_max,digits=4)) pu")
            if !isnothing(vbase_kv)
                v_min_bus = bus_labels[argmin(vm_all)]
                v_max_bus = bus_labels[argmax(vm_all)]
                if haskey(vbase_kv, v_min_bus) && haskey(vbase_kv, v_max_bus)
                    println("    Voltage min/max (kV)  : $(round(pu_to_kv(v_min, vbase_kv[v_min_bus]),digits=3)) / " *
                            "$(round(pu_to_kv(v_max, vbase_kv[v_max_bus]),digits=3)) kV")
                end
            end
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

        # ── PV dispatch ─────────────────────────────────────────
        # NOTE: pv gens are now variable-length (1 for single-phase, 3 for
        # three-phase), so fallback zero-vectors must match each gen's own
        # pmax length rather than assuming 3.
        pv_gens = [(i, g) for (i, g) in data_math["gen"] if get(g, "type", "") == "PV"]
        if !isempty(pv_gens) && haskey(sol, "gen")
            pg_vals = [
                sum(get(sol["gen"][i], "pg", zeros(length(g["pmax"]))))
                for (i, g) in pv_gens if haskey(sol["gen"], i)
            ]
            pg_rated    = [sum(g["pmax"]) for (i, g) in pv_gens]
            pv_output   = sum(pg_vals)
            pv_capacity = sum(pg_rated)
            pv_util     = pv_output / max(1e-9, pv_capacity) * 100
            pv_curtail  = max(0.0, 100.0 - pv_util)
            println("    PV output / capacity  : $(round(pv_output,sigdigits=4)) / $(round(pv_capacity,sigdigits=4)) pu")
            if !isnothing(sbase_kva)
                println("    PV output / capacity  : $(round(pu_to_kw(pv_output,sbase_kva),digits=1)) / " *
                        "$(round(pu_to_kw(pv_capacity,sbase_kva),digits=1)) kW")
            end
            println("    PV utilisation        : $(round(pv_util,digits=1))%  →  curtailment: $(round(pv_curtail,digits=1))%")
        end

        # ── STATCOM dispatch ────────────────────────────────────
        st_gens = [(i, g) for (i, g) in data_math["gen"] if get(g, "type", "") == "STATCOM"]
        if !isempty(st_gens) && haskey(sol, "gen")
            qg_vecs = [get(sol["gen"][i], "qg", zeros(3)) for (i, g) in st_gens if haskey(sol["gen"], i)]

            qg_signed_sum = sum(sum(q) for q in qg_vecs)          # NET: signs can cancel
            qg_abs_sum    = sum(sum(abs.(q)) for q in qg_vecs)    # GROSS: no cancellation
            qg_rated      = [sum(g["qmax"]) for (i, g) in st_gens]

            st_q_total    = qg_signed_sum
            st_q_cap      = sum(qg_rated)
            st_util_net   = abs(st_q_total) / max(1e-9, st_q_cap) * 100
            st_util_gross = qg_abs_sum      / max(1e-9, st_q_cap) * 100
            direction     = st_q_total >= 0 ? "net injecting ↑V" : "net absorbing ↓V"

            println("    STATCOM Q / capacity  : $(round(st_q_total,sigdigits=4)) / $(round(st_q_cap,sigdigits=4)) pu")
            if !isnothing(sbase_kva)
                println("    STATCOM Q / capacity  : $(round(pu_to_kvar(st_q_total,sbase_kva),digits=1)) / " *
                        "$(round(pu_to_kvar(st_q_cap,sbase_kva),digits=1)) kVAr")
            end
            println("    STATCOM utilisation (net)   : $(round(st_util_net,digits=1))%  $(direction)")
            println("    STATCOM utilisation (gross) : $(round(st_util_gross,digits=1))%  " *
                    "(sum of |Q| across units/phases — doesn't let opposing phases cancel out)")

            # Per-phase Q breakdown, now in pu AND kVAr
            phase_labels = ["Phase A", "Phase B", "Phase C"]
            for (pi, ph) in enumerate(phase_labels)
                q_ph_pu = sum(q[pi] for q in qg_vecs)
                arrow   = q_ph_pu >= 0 ? "↑" : "↓"
                if !isnothing(sbase_kva)
                    println("      $(ph) net Q: $(round(q_ph_pu, sigdigits=3)) pu  " *
                            "($(round(pu_to_kvar(q_ph_pu, sbase_kva), digits=2)) kVAr)  $arrow")
                else
                    println("      $(ph) net Q: $(round(q_ph_pu, sigdigits=3)) pu  $arrow")
                end
            end
        end

    else
        println("    WARNING: $(status)")
        println("    Objective : $(get(result, "objective", "N/A"))")
    end

    return (
        status        = string(status),
        v_min         = v_min,
        v_mean        = v_mean,
        v_max         = v_max,
        n_over        = n_over,
        n_under       = n_under,
        pv_util       = pv_util,
        pv_curtail    = pv_curtail,
        pv_output     = pv_output,
        pv_capacity   = pv_capacity,
        st_q_total    = st_q_total,
        st_q_cap      = st_q_cap,
        st_util_net   = st_util_net,
        st_util_gross = st_util_gross,
        vm_per_bus    = vm_per_bus,
        pv_output_kw    = isnothing(sbase_kva) ? NaN : pu_to_kw(pv_output, sbase_kva),
        st_q_total_kvar = isnothing(sbase_kva) ? NaN : pu_to_kvar(st_q_total, sbase_kva),
    )
end

# ═══════════════════════════════════════════════════════════════
# ΔHC SUMMARY TABLE — unchanged
# ═══════════════════════════════════════════════════════════════
function print_delta_hc_table(rows)
    println("\n  ΔHC Summary Table")
    println("  " * "─"^72)
    @printf("  %-32s  %10s  %12s  %14s\n",
            "Configuration", "Util (%)", "Curtail (%)", "ΔHC vs base")
    println("  " * "─"^72)
    baseline = rows[1].util
    for r in rows
        delta     = r.util - baseline
        delta_str = abs(delta) < 0.05 ? "           —" : @sprintf("%+.1f pp", delta)
        @printf("  %-32s  %10s  %12s  %14s\n",
                r.label,
                @sprintf("%.1f%%", r.util),
                @sprintf("%.1f%%", r.curtail),
                delta_str)
    end
    println("  " * "─"^72)
end

# ═══════════════════════════════════════════════════════════════
# PLOTS — trimmed to what Cases 1-3 use
# ═══════════════════════════════════════════════════════════════
mkpath("./plots")

# ── Plot 1: Voltage profile ──────────────────────────────────────────────────
function plot_voltage_profile(r_no_st, r_with_st, pv_kw, best_statcom_kvar;
        savepath="./plots/voltage_profile.pdf")

    common_buses = intersect(keys(r_no_st.vm_per_bus), keys(r_with_st.vm_per_bus))
    sorted_buses = sort(collect(common_buses), by = b -> r_no_st.vm_per_bus[b])
    v_no = [r_no_st.vm_per_bus[b]   for b in sorted_buses]
    v_st = [r_with_st.vm_per_bus[b] for b in sorted_buses]
    xs   = 1:length(sorted_buses)

    fig = Figure(size=(900, 480), backgroundcolor=:white)
    ax  = Axis(fig[1,1],
        xlabel       = "Bus index (sorted by ascending voltage, no STATCOM)",
        ylabel       = "Voltage magnitude (pu)",
        title        = "Voltage profile — single-phase PV, $(pv_kw) kW/unit, with and without STATCOM",
        titlesize    = 15, xlabelsize = 12, ylabelsize = 12,
        ygridvisible = true, xgridvisible = false,
        ygridcolor   = (:black, 0.08),
    )

    hspan!(ax, 1.10, 1.25, color=(:red,  0.07))
    hspan!(ax, 0.75, 0.90, color=(:blue, 0.07))

    lines!(ax, xs, v_no, color=:steelblue,  linewidth=2.5, label="No STATCOM")
    lines!(ax, xs, v_st, color=:darkorange, linewidth=2.5,
           label="With STATCOM ($(best_statcom_kvar) kVAr/unit, 52 units)")

    hlines!(ax, [1.10], color=:red,  linestyle=:dash, linewidth=2.0,
            label="1.10 pu statutory limit")
    hlines!(ax, [0.90], color=:blue, linestyle=:dash, linewidth=2.0,
            label="0.90 pu statutory limit")
    hlines!(ax, [1.00], color=(:black, 0.2), linestyle=:dot, linewidth=1.2)

    ylims!(ax, 0.92, 1.15)
    axislegend(ax, position=:lt, framevisible=true, labelsize=11)
    save(savepath, fig)
    println("  → Saved: $savepath")
    return fig
end

# ── Plot 2: STATCOM rating sweep ─────────────────────────────────────────────
function plot_rating_sweep(rating_results; savepath="./plots/rating_sweep.pdf")

    ratings  = Float64[d.rating_kvar for d in rating_results]
    utils    = Float64[d.util        for d in rating_results]
    baseline = rating_results[1].baseline

    fig = Figure(size=(820, 500), backgroundcolor=:white)
    ax  = Axis(fig[1,1],
        xlabel       = "STATCOM rating (kVAr per unit, total nameplate)",
        ylabel       = "PV utilisation (%)",
        title        = "PV hosting capacity vs. STATCOM rating (52 units, single-phase PV)",
        titlesize    = 15, xlabelsize = 13, ylabelsize = 13,
        xscale       = log10,
        xticks       = (ratings, ["$(Int(r))" for r in ratings]),
        xticklabelrotation = π/6,
        ygridvisible = true, xgridvisible = false,
        ygridcolor   = (:black, 0.08),
    )

    hlines!(ax, [baseline], color=:red, linestyle=:dash, linewidth=2.0,
            label="No STATCOM baseline ($(round(baseline, digits=1))%)")
    hlines!(ax, [100.0], color=(:forestgreen, 0.6), linestyle=:dash, linewidth=2.0,
            label="100% utilisation")

    band!(ax, ratings, fill(baseline, length(ratings)), utils,
          color=(:darkorange, 0.14))
    lines!(ax,   ratings, utils, color=:darkorange, linewidth=3.0)
    scatter!(ax, ratings, utils, color=:darkorange, markersize=10,
             label="52 STATCOMs, uniform spacing")

    for (r, u) in zip(ratings, utils)
        text!(ax, r, u + 1.2,
              text="$(round(u, digits=1))%",
              fontsize=8, align=(:center, :bottom), color=(:darkorange, 0.85))
    end

    axislegend(ax, position=:lt, framevisible=true, labelsize=11)
    ylims!(ax, min(65, minimum(utils)-5), 110)

    save(savepath, fig)
    println("  → Saved: $savepath")
    return fig
end

# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════
do_case1 = true
do_case2 = true
do_case3 = true

const PV_COST      = -1000.0
const STATCOM_COST =  1.0

# Case 2: sweep of realistic single-phase inverter sizes (kW)
const PV_KW_LEVELS = [1.0, 3.0, 5.0, 7.0, 10.0]

# Case 3: fixed PV size for hosting-capacity stress test, chosen to sit
# inside typical UK single-phase residential inverter range (≤5kW is most
# common; higher values push into stress-test territory deliberately)
const CASE3_PV_KW = 5.0

# Case 3: STATCOM rating sweep, kVAr TOTAL nameplate per unit
const STATCOM_RATINGS_KVAR = [1.0, 5.0, 10.0, 20.0, 30.0, 50.0, 70.0, 100.0, 200.0]

data_eng_ref = PMD.parse_file(data_path, transformations=[PMD.transform_loops!])
dm_ref = load_base_network(data_path)
summarise_network(dm_ref)
report_load_phase_distribution(dm_ref)

SBASE_KVA, VBASE_KV = get_base_values(data_eng_ref, dm_ref)
println("  Base values: sbase = $(round(SBASE_KVA,digits=1)) kVA")

# ───────────────────────────────────────────────────────────────
## Case 1: Natural baseline — unchanged
# ───────────────────────────────────────────────────────────────
if do_case1
    println("\n" * "="^55)
    println(" CASE 1: Natural Baseline (no bounds, no DER)")
    println("="^55)
    dm1 = load_base_network(data_path, enforce_bounds=false)
    r1  = solve_and_report(dm1, "Baseline")
end

# ───────────────────────────────────────────────────────────────
## Case 2: PV penetration sweep — now over real kW ratings, single-phase,
##         phase-matched to the colocated load.
# ───────────────────────────────────────────────────────────────
if do_case2
    println("\n" * "="^55)
    println(" CASE 2: PV Penetration Sweep (no bounds)")
    println(" Single-phase PV, phase-matched to colocated load")
    println("="^55)

    for pv_kw in PV_KW_LEVELS
        println("\n  ── pv_kw = $(pv_kw) kW/unit ──")
        dm = load_base_network(data_path, enforce_bounds=false)
        add_pv!(dm; pv_kw=pv_kw, q_scale=1.25, spacing=1, pv_cost=PV_COST)
        solve_and_report(dm, "PV $(pv_kw) kW/unit  [no bounds]")
    end
end

# ───────────────────────────────────────────────────────────────
## Case 3: STATCOM hosting capacity study — single-phase PV fixed at
##         CASE3_PV_KW, STATCOM rating swept in real kVAr.
# ───────────────────────────────────────────────────────────────
if do_case3
    println("\n" * "="^55)
    println(" CASE 3: STATCOM Hosting Capacity Study")
    println(" PV: $(CASE3_PV_KW) kW/unit, single-phase, phase-matched  |  bounds: 0.90–1.10 pu")
    println(" PV cost: $(PV_COST)  |  STATCOM cost: $(STATCOM_COST)")
    println(" Builder: PMD.build_mc_opf  |  Objective: cost (gen[\"cost\"])")
    println("="^55)

    println("\n  ── 3a: PV only — hosting capacity baseline ──")
    dm3a = load_base_network(data_path, enforce_bounds=true)
    add_pv!(dm3a; pv_kw=CASE3_PV_KW, q_scale=1.25, spacing=1, pv_cost=PV_COST)
    r3a  = solve_and_report(dm3a, "PV only  [no STATCOM]")

    println("\n  ── 3b: STATCOM rating sweep (spacing=10) ──")
    rating_results = NamedTuple[]
    r3b_last = nothing
    for rating in STATCOM_RATINGS_KVAR
        local dm3b = load_base_network(data_path, enforce_bounds=true)
        add_pv!(dm3b; pv_kw=CASE3_PV_KW, q_scale=1.25, spacing=1, pv_cost=PV_COST)
        
        # 1. Capture the returned array to count the units placed
        statcom_ids = add_statcoms!(dm3b; rating_kvar=rating, spacing=10, statcom_cost=STATCOM_COST)
        n_st = length(statcom_ids)
        
        println("\n  ── STATCOM rating = $(rating) kVAr/unit  ($(n_st) units) ──")
        
        # Update the report label so the console logs match the unit count
        r = solve_and_report(dm3b, "PV+STATCOM $(rating) kVAr/unit  $(n_st) units")
        
        # 2. Add n_units to the NamedTuple that gets pushed to the results array
        push!(rating_results, (rating_kvar=rating, util=r.pv_util, baseline=r3a.pv_util, n_units=n_st))
        
        global r3b_last = r
    end
    
    plot_rating_sweep(rating_results)

    if !isnothing(r3b_last)
        plot_voltage_profile(r3a, r3b_last, CASE3_PV_KW, STATCOM_RATINGS_KVAR[end])
    end

    # 3. Call r.n_units dynamically in the table row creation
    table_rows = vcat(
        [(label="No STATCOM (baseline)", util=r3a.pv_util, curtail=r3a.pv_curtail)],
        [(label="STATCOM $(Int(r.rating_kvar)) kVAr ($(r.n_units) units)", util=r.util, curtail=100.0 - r.util)
         for r in rating_results]
    )
    print_delta_hc_table(table_rows)
end

println("\n  Plots saved to ./plots/")
println("  Done.")