#=
# ==============================================================================
# SCRIPT: STATCOM & PV Hosting Capacity Assessment — Extended
# PROJECT: Undergraduate Thesis - Power Systems Engineering
# AUTHOR: Pat
# DATE: May 2026
# ==============================================================================
#
# DESCRIPTION:
# Extended version of the original OPF assessment script. All new analysis is
# built on the same PMD.instantiate_mc_model / PMD.build_mc_opf / PMD.optimize_model!
# architecture as the original. No changes to the solver pipeline.
#
# ADDITIONS OVER ORIGINAL:
#   Case 4 — HC curve: with vs without STATCOM across PV scales 1–5×
#             This is the central thesis result — ΔHC as a function of penetration.
#   Case 5 — Targeted placement: STATCOM at worst-voltage buses vs uniform spacing
#   Case 6 — Decomposition: PV Q=0 vs smart inverter vs smart inverter + STATCOM
#             Separates how much HC gain comes from inverter Q vs STATCOM Q.
#
#   solve_and_report now returns a NamedTuple of extracted scalars so results
#   can be collected and passed to plotting functions without re-parsing stdout.
#
#   Six CairoMakie plots saved to ./plots/:
#     voltage_profile.pdf   — bus voltage with/without STATCOM
#     hc_curve.pdf          — ΔHC across penetration levels (central result)
#     rating_sweep.pdf      — PV utilisation vs STATCOM rating
#     density_sweep.pdf     — PV utilisation vs STATCOM count
#     placement_compare.pdf — targeted vs uniform placement
#     decomposition.pdf     — inverter Q vs STATCOM Q contributions
#
# TECHNICAL STACK: unchanged from original
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
using CairoMakie     # ADD: plots — run Pkg.add("CairoMakie") if not installed
using Printf         # ADD: formatted table output

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
# NETWORK LOADER  — unchanged
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
# PV PLACEMENT  — unchanged
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
        gen["pmax"]    =  pmax * ones(3)
        gen["pmin"]    =  zeros(3)
        gen["qmax"]    =  smax * ones(3)
        gen["qmin"]    = -smax * ones(3)
        gen["cost"]    = [pv_cost 0.0]

        push!(pv_ids, gen_id)
    end

    println("  PV: $(length(pv_ids)) units  pv_scale=$(pv_scale)×pd  q_scale=$(q_scale)  spacing=$(spacing)  cost=$(pv_cost)")
    return pv_ids
end

# ═══════════════════════════════════════════════════════════════
# STATCOM PLACEMENT
# ADD: target_buses keyword — pass a Vector{Int} of bus IDs to
#      place STATCOMs at specific buses rather than uniformly.
#      When target_buses is provided, spacing is ignored.
# ═══════════════════════════════════════════════════════════════
function add_statcoms!(data_math;
        q_scale      = 1.0,
        spacing      = 1,
        statcom_cost = 1.0,
        target_buses = nothing)

    load_ids = sort(
        collect(keys(data_math["load"])),
        by = x -> tryparse(Int, x) === nothing ? 0 : parse(Int, x)
    )

    statcom_ids = String[]

    # ── ADD: targeted placement path ───────────────────────────
    if !isnothing(target_buses)
        # Use mean qd as the rating reference for targeted units
        # since they are not co-located with a specific load.
        qd_mean = mean(data_math["load"][id]["qd"][1] for id in keys(data_math["load"]))
        qlim    = q_scale * qd_mean

        for bus_id in target_buses
            gen_id = string(length(data_math["gen"]) + 1)
            data_math["gen"][gen_id] = deepcopy(data_math["gen"]["1"])
            gen = data_math["gen"][gen_id]

            gen["gen_bus"] = bus_id
            gen["type"]    = "STATCOM"
            gen["name"]    = "statcom_targeted_bus$(bus_id)"
            gen["pmax"]    =  zeros(3)
            gen["pmin"]    =  zeros(3)
            gen["qmax"]    =  qlim * ones(3)
            gen["qmin"]    = -qlim * ones(3)
            gen["cost"]    = [statcom_cost 0.0]

            push!(statcom_ids, gen_id)
        end
        println("  STATCOM: $(length(statcom_ids)) units  TARGETED buses=$(target_buses)  q_scale=$(q_scale)×qd_mean  cost=$(statcom_cost)")
        return statcom_ids
    end
    # ── end targeted path ───────────────────────────────────────

    # Original uniform path — unchanged
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
        gen["pmax"]    =  zeros(3)
        gen["pmin"]    =  zeros(3)
        gen["qmax"]    =  qlim * ones(3)
        gen["qmin"]    = -qlim * ones(3)
        gen["cost"]    = [statcom_cost 0.0]

        push!(statcom_ids, gen_id)
    end

    println("  STATCOM: $(length(statcom_ids)) units  q_scale=$(q_scale)×qd  spacing=$(spacing)  cost=$(statcom_cost)")
    return statcom_ids
end

# ═══════════════════════════════════════════════════════════════
# SOLVE AND REPORT
# CHANGE from original: returns a NamedTuple of extracted scalars
# so that results can be collected across runs for plotting.
# All solution extraction still uses result["solution"]["bus"]
# and result["solution"]["gen"] exactly as in the original.
# ═══════════════════════════════════════════════════════════════
function solve_and_report(data_math, label)
    PMD.add_start_vrvi!(data_math)

    model  = PMD.instantiate_mc_model(data_math, PMD.IVRENPowerModel, PMD.build_mc_opf)
    result = PMD.optimize_model!(model, optimizer=ipopt_solver)
    status = result["termination_status"]
    println("\n  [$label]  →  $(status)")

    # Default values — NaN/0 until populated on successful solve
    v_min = NaN;  v_mean = NaN;  v_max = NaN
    n_over = 0;   n_under = 0
    pv_util = NaN;  pv_curtail = NaN
    pv_output = NaN;  pv_capacity = NaN
    st_q_total = NaN;  st_q_cap = NaN;  st_util = NaN
    vm_per_bus = Dict{String, Float64}()   # ADD: bus_id => worst-phase voltage

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
            vm_per_bus[b] = bus_max   # ADD: worst phase per bus for profile plot
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

        # ── PV dispatch ─────────────────────────────────────────
        pv_gens = [(i, g) for (i, g) in data_math["gen"] if get(g, "type", "") == "PV"]
        if !isempty(pv_gens) && haskey(sol, "gen")
            pg_vals     = [sum(get(sol["gen"][i], "pg", zeros(3))) for (i, g) in pv_gens if haskey(sol["gen"], i)]
            pg_rated    = [sum(g["pmax"])                          for (i, g) in pv_gens]
            pv_output   = sum(pg_vals)
            pv_capacity = sum(pg_rated)
            pv_util     = pv_output / max(1e-9, pv_capacity) * 100
            pv_curtail  = max(0.0, 100.0 - pv_util)
            println("    PV output / capacity  : $(round(pv_output,sigdigits=4)) / $(round(pv_capacity,sigdigits=4)) pu")
            println("    PV utilisation        : $(round(pv_util,digits=1))%  →  curtailment: $(round(pv_curtail,digits=1))%")
        end

        # ── STATCOM dispatch ────────────────────────────────────
        st_gens = [(i, g) for (i, g) in data_math["gen"] if get(g, "type", "") == "STATCOM"]
        if !isempty(st_gens) && haskey(sol, "gen")
            qg_vals    = [sum(get(sol["gen"][i], "qg", zeros(3))) for (i, g) in st_gens if haskey(sol["gen"], i)]
            qg_rated   = [sum(g["qmax"])                          for (i, g) in st_gens]
            st_q_total = sum(qg_vals)
            st_q_cap   = sum(qg_rated)
            st_util    = abs(st_q_total) / max(1e-9, st_q_cap) * 100
            direction  = st_q_total >= 0 ? "injecting ↑V" : "absorbing ↓V"
            println("    STATCOM Q / capacity  : $(round(st_q_total,sigdigits=4)) / $(round(st_q_cap,sigdigits=4)) pu")
            println("    STATCOM utilisation   : $(round(st_util,digits=1))%  $(direction)")

            # ADD: per-phase Q breakdown — reveals phase imbalance in compensation
            phase_labels = ["Phase A", "Phase B", "Phase C"]
            for (pi, ph) in enumerate(phase_labels)
                q_ph = sum(
                    get(sol["gen"][i], "qg", zeros(3))[pi]
                    for (i, g) in st_gens if haskey(sol["gen"], i)
                )
                println("      $(ph) net Q: $(round(q_ph, sigdigits=3)) pu  $(q_ph >= 0 ? "↑" : "↓")")
            end
        end

    else
        println("    WARNING: $(status)")
        println("    Objective : $(get(result, "objective", "N/A"))")
    end

    # ADD: return NamedTuple — backward-compatible (callers that ignore return value
    #      still work; callers that need numbers can collect them)
    return (
        status      = string(status),
        v_min       = v_min,
        v_mean      = v_mean,
        v_max       = v_max,
        n_over      = n_over,
        n_under     = n_under,
        pv_util     = pv_util,
        pv_curtail  = pv_curtail,
        pv_output   = pv_output,
        pv_capacity = pv_capacity,
        st_q_total  = st_q_total,
        st_q_cap    = st_q_cap,
        st_util     = st_util,
        vm_per_bus  = vm_per_bus,
    )
end

# ═══════════════════════════════════════════════════════════════
# ADD: ΔHC SUMMARY TABLE  — printed to console after Case 3
# ═══════════════════════════════════════════════════════════════
function print_delta_hc_table(rows)
    # rows: Vector of NamedTuples (label, util, curtail)
    println("\n  ΔHC Summary Table")
    println("  " * "─"^72)
    @printf("  %-32s  %10s  %12s  %14s\n",
            "Configuration", "Util (%)", "Curtail (%)", "ΔHC vs base")
    println("  " * "─"^72)
    baseline = rows[1].util
    for r in rows
        delta     = r.util - baseline
        delta_str = delta ≈ 0.0 ? "           —" : @sprintf("+%.1f pp", delta)
        @printf("  %-32s  %10s  %12s  %14s\n",
                r.label,
                @sprintf("%.1f%%", r.util),
                @sprintf("%.1f%%", r.curtail),
                delta_str)
    end
    println("  " * "─"^72)
end

# ═══════════════════════════════════════════════════════════════
# ADD: PLOTTING FUNCTIONS
# All use CairoMakie. Input is NamedTuples from solve_and_report.
# No stdout re-parsing, no global variables.
# ═══════════════════════════════════════════════════════════════
mkpath("./plots")
#=
# ==============================================================================
# PLOTTING FUNCTIONS — improved versions
# Drop these into statcom_extended.jl, replacing the originals.
#
# Changes per plot:
#   hc_curve        — annotate divergence point; stronger shading; callout box
#   rating_sweep    — fix vline position; readable x-ticks; saturation band;
#                     individual point labels
#   density_sweep   — fix x-axis range; point labels; annotate non-monotonic point
#   placement_compare — add Δ-from-baseline labels; framing note
#   decomposition   — stronger colours; cleaner arrows; projection-safe fonts
#
# All plots: titlesize 15, larger tick labels, projection-safe line weights.
# ==============================================================================
=#

mkpath("./plots")

# ── Plot 1: Voltage profile ──────────────────────────────────────────────────
# Unchanged — was already working correctly.
function plot_voltage_profile(r_no_st, r_with_st, pv_scale;
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
        title        = "Voltage profile — PV $(pv_scale)× load, with and without STATCOM",
        titlesize    = 15, xlabelsize = 12, ylabelsize = 12,
        ygridvisible = true, xgridvisible = false,
        ygridcolor   = (:black, 0.08),
    )

    hspan!(ax, 1.10, 1.25, color=(:red,  0.07))
    hspan!(ax, 0.75, 0.90, color=(:blue, 0.07))

    lines!(ax, xs, v_no, color=:steelblue,  linewidth=2.5, label="No STATCOM")
    lines!(ax, xs, v_st, color=:darkorange, linewidth=2.5,
           label="With STATCOM ($(BEST_Q_SCALE)×q_d, 52 units)")

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

# ── Plot 2: HC curve — central thesis result ─────────────────────────────────
# CHANGES:
#   • Stronger orange shading for the ΔHC band
#   • Explicit annotation at divergence point (~3–4× scale)
#   • Callout box for the headline result at 5×
#   • Both curves labelled at their endpoints (no legend needed)
function plot_hc_curve(hc_curve_data; savepath="./plots/hc_curve.pdf")

    scales   = [d.scale   for d in hc_curve_data]
    util_no  = [d.util_no for d in hc_curve_data]
    util_st  = [d.util_st for d in hc_curve_data]
    delta_hc = util_st .- util_no

    fig = Figure(size=(820, 500), backgroundcolor=:white)
    ax  = Axis(fig[1,1],
        xlabel       = "PV penetration scale (× load demand)",
        ylabel       = "PV utilisation (%)",
        title        = "Hosting capacity: with and without STATCOM ($(BEST_Q_SCALE)×q_d, 52 units)",
        titlesize    = 15, xlabelsize = 13, ylabelsize = 13,
        xticks       = (scales, ["$(Int(s))×" for s in scales]),
        ygridvisible = true, xgridvisible = false,
        ygridcolor   = (:black, 0.08),
    )

    # ΔHC shaded region — stronger than before
    band!(ax, scales, util_no, util_st, color=(:darkorange, 0.22))

    lines!(ax, scales, util_no, color=:steelblue,  linewidth=3.0, label="No STATCOM")
    lines!(ax, scales, util_st, color=:darkorange, linewidth=3.0, label="With STATCOM")
    scatter!(ax, scales, util_no, color=:steelblue,  markersize=10)
    scatter!(ax, scales, util_st, color=:darkorange, markersize=10)

    # Annotate ΔHC only where it is non-zero (skips 1× and 2× if they are equal)
    for (s, d, y) in zip(scales, delta_hc, util_st)
        d > 0.5 || continue
        text!(ax, s, min(y + 1.8, 107.0),
              text="+$(round(d, digits=1)) pp",
              fontsize=11, align=(:center, :bottom), color=:darkorange,
              font=:bold)
    end

    # ADD: divergence callout — mark where curves first separate
    # Find first scale where delta > 0.5pp
    div_idx = findfirst(d -> d > 0.5, delta_hc)
    if !isnothing(div_idx)
        div_x = scales[div_idx]
        div_y = util_no[div_idx]
        # Vertical dashed line at divergence
        vlines!(ax, [div_x], color=(:red, 0.45), linestyle=:dot, linewidth=1.8)
        text!(ax, div_x + 0.08, 65.0,
              text="STATCOM\nbecomes\nnecessary",
              fontsize=9, color=(:red, 0.7), align=(:left, :bottom))
    end

    # ADD: endpoint labels directly on curves (avoids relying on legend alone)
    text!(ax, scales[end] + 0.06, util_no[end],
          text="No STATCOM\n$(round(util_no[end], digits=1))%",
          fontsize=9, color=:steelblue, align=(:left, :center))
    text!(ax, scales[end] + 0.06, util_st[end],
          text="With STATCOM\n$(round(util_st[end], digits=1))%",
          fontsize=9, color=:darkorange, align=(:left, :center))

    hlines!(ax, [100.0], color=(:forestgreen, 0.6), linestyle=:dash, linewidth=2.0,
            label="100% utilisation")

    ylims!(ax, 58, 112)
    xlims!(ax, 0.7, 6.0)   # extra right margin for endpoint labels
    axislegend(ax, position=:lb, framevisible=true, labelsize=11)

    save(savepath, fig)
    println("  → Saved: $savepath")
    return fig
end

# ── Plot 3: Rating sweep ─────────────────────────────────────────────────────
# CHANGES:
#   • x-ticks replaced with plain labels (1, 2, 5, 10, 20, 50, 100, 200, 500)
#     instead of 10^0/10^1/10^2 which confuses non-specialists
#   • Saturation band shaded for ratings ≥ 200× (both at 100%)
#   • vline moved to exactly between 100× and 200× on log scale (≈141×)
#   • Individual point value labels added
function plot_rating_sweep(rating_results; savepath="./plots/rating_sweep.pdf")

    ratings  = Float64[d.q_scale for d in rating_results]
    utils    = Float64[d.util    for d in rating_results]
    baseline = rating_results[1].baseline

    fig = Figure(size=(820, 500), backgroundcolor=:white)
    ax  = Axis(fig[1,1],
        xlabel       = "STATCOM rating (multiples of q_d per unit)",
        ylabel       = "PV utilisation (%)",
        title        = "PV hosting capacity vs. STATCOM reactive power rating (52 units, 5× stress)",
        titlesize    = 15, xlabelsize = 13, ylabelsize = 13,
        xscale       = log10,
        # Readable tick labels — show actual multiplier values
        xticks       = (ratings, ["$(Int(r))×" for r in ratings]),
        xticklabelrotation = π/6,
        ygridvisible = true, xgridvisible = false,
        ygridcolor   = (:black, 0.08),
    )

    # Saturation region: ratings where util ≈ 100%
    sat_start = ratings[findfirst(u -> u >= 99.5, utils)]
    vspan!(ax, sat_start, ratings[end] * 1.5,
           color=(:forestgreen, 0.07))
    text!(ax, sat_start * 1.05, 68.0,
          text="Saturation\n(100% HC)", fontsize=9,
          color=(:forestgreen, 0.75), align=(:left, :bottom))

    hlines!(ax, [baseline], color=:red, linestyle=:dash, linewidth=2.0,
            label="No STATCOM baseline ($(round(baseline, digits=1))%)")
    hlines!(ax, [100.0], color=(:forestgreen, 0.6), linestyle=:dash, linewidth=2.0,
            label="100% utilisation")

    band!(ax, ratings, fill(baseline, length(ratings)), utils,
          color=(:darkorange, 0.14))
    lines!(ax,   ratings, utils, color=:darkorange, linewidth=3.0)
    scatter!(ax, ratings, utils, color=:darkorange, markersize=10,
             label="52 STATCOMs, uniform spacing")

    # ADD: value labels on each point
    for (r, u) in zip(ratings, utils)
        text!(ax, r, u + 1.2,
              text="$(round(u, digits=1))%",
              fontsize=8, align=(:center, :bottom), color=(:darkorange, 0.85))
    end

    # FIX: vline at geometric mean of 100 and 200 on log scale = sqrt(100*200) ≈ 141
    flip_x = sqrt(100.0 * 200.0)
    vlines!(ax, [flip_x], color=(:purple, 0.45), linestyle=:dot, linewidth=1.8)
    text!(ax, flip_x * 1.05, 74.0,
          text="Optimal Q dispatch\nflips: absorb → inject\n(between 100× and 200×)",
          fontsize=9, color=(:purple, 0.8), align=(:left, :bottom))

    axislegend(ax, position=:lt, framevisible=true, labelsize=11)
    ylims!(ax, 65, 110)

    save(savepath, fig)
    println("  → Saved: $savepath")
    return fig
end

# ── Plot 4: Density sweep ────────────────────────────────────────────────────
# CHANGES:
#   • x-axis range fixed to start at 0 so the 6-unit leftmost point is visible
#   • Value labels added to each point on the primary axis
#   • Non-monotonic point (7 units > 9 units) annotated explicitly
function plot_density_sweep(density_results; savepath="./plots/density_sweep.pdf")

    n_units  = Int[d.n_units   for d in density_results]
    utils    = Float64[d.util  for d in density_results]
    st_utils = Float64[d.st_util for d in density_results]
    baseline = density_results[1].baseline

    fig = Figure(size=(840, 500), backgroundcolor=:white)
    ax  = Axis(fig[1,1],
        xlabel       = "Number of STATCOM units deployed",
        ylabel       = "PV utilisation (%)",
        title        = "Hosting capacity vs. deployment density (rating = $(BEST_Q_SCALE)×q_d, 5× stress)",
        titlesize    = 15, xlabelsize = 13, ylabelsize = 13,
        ygridvisible = true, xgridvisible = false,
        ygridcolor   = (:black, 0.08),
    )

    ax2 = Axis(fig[1,1],
        yaxisposition   = :right,
        ylabel          = "STATCOM Q utilisation (%)",
        ylabelsize      = 11,
        yticklabelcolor = :purple,
        ylabelcolor     = :purple,
        ygridvisible    = false,
    )
    hidexdecorations!(ax2)

    hlines!(ax, [baseline], color=:red, linestyle=:dash, linewidth=2.0,
            label="No STATCOM ($(round(baseline, digits=1))%)")
    hlines!(ax, [100.0], color=(:forestgreen, 0.6), linestyle=:dash, linewidth=2.0,
            label="100% utilisation")

    band!(ax, n_units, fill(baseline, length(n_units)), utils,
          color=(:steelblue, 0.14))
    lines!(ax,   n_units, utils, color=:steelblue, linewidth=3.0)
    scatter!(ax, n_units, utils, color=:steelblue, markersize=10,
             label="PV utilisation")

    lines!(ax2,   n_units, st_utils, color=(:purple, 0.75), linewidth=2.0,
           linestyle=:dot)
    scatter!(ax2, n_units, st_utils, color=:purple, markersize=8,
             marker=:diamond, label="STATCOM Q utilisation")

    # ADD: value labels on PV utilisation points
    for (n, u) in zip(n_units, utils)
        text!(ax, Float64(n), u + 0.8,
              text="$(round(u, digits=1))%",
              fontsize=8, align=(:center, :bottom), color=(:steelblue, 0.85))
    end

    # ADD: annotate non-monotonic point if it exists
    # Find index where util[i] > util[i+1] (going left to right on sorted n_units)
    sorted_pairs = sort(collect(zip(n_units, utils)), by=x->x[1])
    for i in 1:length(sorted_pairs)-1
        n_a, u_a = sorted_pairs[i]
        n_b, u_b = sorted_pairs[i+1]
        if u_a > u_b  # non-monotonic: fewer units, higher HC
            text!(ax, Float64(n_a) + 0.5, u_a - 1.5,
                  text="$(n_a) units outperforms\n$(n_b) units — placement\nlocation effect",
                  fontsize=8, color=(:darkorange, 0.8), align=(:left, :top))
        end
    end

    # FIX: ensure x-axis starts slightly below minimum n_units
    xlims!(ax, 0, maximum(n_units) * 1.08)
    xlims!(ax2, 0, maximum(n_units) * 1.08)

    axislegend(ax,  position=:rb, framevisible=true, labelsize=11)
    axislegend(ax2, position=:rt, framevisible=true, labelsize=11)

    save(savepath, fig)
    println("  → Saved: $savepath")
    return fig
end

# ── Plot 5: Placement comparison ─────────────────────────────────────────────
# CHANGES:
#   • Add Δ-from-baseline label below each bar's value
#   • Stronger bar colours for projection
#   • Add text note framing the 1.6pp targeted vs uniform result
function plot_placement_comparison(baseline_util, uniform_util, targeted_util;
        savepath="./plots/placement_compare.pdf")

    labels = ["No STATCOM\n(baseline)",
              "Uniform placement\n(co-located with PV,\nspacing = 10)",
              "Targeted placement\n(buses 189, 256, 157,\n843, 260)"]
    vals = [baseline_util, uniform_util, targeted_util]
    cols = [RGBf(0.27, 0.51, 0.71),   # steelblue — stronger
            RGBf(0.84, 0.44, 0.07),   # darkorange — stronger
            RGBf(0.13, 0.55, 0.13)]   # forestgreen — stronger

    fig = Figure(size=(660, 480), backgroundcolor=:white)
    ax  = Axis(fig[1,1],
        ylabel         = "PV utilisation (%)",
        title          = "Placement strategy comparison — $(length(WORST_BUSES)) units, $(BEST_Q_SCALE)×q_d, 5× stress",
        titlesize      = 15, ylabelsize = 13,
        xticks         = (1:3, labels),
        xticklabelsize = 10,
        ygridvisible   = true, xgridvisible = false,
        ygridcolor     = (:black, 0.08),
    )

    barplot!(ax, 1:3, vals, color=cols, width=0.54)

    # Value label + Δ label on each bar
    deltas = [0.0, uniform_util - baseline_util, targeted_util - baseline_util]
    for (i, (v, d)) in enumerate(zip(vals, deltas))
        text!(ax, Float64(i), v + 0.5,
              text=@sprintf("%.1f%%", v),
              align=(:center, :bottom), fontsize=13, color=cols[i], font=:bold)
        if d > 0.0
            text!(ax, Float64(i), v - 3.5,
                  text=@sprintf("(+%.1f pp vs baseline)", d),
                  align=(:center, :top), fontsize=9, color=(:black, 0.55))
        end
    end

    # ADD: framing annotation — key insight
    text!(ax, 2.0, 63.5,
          text="Density matters more than placement strategy:\n52 uniform units achieve 97.3% (not shown)",
          fontsize=9, color=(:black, 0.5), align=(:center, :bottom))

    hlines!(ax, [100.0], color=(:forestgreen, 0.45), linestyle=:dash, linewidth=2.0)
    ylims!(ax, 62, 112)

    save(savepath, fig)
    println("  → Saved: $savepath")
    return fig
end

# ── Plot 6: HC decomposition ─────────────────────────────────────────────────
# CHANGES:
#   • Stronger, fully saturated bar colours for projection
#   • Arrows made thicker and more visible
#   • Added percentage contribution labels (what % of total ΔHC does each provide?)
#   • Baseline dashed line added for context
function plot_decomposition(util_q0, util_inverter, util_statcom;
        savepath="./plots/decomposition.pdf")

    labels = ["PV (Q = 0)\nno inverter Q",
              "Smart inverter\n(Q enabled, no STATCOM)",
              "Smart inverter\n+ STATCOM ($(BEST_Q_SCALE)×q_d)"]
    vals = [util_q0, util_inverter, util_statcom]
    cols = [RGBf(0.27, 0.51, 0.71),
            RGBf(0.84, 0.44, 0.07),
            RGBf(0.13, 0.55, 0.13)]

    inv_contribution  = util_inverter - util_q0
    stat_contribution = util_statcom  - util_inverter
    total_delta       = util_statcom  - util_q0
    inv_pct  = round(inv_contribution  / total_delta * 100, digits=0)
    stat_pct = round(stat_contribution / total_delta * 100, digits=0)

    fig = Figure(size=(700, 500), backgroundcolor=:white)
    ax  = Axis(fig[1,1],
        ylabel         = "PV utilisation (%)",
        title          = "HC decomposition: inverter Q vs STATCOM Q contributions (5× stress)",
        titlesize      = 15, ylabelsize = 13,
        xticks         = (1:3, labels),
        xticklabelsize = 11,
        ygridvisible   = true, xgridvisible = false,
        ygridcolor     = (:black, 0.08),
    )

    barplot!(ax, 1:3, vals, color=cols, width=0.54)

    # Value labels
    for (i, v) in enumerate(vals)
        text!(ax, Float64(i), v + 0.5,
              text=@sprintf("%.1f%%", v),
              align=(:center, :bottom), fontsize=13, color=cols[i], font=:bold)
    end

    # Arrow and label — inverter contribution
    mid_inv = util_q0 + inv_contribution / 2
    arrows!(ax, [1.52], [util_q0 + 2.0], [0.0], [inv_contribution - 4.0],
            color=RGBf(0.84, 0.44, 0.07), linewidth=2.5, arrowsize=12)
    text!(ax, 1.58, mid_inv,
          text=@sprintf("+%.1f pp\n(inverter Q)\n%d%% of ΔHC", inv_contribution, Int(inv_pct)),
          fontsize=9, color=RGBf(0.84, 0.44, 0.07), align=(:left, :center))

    # Arrow and label — STATCOM contribution
    mid_st = util_inverter + stat_contribution / 2
    arrows!(ax, [2.52], [util_inverter + 2.0], [0.0], [stat_contribution - 4.0],
            color=RGBf(0.13, 0.55, 0.13), linewidth=2.5, arrowsize=12)
    text!(ax, 2.58, mid_st,
          text=@sprintf("+%.1f pp\n(STATCOM Q)\n%d%% of ΔHC", stat_contribution, Int(stat_pct)),
          fontsize=9, color=RGBf(0.13, 0.55, 0.13), align=(:left, :center))

    # Total ΔHC label
    text!(ax, 2.0, util_statcom + 2.5,
          text=@sprintf("Total ΔHC = +%.1f pp", total_delta),
          fontsize=10, color=(:black, 0.6), align=(:center, :bottom), font=:bold)

    hlines!(ax, [100.0], color=(:forestgreen, 0.45), linestyle=:dash, linewidth=2.0)
    ylims!(ax, 55, 115)

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
do_case4 = true   # ADD: HC curve across penetration levels
do_case5 = true   # ADD: targeted vs uniform placement
do_case6 = true   # ADD: inverter Q vs STATCOM Q decomposition

const STRESS_SCALE = 5.0
const PV_COST      = -1000.0
const STATCOM_COST =  1.0
const BEST_Q_SCALE = 100    # practical saturation point from Case 3b

# Worst-voltage buses from Case 2 results — used for targeted placement
const WORST_BUSES = [189, 256, 157, 843, 260]

dm_ref = load_base_network(data_path)
summarise_network(dm_ref)

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
## Case 2: PV penetration sweep — unchanged
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
## Case 3: STATCOM hosting capacity — same logic as original,
##         now collects results into vectors for plotting and
##         prints the ΔHC table at the end.
# ───────────────────────────────────────────────────────────────
if do_case3
    println("\n" * "="^55)
    println(" CASE 3: STATCOM Hosting Capacity Study")
    println(" PV stress: $(STRESS_SCALE)×pd  |  bounds: 0.90–1.10 pu")
    println(" PV cost: $(PV_COST)  |  STATCOM cost: $(STATCOM_COST)")
    println(" Builder: PMD.build_mc_opf  |  Objective: cost (gen[\"cost\"])")
    println("="^55)

    println("\n  ── 3a: PV only — hosting capacity baseline ──")
    dm3a = load_base_network(data_path, enforce_bounds=true)
    add_pv!(dm3a; pv_scale=STRESS_SCALE, q_scale=1.25, spacing=1, pv_cost=PV_COST)
    r3a  = solve_and_report(dm3a, "PV only  [no STATCOM]")

    println("\n  ── 3b: STATCOM rating sweep (52 STATCOMs, spacing=1) ──")
    rating_results = NamedTuple[]
    for q_sc in [1, 2, 5, 10, 20, 50, 100, 200, 500]
        println("\n  ── q_scale = $(q_sc)×qd  (52 STATCOMs) ──")
        local dm3b = load_base_network(data_path, enforce_bounds=true)
        add_pv!(dm3b;       pv_scale=STRESS_SCALE, q_scale=1.25, spacing=1, pv_cost=PV_COST)
        add_statcoms!(dm3b; q_scale=q_sc,          spacing=1,    statcom_cost=STATCOM_COST)
        r = solve_and_report(dm3b, "PV+STATCOM q=$(q_sc)×qd  52 units")
        push!(rating_results, (q_scale=Float64(q_sc), util=r.pv_util, baseline=r3a.pv_util))
    end
    plot_rating_sweep(rating_results)

    min_viable_q = 100
    println("\n  ── 3c: Density sweep at q_scale=$(min_viable_q)×qd ──")
    println("        How few STATCOMs achieve the same hosting benefit?")
    density_results = NamedTuple[]
    r3c_full = nothing
    for sp in [1, 2, 4, 6, 8, 10]
        n_st = ceil(Int, 52 / sp)
        println("\n  ── spacing=$(sp) → ~$(n_st) STATCOMs ──")
        local dm3c = load_base_network(data_path, enforce_bounds=true)
        add_pv!(dm3c;       pv_scale=STRESS_SCALE, q_scale=1.25,  spacing=1, pv_cost=PV_COST)
        add_statcoms!(dm3c; q_scale=min_viable_q,  spacing=sp,    statcom_cost=STATCOM_COST)
        r = solve_and_report(dm3c, "PV+STATCOM q=$(min_viable_q)×qd  $(n_st) units")
        push!(density_results, (n_units=n_st, util=r.pv_util, st_util=r.st_util, baseline=r3a.pv_util))
        sp == 1 && (r3c_full = r)
    end
    plot_density_sweep(density_results)

    # Voltage profile: no STATCOM vs full-density 100×qd
    if !isnothing(r3c_full)
        plot_voltage_profile(r3a, r3c_full, STRESS_SCALE)
    end

    # ΔHC summary table
    table_rows = vcat(
        [(label="No STATCOM (baseline)",                util=r3a.pv_util, curtail=r3a.pv_curtail)],
        [(label="STATCOM $(Int(r.q_scale))×qd (52 units)", util=r.util,  curtail=100.0 - r.util)
         for r in rating_results]
    )
    print_delta_hc_table(table_rows)
end

# ───────────────────────────────────────────────────────────────
## ADD Case 4: HC curve across PV penetration levels
##
## Central thesis result. Shows ΔHC as a function of penetration:
## at which scale does the STATCOM start to matter, and how large
## is the gain at each level?
##
## Runs pairs of (no STATCOM, with STATCOM) solves at each scale.
# ───────────────────────────────────────────────────────────────
if do_case4
    println("\n" * "="^55)
    println(" CASE 4: HC Curve — STATCOM vs No-STATCOM")
    println(" PV scales 1–5×pd  |  bounds enforced")
    println(" STATCOM: $(BEST_Q_SCALE)×qd, 52 units, spacing=1")
    println("="^55)

    hc_curve_data = NamedTuple[]

    for scale in [1.0, 2.0, 3.0, 4.0, 5.0]
        println("\n  ── pv_scale = $(scale)×pd ──")

        dm_no = load_base_network(data_path, enforce_bounds=true)
        add_pv!(dm_no; pv_scale=scale, q_scale=1.25, spacing=1, pv_cost=PV_COST)
        r_no  = solve_and_report(dm_no, "PV $(scale)×  [no STATCOM]")

        dm_st = load_base_network(data_path, enforce_bounds=true)
        add_pv!(dm_st;       pv_scale=scale,        q_scale=1.25, spacing=1, pv_cost=PV_COST)
        add_statcoms!(dm_st; q_scale=BEST_Q_SCALE,  spacing=1,    statcom_cost=STATCOM_COST)
        r_st  = solve_and_report(dm_st, "PV $(scale)× + STATCOM $(BEST_Q_SCALE)×qd")

        push!(hc_curve_data, (scale=scale, util_no=r_no.pv_util, util_st=r_st.pv_util))
    end

    plot_hc_curve(hc_curve_data)

    println("\n  HC Curve Summary:")
    @printf("  %-8s  %-12s  %-14s  %s\n", "Scale", "No STATCOM", "With STATCOM", "ΔHC")
    println("  " * "─"^50)
    for d in hc_curve_data
        delta = d.util_st - d.util_no
        @printf("  %-8s  %-12s  %-14s  %s\n",
                "$(d.scale)×pd",
                @sprintf("%.1f%%", d.util_no),
                @sprintf("%.1f%%", d.util_st),
                @sprintf("+%.1f pp", delta))
    end
end

# ───────────────────────────────────────────────────────────────
## ADD Case 5: Targeted vs uniform placement
##
## Compares two strategies with the same device count and rating:
##   Uniform: one STATCOM every 10 load buses (~5–6 units)
##   Targeted: STATCOM at the 5 highest-voltage buses from Case 2
##
## If targeted outperforms uniform, it suggests the OPF-identified
## violation locations are also optimal placement locations.
# ───────────────────────────────────────────────────────────────
if do_case5
    println("\n" * "="^55)
    println(" CASE 5: Targeted vs Uniform Placement")
    println(" Both: $(length(WORST_BUSES)) units, $(BEST_Q_SCALE)×qd rating")
    println(" Targeted buses: $(WORST_BUSES)")
    println("="^55)

    println("\n  ── 5a: No STATCOM reference ──")
    dm5_ref = load_base_network(data_path, enforce_bounds=true)
    add_pv!(dm5_ref; pv_scale=STRESS_SCALE, q_scale=1.25, spacing=1, pv_cost=PV_COST)
    r5_ref  = solve_and_report(dm5_ref, "No STATCOM reference")

    println("\n  ── 5b: Uniform (spacing=10, ~$(length(WORST_BUSES)) units) ──")
    dm5_uni = load_base_network(data_path, enforce_bounds=true)
    add_pv!(dm5_uni;       pv_scale=STRESS_SCALE, q_scale=1.25,   spacing=1,  pv_cost=PV_COST)
    add_statcoms!(dm5_uni; q_scale=BEST_Q_SCALE,  spacing=10,     statcom_cost=STATCOM_COST)
    r5_uni  = solve_and_report(dm5_uni, "Uniform  ~$(length(WORST_BUSES)) units")

    println("\n  ── 5c: Targeted at worst-voltage buses ──")
    dm5_tgt = load_base_network(data_path, enforce_bounds=true)
    add_pv!(dm5_tgt; pv_scale=STRESS_SCALE, q_scale=1.25, spacing=1, pv_cost=PV_COST)
    add_statcoms!(dm5_tgt; q_scale=BEST_Q_SCALE, statcom_cost=STATCOM_COST,
                  target_buses=WORST_BUSES)
    r5_tgt  = solve_and_report(dm5_tgt, "Targeted  $(length(WORST_BUSES)) units")

    plot_placement_comparison(r5_ref.pv_util, r5_uni.pv_util, r5_tgt.pv_util)

    println("\n  Placement summary:")
    @printf("  %-30s  %.1f%%\n",                    "No STATCOM",         r5_ref.pv_util)
    @printf("  %-30s  %.1f%%  (Δ = +%.1f pp)\n",    "Uniform (~$(length(WORST_BUSES)) units)", r5_uni.pv_util, r5_uni.pv_util - r5_ref.pv_util)
    @printf("  %-30s  %.1f%%  (Δ = +%.1f pp)\n",    "Targeted ($(length(WORST_BUSES)) units)",  r5_tgt.pv_util, r5_tgt.pv_util - r5_ref.pv_util)
end

# ───────────────────────────────────────────────────────────────
## ADD Case 6: HC decomposition — inverter Q vs STATCOM Q
##
## Quantifies the separate contributions of smart inverter Q
## capability and STATCOM Q to the total ΔHC.
##
## The OPF already uses inverter Q (q_scale=1.25 gives headroom).
## This case makes that contribution explicit by comparing:
##   6a: Q locked to zero (pure active injection, no Q from inverter)
##   6b: Full inverter Q enabled (standard Case 3a scenario)
##   6c: Full inverter Q + STATCOM (full system)
# ───────────────────────────────────────────────────────────────
if do_case6
    println("\n" * "="^55)
    println(" CASE 6: HC Decomposition — Inverter Q vs STATCOM Q")
    println(" PV stress: $(STRESS_SCALE)×pd  |  bounds enforced")
    println("="^55)

    println("\n  ── 6a: PV active-only (Q disabled) ──")
    dm6a  = load_base_network(data_path, enforce_bounds=true)
    pv_6a = add_pv!(dm6a; pv_scale=STRESS_SCALE, q_scale=1.25, spacing=1, pv_cost=PV_COST)
    for id in pv_6a
        dm6a["gen"][id]["qmax"] = zeros(3)
        dm6a["gen"][id]["qmin"] = zeros(3)
    end
    r6a = solve_and_report(dm6a, "PV active-only  [Q = 0]")

    println("\n  ── 6b: Smart inverter only (Q enabled, no STATCOM) ──")
    dm6b = load_base_network(data_path, enforce_bounds=true)
    add_pv!(dm6b; pv_scale=STRESS_SCALE, q_scale=1.25, spacing=1, pv_cost=PV_COST)
    r6b  = solve_and_report(dm6b, "Smart inverter only  [no STATCOM]")

    println("\n  ── 6c: Smart inverter + STATCOM ($(BEST_Q_SCALE)×qd, 52 units) ──")
    dm6c = load_base_network(data_path, enforce_bounds=true)
    add_pv!(dm6c;       pv_scale=STRESS_SCALE, q_scale=1.25,   spacing=1, pv_cost=PV_COST)
    add_statcoms!(dm6c; q_scale=BEST_Q_SCALE,  spacing=1,      statcom_cost=STATCOM_COST)
    r6c  = solve_and_report(dm6c, "Smart inverter + STATCOM $(BEST_Q_SCALE)×qd")

    plot_decomposition(r6a.pv_util, r6b.pv_util, r6c.pv_util)

    println("\n  HC decomposition:")
    @printf("  PV (Q = 0, no inverter Q)       : %.1f%%\n", r6a.pv_util)
    @printf("  + Smart inverter Q              : %.1f%%  → +%.1f pp  (inverter contribution)\n",
            r6b.pv_util, r6b.pv_util - r6a.pv_util)
    @printf("  + STATCOM on top of inverter    : %.1f%%  → +%.1f pp  (STATCOM contribution)\n",
            r6c.pv_util, r6c.pv_util - r6b.pv_util)
    @printf("  Total ΔHC vs Q=0 baseline       : +%.1f pp\n",
            r6c.pv_util - r6a.pv_util)
end

println("\n  All plots saved to ./plots/")
println("  Done.")