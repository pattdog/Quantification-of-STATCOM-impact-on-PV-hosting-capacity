using Pkg
Pkg.activate("./")
using rosetta_distribution_opf
import PowerModelsDistribution
import InfrastructureModels
using Ipopt
using JuMP
using Statistics
using CairoMakie          # ADD: for all plots
using DataFrames          # ADD: for ΔHC summary tables
using Printf             # ADD THIS LINE
import Logging

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
# NETWORK LOADER  (unchanged)
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
# PV PLACEMENT  (unchanged)
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
# STATCOM PLACEMENT — EXTENDED: add targeted placement option
# ═══════════════════════════════════════════════════════════════
function add_statcoms!(data_math;
        q_scale      = 1.0,
        spacing      = 1,
        statcom_cost = 1.0,
        target_buses = nothing)   # ADD: pass specific bus IDs for targeted placement

    load_ids = sort(
        collect(keys(data_math["load"])),
        by = x -> tryparse(Int, x) === nothing ? 0 : parse(Int, x)
    )

    statcom_ids = String[]

    # ADD: targeted placement path — place at specific buses instead of uniformly
    if !isnothing(target_buses)
        qd_vals  = [load["qd"][1] for (i, load) in data_math["load"]]
        qd_mean  = mean(qd_vals)
        qlim     = q_scale * qd_mean

        for bus_id in target_buses
            gen_id = string(length(data_math["gen"]) + 1)
            data_math["gen"][gen_id] = deepcopy(data_math["gen"]["1"])
            gen = data_math["gen"][gen_id]

            gen["gen_bus"] = parse(Int, string(bus_id))
            gen["type"]    = "STATCOM"
            gen["name"]    = "statcom_targeted_bus$(bus_id)"
            gen["pmax"]    =  zeros(3)
            gen["pmin"]    =  zeros(3)
            gen["qmax"]    =  qlim * ones(3)
            gen["qmin"]    = -qlim * ones(3)
            gen["cost"]    = [statcom_cost 0.0]

            push!(statcom_ids, gen_id)
        end
        println("  STATCOM: $(length(statcom_ids)) units  TARGETED at buses $(target_buses)  q_scale=$(q_scale)×qd_mean  cost=$(statcom_cost)")
        return statcom_ids
    end

    # Original uniform path
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
# REPORT RESULTS — EXTENDED: return a Dict for plotting
# ═══════════════════════════════════════════════════════════════
function report_results(data_math, vr_val, vi_val, pg_val, qg_val, ref, status)

    # ADD: result dict — lets us collect data across runs without re-parsing stdout
    result = Dict{String, Any}(
        "status"        => string(status),
        "v_min"         => NaN, "v_mean" => NaN, "v_max" => NaN,
        "n_overvolt"    => 0,   "n_undervolt" => 0,
        "pv_util"       => NaN, "pv_curtail" => NaN,
        "pv_output"     => NaN, "pv_capacity" => NaN,
        "statcom_q"     => NaN, "statcom_q_cap" => NaN, "statcom_util" => NaN,
        # ADD: per-bus and per-phase data for voltage profile and Q heatmap plots
        "vm_per_bus"    => Dict{Int, Vector{Float64}}(),
        "qg_per_statcom"=> Dict{Int, Vector{Float64}}(),
    )

    if string(status) in ["LOCALLY_SOLVED", "OPTIMAL"]

        vm_all     = Float64[]
        bus_labels = Int[]
        vm_per_bus = Dict{Int, Vector{Float64}}()   # ADD: bus_id => [vA, vB, vC]

        for (b, bus) in ref[:bus]
            bus["bus_type"] == 3 && continue
            phases = Float64[]
            for p in 1:3
                v = abs(vr_val[p,b] + im * vi_val[p,b])
                push!(vm_all, v)
                push!(bus_labels, b)
                push!(phases, v)
            end
            vm_per_bus[b] = phases
        end

        result["vm_per_bus"] = vm_per_bus

        if !isempty(vm_all)
            v_min   = minimum(vm_all)
            v_mean  = mean(vm_all)
            v_max   = maximum(vm_all)
            n_over  = count(v -> v > 1.10, vm_all)
            n_under = count(v -> v < 0.90, vm_all)

            result["v_min"] = v_min;  result["v_mean"] = v_mean
            result["v_max"] = v_max;  result["n_overvolt"] = n_over
            result["n_undervolt"] = n_under

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
        end

        # PV dispatch
        pv_gen_ids = [parse(Int, i) for (i, g) in data_math["gen"] if get(g, "type", "") == "PV"]
        if !isempty(pv_gen_ids)
            pg_total = sum(sum(pg_val[:, i]) for i in pv_gen_ids)
            pg_cap   = sum(sum(data_math["gen"][string(i)]["pmax"]) for i in pv_gen_ids)
            util     = pg_total / max(1e-9, pg_cap) * 100
            curtail  = max(0.0, 100.0 - util)
            result["pv_output"] = pg_total;  result["pv_capacity"] = pg_cap
            result["pv_util"] = util;         result["pv_curtail"] = curtail
            println("    PV output / capacity  : $(round(pg_total,sigdigits=4)) / $(round(pg_cap,sigdigits=4)) pu")
            println("    PV utilisation        : $(round(util,digits=1))%  →  curtailment: $(round(curtail,digits=1))%")
        end

        # STATCOM dispatch — ADD: per-STATCOM per-phase Q for heatmap
        st_gen_ids = [parse(Int, i) for (i, g) in data_math["gen"] if get(g, "type", "") == "STATCOM"]
        if !isempty(st_gen_ids)
            qg_total  = sum(sum(qg_val[:, i]) for i in st_gen_ids)
            qg_cap    = sum(sum(data_math["gen"][string(i)]["qmax"]) for i in st_gen_ids)
            util      = abs(qg_total) / max(1e-9, qg_cap) * 100
            direction = qg_total >= 0 ? "injecting ↑V" : "absorbing ↓V"

            # ADD: store per-STATCOM Q and bus for spatial plot
            for sid in st_gen_ids
                bus_id = data_math["gen"][string(sid)]["gen_bus"]
                result["qg_per_statcom"][bus_id] = qg_val[:, sid]
            end

            result["statcom_q"] = qg_total;   result["statcom_q_cap"] = qg_cap
            result["statcom_util"] = util
            println("    STATCOM Q / capacity  : $(round(qg_total,sigdigits=4)) / $(round(qg_cap,sigdigits=4)) pu")
            println("    STATCOM utilisation   : $(round(util,digits=1))%  $(direction)")

            # ADD: per-phase STATCOM dispatch breakdown
            phase_labels = ["Phase A", "Phase B", "Phase C"]
            for (pi, ph) in enumerate(phase_labels)
                q_ph = sum(qg_val[pi, sid] for sid in st_gen_ids)
                println("      $(ph) net Q: $(round(q_ph, sigdigits=3)) pu  $(q_ph >= 0 ? "↑" : "↓")")
            end
        end

    else
        println("    WARNING: $(status)")
    end

    return result   # ADD: return for collection
end

# ═══════════════════════════════════════════════════════════════
# SOLVE AND REPORT — unchanged except returns result dict
# ═══════════════════════════════════════════════════════════════
function solve_and_report(data_math, label; obj="cost")

    global ref = IM.build_ref(
        data_math,
        PMD.ref_add_core!,
        PMD._pmd_global_keys,
        PMD.pmd_it_name
    )[:it][:pmd][:nw][0]

    global model = JuMP.Model(ipopt_solver)
    global objective = obj

    _old_logger = Logging.global_logger()
    Logging.global_logger(Logging.ConsoleLogger(stderr, Logging.Error))
    _old_stdout = stdout
    (rd, wr) = redirect_stdout()

    include("./core/variables.jl")
    include("./core/constraints.jl")
    include("./core/objectives.jl")

    redirect_stdout(_old_stdout)
    close(wr)
    Logging.global_logger(_old_logger)

    JuMP.optimize!(model)
    status = JuMP.termination_status(model)
    println("\n  [$label]  →  $(status)")

    vr_val = JuMP.value.(vr)
    vi_val = JuMP.value.(vi)
    pg_val = JuMP.value.(pg)
    qg_val = JuMP.value.(qg)

    result = report_results(data_math, vr_val, vi_val, pg_val, qg_val, ref, status)
    return result  # ADD: return dict
end

# ═══════════════════════════════════════════════════════════════
# ADD: PLOTTING FUNCTIONS
# All plots saved to ./plots/
# ═══════════════════════════════════════════════════════════════
mkpath("./plots")

# ── Plot 1: Voltage profile — sorted bus voltages, with/without STATCOM ────
# Call after you have two result dicts: baseline_r and statcom_r
function plot_voltage_profile(baseline_r, statcom_r, pv_scale; savepath="./plots/voltage_profile.pdf")
    fig = Figure(size=(900, 480), backgroundcolor=:white)
    ax  = Axis(fig[1,1],
        xlabel = "Bus index (sorted by baseline voltage)",
        ylabel = "Voltage magnitude (pu)",
        title  = "Voltage profile — PV $(pv_scale)× load, with and without STATCOM",
        titlesize = 14,
        xlabelsize = 12, ylabelsize = 12,
        ygridvisible = true, xgridvisible = false,
        ygridcolor = (:black, 0.07),
    )

    # Extract max voltage per bus (worst phase) for each result
    function sorted_vmax(result_dict)
        bvs = result_dict["vm_per_bus"]
        isempty(bvs) && return Int[], Float64[]
        bus_ids  = sort(collect(keys(bvs)))
        v_maxes  = [maximum(bvs[b]) for b in bus_ids]
        order    = sortperm(v_maxes)   # sort by baseline voltage ascending
        return bus_ids[order], v_maxes[order]
    end

    b_ids, v_base = sorted_vmax(baseline_r)
    _, v_stat     = sorted_vmax(statcom_r)

    xs = 1:length(b_ids)

    # Shaded violation zone
    hspan!(ax, 1.10, 1.20, color=(:red, 0.06))
    hspan!(ax, 0.80, 0.90, color=(:blue, 0.06))

    lines!(ax, xs, v_base, color=(:steelblue, 0.8), linewidth=2.0, label="No STATCOM")
    lines!(ax, xs, v_stat, color=(:darkorange, 0.9), linewidth=2.0, label="With STATCOM (optimal Q dispatch)")

    # Limit lines
    hlines!(ax, [1.10], color=:red,  linestyle=:dash, linewidth=1.5, label="1.10 pu upper limit")
    hlines!(ax, [0.90], color=:blue, linestyle=:dash, linewidth=1.5, label="0.90 pu lower limit")
    hlines!(ax, [1.00], color=(:black, 0.2), linestyle=:dot, linewidth=1.0)

    axislegend(ax, position=:lt, framevisible=true, labelsize=11)
    ylims!(ax, 0.92, 1.15)

    save(savepath, fig)
    println("  → Saved: $savepath")
    return fig
end

# ── Plot 2: HC curve — with and without STATCOM across PV penetration ───────
# hc_data format: Vector of NamedTuples (scale, util_no_st, util_with_st)
function plot_hc_curve(hc_data; savepath="./plots/hc_curve.pdf")
    fig = Figure(size=(800, 460), backgroundcolor=:white)
    ax  = Axis(fig[1,1],
        xlabel = "PV penetration scale (× load demand)",
        ylabel = "PV utilisation (%)",
        title  = "Hosting Capacity: with and without STATCOM",
        titlesize = 14,
        xlabelsize = 12, ylabelsize = 12,
        ygridvisible = true, xgridvisible = false,
        ygridcolor   = (:black, 0.07),
    )

    scales    = [d.scale    for d in hc_data]
    util_no   = [d.util_no  for d in hc_data]
    util_st   = [d.util_st  for d in hc_data]
    delta_hc  = util_st .- util_no

    # Fill between curves — the ΔHC region
    band!(ax, scales, util_no, util_st, color=(:darkorange, 0.15))

    lines!(ax, scales, util_no, color=:steelblue,   linewidth=2.5, label="No STATCOM")
    lines!(ax, scales, util_st, color=:darkorange,  linewidth=2.5, label="With STATCOM (100×qd, 52 units)")

    scatter!(ax, scales, util_no, color=:steelblue,  markersize=8)
    scatter!(ax, scales, util_st, color=:darkorange, markersize=8)

    # Annotate ΔHC at each point
    for (s, d, y) in zip(scales, delta_hc, util_st)
        d > 0.5 || continue
        text!(ax, s, y + 1.5, text="+$(round(d,digits=1))pp", fontsize=9,
              align=(:center,:bottom), color=:darkorange)
    end

    hlines!(ax, [100.0], color=(:green, 0.5), linestyle=:dash, linewidth=1.5, label="100% utilisation")

    ylims!(ax, 60, 108)
    xlims!(ax, 0.8, 5.5)
    axislegend(ax, position=:lb, framevisible=true, labelsize=11)

    save(savepath, fig)
    println("  → Saved: $savepath")
    return fig
end

# ── Plot 3: Rating sweep — PV utilisation vs STATCOM Q rating ───────────────
function plot_rating_sweep(rating_data; savepath="./plots/rating_sweep.pdf")
    fig = Figure(size=(800, 460), backgroundcolor=:white)
    ax  = Axis(fig[1,1],
        xlabel = "STATCOM rating (multiples of q_d per unit)",
        ylabel = "PV utilisation (%)",
        title  = "PV hosting capacity vs. STATCOM reactive power rating",
        titlesize = 14,
        xlabelsize = 12, ylabelsize = 12,
        xscale = log10,
        ygridvisible = true, xgridvisible = false,
        ygridcolor   = (:black, 0.07),
    )

    ratings  = [d.q_scale for d in rating_data]
    utils    = [d.util    for d in rating_data]
    baseline = rating_data[1].baseline  # no-STATCOM reference

    # Baseline reference band
    hlines!(ax, [baseline], color=:red, linestyle=:dash, linewidth=1.5,
            label="No STATCOM baseline ($(round(baseline,digits=1))%)")
    hlines!(ax, [100.0],    color=(:green,0.5), linestyle=:dash, linewidth=1.5,
            label="100% utilisation")

    band!(ax, ratings, fill(baseline, length(ratings)), utils, color=(:darkorange, 0.12))
    lines!(ax, ratings, utils, color=:darkorange, linewidth=2.5)
    scatter!(ax, ratings, utils, color=:darkorange, markersize=9, label="STATCOM (52 units)")

    # Annotate Q→V sign flip (absorb → inject)
    # From results: flip occurs between q100 and q200
    vlines!(ax, [100.0], color=(:purple, 0.4), linestyle=:dot, linewidth=1.2)
    text!(ax, 105.0, 78.0,
          text="Q mode flips\nabsorb→inject\nbeyond 100×", fontsize=9,
          color=:purple, align=(:left,:bottom))

    axislegend(ax, position=:lt, framevisible=true, labelsize=11)
    ylims!(ax, 65, 108)

    save(savepath, fig)
    println("  → Saved: $savepath")
    return fig
end

# ── Plot 4: Density sweep — utilisation vs number of STATCOMs ───────────────
function plot_density_sweep(density_data; savepath="./plots/density_sweep.pdf")
    fig = Figure(size=(800, 460), backgroundcolor=:white)
    ax  = Axis(fig[1,1],
        xlabel = "Number of STATCOM units deployed",
        ylabel = "PV utilisation (%)",
        title  = "Hosting capacity vs. STATCOM deployment density (rating = 100×q_d)",
        titlesize = 14,
        xlabelsize = 12, ylabelsize = 12,
        ygridvisible = true, xgridvisible = false,
        ygridcolor   = (:black, 0.07),
    )

    n_units  = [d.n_units  for d in density_data]
    utils    = [d.util     for d in density_data]
    st_utils = [d.st_util  for d in density_data]  # STATCOM Q utilisation %
    baseline = density_data[1].baseline

    hlines!(ax, [baseline], color=:red, linestyle=:dash, linewidth=1.5,
            label="No STATCOM ($(round(baseline,digits=1))%)")
    hlines!(ax, [100.0], color=(:green,0.5), linestyle=:dash, linewidth=1.5)

    band!(ax, n_units, fill(baseline, length(n_units)), utils, color=(:steelblue,0.12))
    lines!(ax, n_units, utils, color=:steelblue, linewidth=2.5)
    scatter!(ax, n_units, utils, color=:steelblue, markersize=9, label="PV utilisation")

    # Overlay STATCOM utilisation on secondary axis
    ax2 = Axis(fig[1,1],
        yaxisposition = :right,
        ylabel = "STATCOM Q utilisation (%)",
        ylabelsize = 11,
        yticklabelcolor = :purple,
        ylabelcolor = :purple,
        ygridvisible = false,
    )
    hidexdecorations!(ax2)
    lines!(ax2, n_units, st_utils, color=(:purple, 0.7), linewidth=1.5,
           linestyle=:dot, label="STATCOM Q utilisation")
    scatter!(ax2, n_units, st_utils, color=:purple, markersize=7, marker=:diamond)

    axislegend(ax,  position=:rb, framevisible=true, labelsize=10)
    axislegend(ax2, position=:rt, framevisible=true, labelsize=10)

    save(savepath, fig)
    println("  → Saved: $savepath")
    return fig
end

# ── Plot 5: Targeted vs uniform placement comparison bar ─────────────────────
function plot_placement_comparison(uniform_util, targeted_util, baseline_util;
        savepath="./plots/placement_comparison.pdf")
    fig = Figure(size=(600, 420), backgroundcolor=:white)
    ax  = Axis(fig[1,1],
        ylabel    = "PV utilisation (%)",
        title     = "Placement strategy comparison (same STATCOM count and rating)",
        titlesize = 14,
        ylabelsize = 12,
        ygridvisible = true, xgridvisible = false,
        xticks = (1:3, ["No STATCOM\n(baseline)", "Uniform placement\n(co-located with PV)", "Targeted placement\n(worst-voltage buses)"]),
        xticklabelsize = 11,
    )

    vals = [baseline_util, uniform_util, targeted_util]
    colors = [:steelblue, :darkorange, :darkgreen]
    barplot!(ax, 1:3, vals, color=colors, width=0.55)

    # Value labels
    for (i, v) in enumerate(vals)
        text!(ax, i, v + 0.5, text="$(round(v,digits=1))%",
              align=(:center,:bottom), fontsize=12, color=colors[i])
    end

    hlines!(ax, [100.0], color=(:green,0.4), linestyle=:dash, linewidth=1.5)
    ylims!(ax, 60, 108)

    save(savepath, fig)
    println("  → Saved: $savepath")
    return fig
end

# ── Plot 6: ΔHC summary table ────────────────────────────────────────────────
function print_delta_hc_table(results_vec)
    println("\n  ΔHC Summary Table")
    println("  " * "─"^72)
    println("  $(lpad("Configuration",30))  $(lpad("Utilisation",12))  $(lpad("Curtailment",12))  $(lpad("ΔHC vs baseline",16))")
    println("  " * "─"^72)
    baseline = results_vec[1].util
    for r in results_vec
        delta = r.util - baseline
        delta_str = delta == 0.0 ? "  —" : @sprintf("+%.1fpp", delta)
        println("  $(lpad(r.label,30))  $(lpad(@sprintf("%.1f%%",r.util),12))  $(lpad(@sprintf("%.1f%%",r.curtail),12))  $(lpad(delta_str,16))")
    end
    println("  " * "─"^72)
end

# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════
using Printf

do_case1 = true
do_case2 = true
do_case3 = true
do_case4 = true   # ADD: HC curve across PV penetration with STATCOM
do_case5 = true   # ADD: targeted placement vs uniform
do_case6 = true   # ADD: smart inverter only baseline

const STRESS_SCALE  = 5.0
const PV_COST       = -1000.0
const STATCOM_COST  =  1.0
const BEST_Q_SCALE  = 100      # from Case 3b: practical saturation point

# Worst-voltage buses from Case 2 results (used for targeted placement)
const WORST_BUSES   = [189, 256, 157, 843, 260]

dm_ref = load_base_network(data_path)
summarise_network(dm_ref)

# ───────────────────────────────────────────────────────────────
## Case 1: Natural baseline  (unchanged)
# ───────────────────────────────────────────────────────────────
if do_case1
    println("\n" * "="^55)
    println(" CASE 1: Natural Baseline (no bounds, no DER)")
    println("="^55)
    dm1 = load_base_network(data_path, enforce_bounds=false)
    r1  = solve_and_report(dm1, "Baseline"; obj="cost")
end

# ───────────────────────────────────────────────────────────────
## Case 2: PV penetration sweep — no bounds  (unchanged)
# ───────────────────────────────────────────────────────────────
if do_case2
    println("\n" * "="^55)
    println(" CASE 2: PV Penetration Sweep (no bounds)")
    println("="^55)

    for scale in [1.0, 2.0, 3.0, 4.0, 5.0]
        println("\n  ── pv_scale = $(scale)×pd ──")
        dm = load_base_network(data_path, enforce_bounds=false)
        add_pv!(dm; pv_scale=scale, q_scale=1.25, spacing=1, pv_cost=-1000.0)
        solve_and_report(dm, "PV scale=$(scale)×pd"; obj="cost")
    end
end

# ───────────────────────────────────────────────────────────────
## Case 3: STATCOM hosting capacity — bounds enforced  (unchanged)
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
    r3a  = solve_and_report(dm3a, "PV only  [no STATCOM]"; obj="cost")

    println("\n  ── 3b: STATCOM rating sweep (52 STATCOMs, spacing=1) ──")
    rating_results = NamedTuple[]
    for q_sc in [1, 2, 5, 10, 20, 50, 100, 200, 500]
        println("\n  ── q_scale = $(q_sc)×qd ──")
        local dm3b = load_base_network(data_path, enforce_bounds=true)
        add_pv!(dm3b;       pv_scale=STRESS_SCALE, q_scale=1.25, spacing=1, pv_cost=PV_COST)
        add_statcoms!(dm3b; q_scale=q_sc,          spacing=1,    statcom_cost=STATCOM_COST)
        r = solve_and_report(dm3b, "PV+STATCOM q=$(q_sc)×qd  52 units"; obj="cost")
        push!(rating_results, (q_scale=Float64(q_sc), util=r["pv_util"], baseline=r3a["pv_util"]))
    end
    # ADD: produce rating sweep plot
    plot_rating_sweep(rating_results)

    println("\n  ── 3c: Density sweep at q_scale=$(BEST_Q_SCALE)×qd ──")
    density_results = NamedTuple[]
    r3c_52 = nothing
    for sp in [1, 2, 4, 6, 8, 10]
        n_st = ceil(Int, 52 / sp)
        println("\n  ── spacing=$(sp) → ~$(n_st) STATCOMs ──")
        local dm3c = load_base_network(data_path, enforce_bounds=true)
        add_pv!(dm3c;       pv_scale=STRESS_SCALE, q_scale=1.25,    spacing=1, pv_cost=PV_COST)
        add_statcoms!(dm3c; q_scale=BEST_Q_SCALE,  spacing=sp,      statcom_cost=STATCOM_COST)
        r = solve_and_report(dm3c, "PV+STATCOM q=$(BEST_Q_SCALE)×qd  $(n_st) units"; obj="cost")
        push!(density_results, (n_units=n_st, util=r["pv_util"], st_util=r["statcom_util"], baseline=r3a["pv_util"]))
        sp == 1 && (r3c_52 = r)
    end
    # ADD: produce density sweep plot + voltage profile comparison
    plot_density_sweep(density_results)
    if !isnothing(r3c_52)
        plot_voltage_profile(r3a, r3c_52, STRESS_SCALE)
    end

    # ADD: ΔHC summary table
    table_rows = vcat(
        [(label="No STATCOM (baseline)", util=r3a["pv_util"], curtail=r3a["pv_curtail"])],
        [(label="STATCOM $(r.q_scale)×qd (52 units)", util=r.util, curtail=100.0-r.util)
         for r in rating_results]
    )
    print_delta_hc_table(table_rows)
end

# ───────────────────────────────────────────────────────────────
## ADD Case 4: HC curve — with and without STATCOM across PV scales
# This is the CENTRAL thesis result — ΔHC as a function of penetration
# ───────────────────────────────────────────────────────────────
if do_case4
    println("\n" * "="^55)
    println(" CASE 4: HC Curve — STATCOM vs No-STATCOM")
    println(" PV scales 1–5×pd  |  bounds enforced  |  STATCOM: $(BEST_Q_SCALE)×qd  52 units")
    println("="^55)

    hc_curve_data = NamedTuple[]
    for scale in [1.0, 2.0, 3.0, 4.0, 5.0]
        println("\n  ── pv_scale = $(scale)×pd ──")

        # Without STATCOM
        dm_no = load_base_network(data_path, enforce_bounds=true)
        add_pv!(dm_no; pv_scale=scale, q_scale=1.25, spacing=1, pv_cost=PV_COST)
        r_no  = solve_and_report(dm_no, "PV $(scale)×  [no STATCOM]"; obj="cost")

        # With STATCOM (best configuration from Case 3b)
        dm_st = load_base_network(data_path, enforce_bounds=true)
        add_pv!(dm_st;       pv_scale=scale,       q_scale=1.25,     spacing=1, pv_cost=PV_COST)
        add_statcoms!(dm_st; q_scale=BEST_Q_SCALE,  spacing=1,        statcom_cost=STATCOM_COST)
        r_st  = solve_and_report(dm_st, "PV $(scale)× + STATCOM $(BEST_Q_SCALE)×qd"; obj="cost")

        push!(hc_curve_data, (
            scale    = scale,
            util_no  = r_no["pv_util"],
            util_st  = r_st["pv_util"],
        ))
    end
    plot_hc_curve(hc_curve_data)

    println("\n  HC Curve Summary:")
    println("  $(lpad("Scale",8))  $(lpad("No STATCOM",12))  $(lpad("With STATCOM",14))  $(lpad("ΔHC",8))")
    println("  " * "─"^50)
    for d in hc_curve_data
        delta = d.util_st - d.util_no
        println("  $(lpad("$(d.scale)×pd",8))  $(lpad(@sprintf("%.1f%%",d.util_no),12))  $(lpad(@sprintf("%.1f%%",d.util_st),14))  $(lpad(@sprintf("+%.1fpp",delta),8))")
    end
end

# ───────────────────────────────────────────────────────────────
## ADD Case 5: Targeted placement — STATCOM at worst-voltage buses
# Tests whether placing fewer STATCOMs at the right locations
# beats uniform co-location
# ───────────────────────────────────────────────────────────────
if do_case5
    println("\n" * "="^55)
    println(" CASE 5: Targeted vs Uniform Placement")
    println(" STATCOM at worst-voltage buses: $(WORST_BUSES)")
    println(" Same count (5 units), same rating $(BEST_Q_SCALE)×qd_mean")
    println("="^55)

    # Uniform: 5 STATCOMs, spacing=10 (matches ~same count)
    println("\n  ── 5a: Uniform placement (spacing=10, ~6 units) ──")
    dm5a = load_base_network(data_path, enforce_bounds=true)
    add_pv!(dm5a;       pv_scale=STRESS_SCALE, q_scale=1.25,    spacing=1,  pv_cost=PV_COST)
    add_statcoms!(dm5a; q_scale=BEST_Q_SCALE,  spacing=10,      statcom_cost=STATCOM_COST)
    r5a = solve_and_report(dm5a, "Uniform placement  ~6 units"; obj="cost")

    # Targeted: STATCOM placed at the 5 highest-voltage buses from Case 2
    println("\n  ── 5b: Targeted placement at worst-voltage buses ──")
    dm5b = load_base_network(data_path, enforce_bounds=true)
    add_pv!(dm5b; pv_scale=STRESS_SCALE, q_scale=1.25, spacing=1, pv_cost=PV_COST)
    add_statcoms!(dm5b; q_scale=BEST_Q_SCALE, statcom_cost=STATCOM_COST,
                  target_buses=WORST_BUSES)
    r5b = solve_and_report(dm5b, "Targeted placement  5 units"; obj="cost")

    # No STATCOM reference
    dm5c = load_base_network(data_path, enforce_bounds=true)
    add_pv!(dm5c; pv_scale=STRESS_SCALE, q_scale=1.25, spacing=1, pv_cost=PV_COST)
    r5c = solve_and_report(dm5c, "No STATCOM reference"; obj="cost")

    plot_placement_comparison(r5a["pv_util"], r5b["pv_util"], r5c["pv_util"])

    println("\n  Placement comparison:")
    println("    No STATCOM         : $(round(r5c["pv_util"],digits=1))% utilisation")
    println("    Uniform (~6 units) : $(round(r5a["pv_util"],digits=1))% utilisation  Δ=+$(round(r5a["pv_util"]-r5c["pv_util"],digits=1))pp")
    println("    Targeted (5 units) : $(round(r5b["pv_util"],digits=1))% utilisation  Δ=+$(round(r5b["pv_util"]-r5c["pv_util"],digits=1))pp")
end

# ───────────────────────────────────────────────────────────────
## ADD Case 6: Smart inverter only — no STATCOM but PV Q enabled
# Isolates how much work the PV inverters do vs the STATCOM
# PV units already have q_scale=1.25 giving reactive headroom.
# This case lets the OPF use that headroom without a STATCOM.
# The result shows how much ΔHC is inverter-driven vs STATCOM-driven.
# ───────────────────────────────────────────────────────────────
if do_case6
    println("\n" * "="^55)
    println(" CASE 6: Smart Inverter Only Baseline")
    println(" PV with Q capability (q_scale=1.25), no STATCOM")
    println(" Compare Q-constrained (pmin=qmin=0) vs full inverter capability")
    println("="^55)

    # 6a: PV with Q locked to zero (pure active power, no reactive)
    println("\n  ── 6a: PV active-only (Q capability disabled) ──")
    dm6a = load_base_network(data_path, enforce_bounds=true)
    pv_ids_6a = add_pv!(dm6a; pv_scale=STRESS_SCALE, q_scale=1.25, spacing=1, pv_cost=PV_COST)
    for id in pv_ids_6a
        dm6a["gen"][id]["qmax"] = zeros(3)
        dm6a["gen"][id]["qmin"] = zeros(3)
    end
    r6a = solve_and_report(dm6a, "PV active-only  [Q=0]"; obj="cost")

    # 6b: PV with full inverter Q capability (this is already Case 3a,
    #     but we run it again explicitly for the comparison)
    println("\n  ── 6b: PV with smart inverter Q (q_scale=1.25, no STATCOM) ──")
    dm6b = load_base_network(data_path, enforce_bounds=true)
    add_pv!(dm6b; pv_scale=STRESS_SCALE, q_scale=1.25, spacing=1, pv_cost=PV_COST)
    r6b = solve_and_report(dm6b, "Smart inverter only  [no STATCOM]"; obj="cost")

    # 6c: Smart inverter + STATCOM — full system
    println("\n  ── 6c: Smart inverter + STATCOM (full system) ──")
    dm6c = load_base_network(data_path, enforce_bounds=true)
    add_pv!(dm6c;       pv_scale=STRESS_SCALE, q_scale=1.25,   spacing=1, pv_cost=PV_COST)
    add_statcoms!(dm6c; q_scale=BEST_Q_SCALE,  spacing=1,      statcom_cost=STATCOM_COST)
    r6c = solve_and_report(dm6c, "Smart inverter + STATCOM $(BEST_Q_SCALE)×qd"; obj="cost")

    println("\n  Decomposition of HC contributions:")
    println("    PV (Q=0, no inverter Q)       : $(round(r6a["pv_util"],digits=1))% utilisation")
    println("    + Smart inverter Q             : $(round(r6b["pv_util"],digits=1))%  Δ=+$(round(r6b["pv_util"]-r6a["pv_util"],digits=1))pp  ← inverter contribution")
    println("    + STATCOM on top of inverter   : $(round(r6c["pv_util"],digits=1))%  Δ=+$(round(r6c["pv_util"]-r6b["pv_util"],digits=1))pp  ← STATCOM contribution")
    println("    Total ΔHC vs Q=0 baseline      : +$(round(r6c["pv_util"]-r6a["pv_util"],digits=1))pp")
end

println("\n  All plots saved to ./plots/")
println("  Done.")