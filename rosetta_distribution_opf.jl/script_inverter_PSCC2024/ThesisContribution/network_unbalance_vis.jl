#=
==============================================================================
NETWORK UNBALANCE VISUALISATION
==============================================================================
Purpose: visualise the unbalance already present in the network dataset,
BEFORE any PV or STATCOM is added. This isolates the baseline problem your
thesis is responding to -- how unbalanced is customer/load allocation across
the three phases on this feeder, and what does that do to voltage unbalance
under natural operating conditions.

Built on the same validated pipeline as the rest of the thesis scripts (raw
includes of variables.jl / constraints_PBalance.jl / objectives_FIXED.jl,
sbase via make_pu=false + make_per_unit!, add_start_vrvi!, objective="cost").
No PV, no STATCOM -- just the feeder as given, solved once at
load_multiplier = 1.0 with enforce_bounds = false so nothing is artificially
masked by voltage constraints.

Produces four plots, saved to ./plots/:
  load_phase_distribution.pdf  -- count of loads per phase, and total P/Q
                                   demand per phase (the root cause)
  bus_vuf_sorted.pdf           -- VUF at every bus, sorted ascending, with
                                   the ENWL 2% (typical grid-code-style)
                                   reference line
  bus_vuf_histogram.pdf        -- distribution of VUF across all buses
  voltage_profile_by_phase.pdf -- three-phase voltage magnitude profile
                                   along the feeder (sorted by mean voltage),
                                   showing how far the phases diverge from
                                   each other under natural imbalance

Run with:
    julia --project=. network_unbalance_visualisation.jl
==============================================================================
=#

using Logging
Logging.disable_logging(Logging.Warn)

using Pkg
Pkg.activate("./")
using rosetta_distribution_opf
import PowerModelsDistribution
import InfrastructureModels
using Ipopt
using JuMP
using Printf
using Statistics
using CairoMakie

const PMD  = PowerModelsDistribution
const RPMD = rosetta_distribution_opf
const IM   = InfrastructureModels
PMD.silence!()

ipopt_solver = JuMP.optimizer_with_attributes(
    Ipopt.Optimizer,
    "print_level" => 0,
    "sb"          => "yes",
    "max_iter"    => 8000,
    "warm_start_init_point" => "no",
)

data_path = "./rosetta_distribution_opf.jl/data/ENWL_4w_Network1_Feeder1/Master.dss"

kw_to_pu(p_kw, sbase_kva)     = p_kw / sbase_kva
kvar_to_pu(q_kvar, sbase_kva) = q_kvar / sbase_kva
pu_to_kw(p_pu, sbase_kva)     = p_pu * sbase_kva
pu_to_kvar(q_pu, sbase_kva)   = q_pu * sbase_kva

# -----------------------------------------------------------------------
# Network loader -- same validated sbase fix as the rest of the thesis
# scripts. enforce_bounds=false so voltage isn't artificially clamped --
# we want to see the network's natural, unmasked behaviour.
# -----------------------------------------------------------------------
function load_base_network(data_path; load_multiplier=1.0, enforce_bounds=false, sbase_kva=100)
    data_eng  = PMD.parse_file(data_path, transformations=[PMD.transform_loops!])
    data_math = PMD.transform_data_model(
        data_eng, multinetwork=false, kron_reduce=false, phase_project=false,
        make_pu=false
    )
    PMD.make_per_unit!(data_math; sbase=sbase_kva)
    PMD.add_start_vrvi!(data_math)

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

function get_sbase_kva(data_math)
    return data_math["settings"]["sbase"] * data_math["settings"]["power_scale_factor"] / 1000
end

# -----------------------------------------------------------------------
# Load-side imbalance: counts and total demand per phase. This is the ROOT
# CAUSE the rest of the plots trace through to a voltage-side symptom.
# -----------------------------------------------------------------------
function load_phase_stats(data_math, sbase_kva)
    tally   = Dict(1 => 0, 2 => 0, 3 => 0)
    pd_kw   = Dict(1 => 0.0, 2 => 0.0, 3 => 0.0)
    qd_kvar = Dict(1 => 0.0, 2 => 0.0, 3 => 0.0)
    threephase = 0

    for (id, load) in data_math["load"]
        phase_conns = filter(c -> c != 4, load["connections"])
        if length(phase_conns) == 1
            p = phase_conns[1]
            tally[p] += 1
            pd_kw[p]   += pu_to_kw(sum(load["pd"]), sbase_kva)
            qd_kvar[p] += pu_to_kvar(sum(load["qd"]), sbase_kva)
        else
            threephase += 1
        end
    end

    return (tally=tally, pd_kw=pd_kw, qd_kvar=qd_kvar, threephase=threephase)
end

# -----------------------------------------------------------------------
# Build + solve -- raw includes, objective="cost" (no PV/STATCOM present,
# so this just resolves the natural power flow via the substation gen).
# -----------------------------------------------------------------------
function build_and_solve(data_math)
    global ref = IM.build_ref(data_math, PMD.ref_add_core!, PMD._pmd_global_keys, PMD.pmd_it_name)[:it][:pmd][:nw][0]
    global model = JuMP.Model(ipopt_solver)

    include("../core/variables.jl")
    include("../core/constraints_unified.jl")

    global objective = "cost"
    include("../core/objectives_FIXED.jl")
    JuMP.optimize!(model)
    return JuMP.termination_status(model)
end

# -----------------------------------------------------------------------
# Per-bus VUF and voltage magnitude, for every bus with all 3 phases.
# -----------------------------------------------------------------------
function all_bus_vuf_and_vmag()
    alpha = exp(im*2/3*pi)
    T = 1/3 * [1 1 1 ; 1 alpha alpha^2 ; 1 alpha^2 alpha]

    bus_ids   = Int[]
    vufs      = Float64[]
    vm_ph     = Vector{Vector{Float64}}()   # [bus][phase] magnitude

    for (i, bus) in ref[:bus]
        get(bus, "bus_type", 1) == 3 && continue   # skip source/reference bus
        terms = bus["terminals"]
        !(1 in terms && 2 in terms && 3 in terms) && continue

        vr_val = [JuMP.value(vr[p,i]) for p in 1:3]
        vi_val = [JuMP.value(vi[p,i]) for p in 1:3]
        vph    = vr_val .+ im .* vi_val
        v012   = T * vph
        vpos, vneg = abs(v012[2]), abs(v012[3])
        vuf = vpos > 1e-6 ? vneg/vpos : 0.0

        push!(bus_ids, i)
        push!(vufs, vuf)
        push!(vm_ph, abs.(vph))
    end

    return bus_ids, vufs, vm_ph
end

# -----------------------------------------------------------------------
# PLOTS
# -----------------------------------------------------------------------
mkpath("./plots")

function plot_load_phase_distribution(stats; savepath="./plots/load_phase_distribution.pdf")
    labels = ["Phase A", "Phase B", "Phase C"]
    counts = [stats.tally[1], stats.tally[2], stats.tally[3]]
    pd     = [stats.pd_kw[1], stats.pd_kw[2], stats.pd_kw[3]]
    qd     = [stats.qd_kvar[1], stats.qd_kvar[2], stats.qd_kvar[3]]

    fig = Figure(size=(1000, 420), backgroundcolor=:white)

    ax1 = Axis(fig[1,1],
        title="Load count per phase", ylabel="Number of loads",
        xticks=(1:3, labels), titlesize=14,
        ygridvisible=true, xgridvisible=false, ygridcolor=(:black,0.08))
    barplot!(ax1, 1:3, counts, color=[:steelblue, :darkorange, :forestgreen], width=0.5)
    for (i,c) in enumerate(counts)
        text!(ax1, i, c+0.5, text=string(c), align=(:center,:bottom), fontsize=12)
    end

    ax2 = Axis(fig[1,2],
        title="Total active demand per phase", ylabel="kW",
        xticks=(1:3, labels), titlesize=14,
        ygridvisible=true, xgridvisible=false, ygridcolor=(:black,0.08))
    barplot!(ax2, 1:3, pd, color=[:steelblue, :darkorange, :forestgreen], width=0.5)
    for (i,v) in enumerate(pd)
        text!(ax2, i, v+maximum(pd)*0.02, text=@sprintf("%.1f",v), align=(:center,:bottom), fontsize=11)
    end

    ax3 = Axis(fig[1,3],
        title="Total reactive demand per phase", ylabel="kVAr",
        xticks=(1:3, labels), titlesize=14,
        ygridvisible=true, xgridvisible=false, ygridcolor=(:black,0.08))
    barplot!(ax3, 1:3, qd, color=[:steelblue, :darkorange, :forestgreen], width=0.5)
    for (i,v) in enumerate(qd)
        text!(ax3, i, v+maximum(qd)*0.02, text=@sprintf("%.1f",v), align=(:center,:bottom), fontsize=11)
    end

    Label(fig[0,1:3], "Load-side (root cause) phase imbalance -- $(stats.threephase) three-phase loads excluded",
          fontsize=15, font=:bold)

    save(savepath, fig)
    println("  → Saved: $savepath")
    return fig
end

function plot_bus_vuf_sorted(bus_ids, vufs; savepath="./plots/bus_vuf_sorted.pdf")
    order  = sortperm(vufs)
    sorted_vufs = vufs[order] .* 100
    xs = 1:length(sorted_vufs)

    fig = Figure(size=(900, 480), backgroundcolor=:white)
    ax = Axis(fig[1,1],
        xlabel="Bus index (sorted ascending by VUF)",
        ylabel="Voltage Unbalance Factor (%)",
        title="VUF at every bus -- natural network, no PV/STATCOM",
        titlesize=15, xlabelsize=12, ylabelsize=12,
        ygridvisible=true, xgridvisible=false, ygridcolor=(:black,0.08))

    lines!(ax, xs, sorted_vufs, color=:steelblue, linewidth=2.0)
    hlines!(ax, [2.0], color=:red, linestyle=:dash, linewidth=2.0,
            label="2% reference (typical grid-code style limit)")

    n_over = count(v -> v > 2.0, sorted_vufs)
    if n_over > 0
        text!(ax, length(xs)*0.05, maximum(sorted_vufs)*0.95,
              text="$n_over of $(length(xs)) buses exceed 2% VUF",
              fontsize=11, color=:red)
    end

    axislegend(ax, position=:lt, framevisible=true, labelsize=11)
    save(savepath, fig)
    println("  → Saved: $savepath")
    return fig
end

function plot_bus_vuf_histogram(vufs; savepath="./plots/bus_vuf_histogram.pdf")
    fig = Figure(size=(700, 450), backgroundcolor=:white)
    ax = Axis(fig[1,1],
        xlabel="Voltage Unbalance Factor (%)", ylabel="Number of buses",
        title="Distribution of VUF across the feeder",
        titlesize=15, xlabelsize=12, ylabelsize=12,
        ygridvisible=true, xgridvisible=false, ygridcolor=(:black,0.08))
    hist!(ax, vufs.*100, bins=25, color=(:steelblue,0.75), strokewidth=1, strokecolor=:white)
    vlines!(ax, [2.0], color=:red, linestyle=:dash, linewidth=2.0, label="2% reference")
    axislegend(ax, position=:rt, framevisible=true, labelsize=11)
    save(savepath, fig)
    println("  → Saved: $savepath")
    return fig
end

function plot_voltage_profile_by_phase(bus_ids, vm_ph; savepath="./plots/voltage_profile_by_phase.pdf")
    means = [mean(v) for v in vm_ph]
    order = sortperm(means)
    xs = 1:length(order)

    va = [vm_ph[i][1] for i in order]
    vb = [vm_ph[i][2] for i in order]
    vc = [vm_ph[i][3] for i in order]

    fig = Figure(size=(950, 480), backgroundcolor=:white)
    ax = Axis(fig[1,1],
        xlabel="Bus index (sorted by mean voltage)", ylabel="Voltage magnitude (pu)",
        title="Per-phase voltage profile -- natural network, no PV/STATCOM",
        titlesize=15, xlabelsize=12, ylabelsize=12,
        ygridvisible=true, xgridvisible=false, ygridcolor=(:black,0.08))

    lines!(ax, xs, va, color=:steelblue,   linewidth=1.8, label="Phase A")
    lines!(ax, xs, vb, color=:darkorange,  linewidth=1.8, label="Phase B")
    lines!(ax, xs, vc, color=:forestgreen, linewidth=1.8, label="Phase C")
    hlines!(ax, [1.0], color=(:black,0.25), linestyle=:dot, linewidth=1.2)

    axislegend(ax, position=:lb, framevisible=true, labelsize=11)
    save(savepath, fig)
    println("  → Saved: $savepath")
    return fig
end

# -----------------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------------
dm = load_base_network(data_path; load_multiplier=1.0, enforce_bounds=false)
SBASE_KVA = get_sbase_kva(dm)

println("="^80)
println(" NETWORK UNBALANCE VISUALISATION -- baseline, no PV, no STATCOM")
println(" sbase = $SBASE_KVA kVA")
println("="^80)

stats = load_phase_stats(dm, SBASE_KVA)
println("  Load count : A=$(stats.tally[1])  B=$(stats.tally[2])  C=$(stats.tally[3])  (3-phase: $(stats.threephase))")
println("  P demand   : A=$(round(stats.pd_kw[1],digits=1))kW  B=$(round(stats.pd_kw[2],digits=1))kW  C=$(round(stats.pd_kw[3],digits=1))kW")
println("  Q demand   : A=$(round(stats.qd_kvar[1],digits=1))kVAr  B=$(round(stats.qd_kvar[2],digits=1))kVAr  C=$(round(stats.qd_kvar[3],digits=1))kVAr")

println("\n  Solving natural network...")
status = build_and_solve(dm)
println("  Status: $status")

if status in [JuMP.LOCALLY_SOLVED, JuMP.OPTIMAL, JuMP.ALMOST_LOCALLY_SOLVED]
    bus_ids, vufs, vm_ph = all_bus_vuf_and_vmag()
    println("  Buses analysed: $(length(bus_ids))")
    println("  VUF   : mean=$(round(100*mean(vufs),digits=3))%  max=$(round(100*maximum(vufs),digits=3))%  " *
            "(worst bus: $(bus_ids[argmax(vufs)]))")
    n_over2 = count(v -> v > 0.02, vufs)
    println("  Buses exceeding 2% VUF: $n_over2 / $(length(bus_ids))")

    plot_load_phase_distribution(stats)
    plot_bus_vuf_sorted(bus_ids, vufs)
    plot_bus_vuf_histogram(vufs)
    plot_voltage_profile_by_phase(bus_ids, vm_ph)
else
    println("  WARNING: natural network did not solve -- cannot produce voltage-based plots.")
    plot_load_phase_distribution(stats)
end

println("\n  Plots saved to ./plots/")
println("  Done.")