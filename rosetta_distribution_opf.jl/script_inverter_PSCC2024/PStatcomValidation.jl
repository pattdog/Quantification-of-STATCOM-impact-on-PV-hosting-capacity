#=
==============================================================================
LOAD-MULTIPLIER SWEEP: STATCOM P-exchange hosting-capacity study (v5)
==============================================================================
v5 changes from v4:
  1. 10 STATCOMs instead of 1, spread evenly across the feeder using the
     same `spacing` convention as the seminar script's add_pv!/add_statcoms!
     -- step through the sorted load-bus list at a fixed interval so units
     land roughly evenly along the feeder rather than clustering at one bus.
  2. "warm_start_init_point" => "yes" REMOVED from the Ipopt config, per
     discussion: that flag tells Ipopt to skip its own initialization of
     bound multipliers and trust the supplied starting point, but no
     explicit starting values were being supplied for the STATCOM's new
     variables -- a mismatch that can make Ipopt's local search WORSE, not
     better, and was a leading suspect for the false LOCALLY_INFEASIBLE
     results seen in v4 (a mathematical proof exists that B/C's feasible
     region always contains A's -- see chat -- so any B/C infeasibility at
     a load level where A solves is necessarily a solver artifact, not a
     real network limit).
  3. max_iter raised to 50,000 (from 3,000) to give the solver more room
     now that warm-start isn't short-circuiting its own initialization.
  Explicit zero-starting of the STATCOMs' crg/cig (also discussed as a
  candidate fix) was NOT added, per instruction -- if false infeasibility
  persists after this run, that remains the next thing to try.

Dispatch reporting is now AGGREGATED across all 10 units (per-unit-per-phase
would be 30 lines per scenario per load level) -- net and gross P/Q per
phase, plus each unit's individual utilization (S_i/rating) to confirm the
normalized S^2 capability constraint is still being honored per-unit.

Three scenarios per load level:
    A) No STATCOM at all (baseline)
    B) 10x STATCOM, Q-only (statcom_p_exchange = false) -- old model
    C) 10x STATCOM, P-exchange enabled (statcom_p_exchange = true) -- new

Run with:
    julia --project=. statcom_load_sweep.jl
==============================================================================
=#

using Logging
Logging.disable_logging(Logging.Warn)
# ^ Suppresses ALL @warn-level output for the rest of this session. Blunt --
# would also hide any genuinely NEW warning -- but the known soft-scope
# ambiguity warnings from constraints_PBalance.jl/objectives.jl fire on
# every include() (up to 36 times across this sweep) and would otherwise
# bury the actual results. Comment out if debugging something specific.

using Pkg
Pkg.activate("./")
using rosetta_distribution_opf
import PowerModelsDistribution
import InfrastructureModels
using Ipopt
using JuMP
using Printf
using Statistics

const PMD  = PowerModelsDistribution
const RPMD = rosetta_distribution_opf
const IM   = InfrastructureModels
PMD.silence!()


ipopt_solver = JuMP.optimizer_with_attributes(
    Ipopt.Optimizer,
    "print_level" => 0,
    "sb"          => "yes",
    "max_iter"    => 50000,
    "warm_start_init_point" => "no", # claude says keep no
    "tol" => 1e-3,                   # Relax the main tolerance (1e-4 -> 1e-3)
    "acceptable_tol" => 1e-1,        # Be very forgiving if it gets close
    "constr_viol_tol" => 1e-3,       # Allow tiny overlaps in constraints
    
    # Advanced Scaling - Helps with 4-wire numerical issues
    "nlp_scaling_method" => "gradient-based", 
)

data_path = "./rosetta_distribution_opf.jl/data/ENWL_4w_Network1_Feeder1/Master.dss"

# -----------------------------------------------------------------------
# kW/kVAr <-> pu conversions
# -----------------------------------------------------------------------
kw_to_pu(p_kw, sbase_kva)     = p_kw / sbase_kva
kvar_to_pu(q_kvar, sbase_kva) = q_kvar / sbase_kva
pu_to_kw(p_pu, sbase_kva)     = p_pu * sbase_kva
pu_to_kvar(q_pu, sbase_kva)   = q_pu * sbase_kva

# -----------------------------------------------------------------------
# Network loader -- identical to the known-working driver.
# -----------------------------------------------------------------------
function load_base_network(data_path; load_multiplier=1.0, enforce_bounds=true)
    data_eng  = PMD.parse_file(data_path, transformations=[PMD.transform_loops!])
    # --- SBASE OVERRIDE ---
    # PowerModelsDistribution typically stores this in Volt-Amperes (e.g., 1e6 = 1 MVA)
    # Print data_eng["settings"]["sbase"] once to see its current value, 
    # then overwrite it to a smaller base like 100 kVA (100_000.0) or 1 MVA (1_000_000.0)
    # ----------------------
    data_math = PMD.transform_data_model(
        data_eng, multinetwork=false, kron_reduce=false, phase_project=false,
        make_pu=false   # skip PMD's own pu conversion
    )
    PMD.make_per_unit!(data_math; sbase=1000.0)   # do it yourself, with your chosen base

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

function total_load_kw_kvar(data_math, sbase_kva)
    total_pd_pu = sum(sum(l["pd"]) for (i,l) in data_math["load"])
    total_qd_pu = sum(sum(l["qd"]) for (i,l) in data_math["load"])
    return pu_to_kw(total_pd_pu, sbase_kva), pu_to_kvar(total_qd_pu, sbase_kva)
end

# -----------------------------------------------------------------------
# NEW: place N_UNITS STATCOMs evenly across the feeder, same `spacing`
# convention as the seminar script's add_pv!/add_statcoms! -- step through
# the sorted load-bus list at a fixed interval. Returns a Vector{Int} of
# gen ids.
# -----------------------------------------------------------------------
function add_statcoms!(data_math, sbase_kva; n_units, rating_kvar_each, p_exchange::Bool)
    load_ids = sort(
        collect(keys(data_math["load"])),
        by = x -> tryparse(Int, x) === nothing ? 0 : parse(Int, x)
    )
    n_loads = length(load_ids)
    spacing = max(1, div(n_loads, n_units))

    s_rated_pu = kvar_to_pu(rating_kvar_each, sbase_kva)

    gen_ids = Int[]
    placed  = 0
    idx = 1
    while placed < n_units && idx <= n_loads
        load = data_math["load"][load_ids[idx]]
        target_bus = load["load_bus"]

        gen_id = string(length(data_math["gen"]) + 1)
        data_math["gen"][gen_id] = deepcopy(data_math["gen"]["1"])
        gen = data_math["gen"][gen_id]

        gen["gen_bus"] = target_bus
        gen["type"]    = "STATCOM"
        gen["name"]    = "statcom_$(load_ids[idx])"
        gen["cost"]    = [0.001, 0.0]

        if p_exchange
            gen["pmax"] =  fill(s_rated_pu, 3)
            gen["pmin"] = -fill(s_rated_pu, 3)
            gen["qmax"] =  fill(s_rated_pu, 3)
            gen["qmin"] = -fill(s_rated_pu, 3)
            gen["statcom_p_exchange"] = true
            gen["s_rated"]            = fill(s_rated_pu, 3)
            gen["p_loss"]             = 0.0
        else
            gen["pmax"] = zeros(3)
            gen["pmin"] = zeros(3)
            gen["qmax"] =  fill(s_rated_pu, 3)
            gen["qmin"] = -fill(s_rated_pu, 3)
            gen["statcom_p_exchange"] = false
        end

        push!(gen_ids, parse(Int, gen_id))
        placed += 1
        idx += spacing
    end

    return gen_ids
end

# -----------------------------------------------------------------------
# Build + solve one scenario at one load level
# -----------------------------------------------------------------------
function build_and_solve(data_math)
    global ref = IM.build_ref(data_math, PMD.ref_add_core!, PMD._pmd_global_keys, PMD.pmd_it_name)[:it][:pmd][:nw][0]
    global model = JuMP.Model(ipopt_solver)

    include("./core/variables.jl")


    # === ADD: start values for STATCOM crg/cig, if any are present ===
    if @isdefined(STATCOM_GEN_IDS)
        for gid in STATCOM_GEN_IDS
            for p in 1:3
                JuMP.set_start_value(crg[p,gid], 0.0)
                JuMP.set_start_value(cig[p,gid], 0.0)
            end
        end
    end
    # === END ADD ===

    include("./core/constraints_PBalance.jl")

    global objective = "VUF"
    global objective_aggregation = :sum   # or :max
    include("./core/objectives_FIXED.jl")

    JuMP.optimize!(model)
    return JuMP.termination_status(model)
end

# -----------------------------------------------------------------------
# Worst-case VUF across every bus -- true ratio |V_neg|/|V_pos|.
# -----------------------------------------------------------------------
function worst_case_vuf()
    alpha = exp(im*2/3*pi)
    T = 1/3 * [1 1 1 ; 1 alpha alpha^2 ; 1 alpha^2 alpha]

    worst_vuf  = 0.0
    worst_bus  = nothing

    for (i, bus) in ref[:bus]
        terms = bus["terminals"]
        if !(1 in terms && 2 in terms && 3 in terms)
            continue
        end
        vr_val = [JuMP.value(vr[p,i]) for p in 1:3]
        vi_val = [JuMP.value(vi[p,i]) for p in 1:3]
        vph = vr_val .+ im .* vi_val
        v012 = T * vph
        vpos = abs(v012[2])
        vneg = abs(v012[3])
        if vpos > 1e-6
            vuf = vneg / vpos
            if vuf > worst_vuf
                worst_vuf = vuf
                worst_bus = i
            end
        end
    end

    return worst_vuf, worst_bus
end

# -----------------------------------------------------------------------
# Worst-case |V| magnitude across every bus (pu).
# -----------------------------------------------------------------------
function worst_case_vmag()
    min_v, min_bus, min_phase = Inf, nothing, nothing
    max_v, max_bus, max_phase = -Inf, nothing, nothing

    for (i, bus) in ref[:bus]
        get(bus, "bus_type", 1) == 3 && continue
        terms = bus["terminals"]
        for p in 1:3
            p in terms || continue
            vr_val = JuMP.value(vr[p,i])
            vi_val = JuMP.value(vi[p,i])
            vm = sqrt(vr_val^2 + vi_val^2)
            if vm < min_v
                min_v, min_bus, min_phase = vm, i, p
            end
            if vm > max_v
                max_v, max_bus, max_phase = vm, i, p
            end
        end
    end

    return (min_v=min_v, min_bus=min_bus, min_phase=min_phase,
            max_v=max_v, max_bus=max_bus, max_phase=max_phase)
end

# -----------------------------------------------------------------------
# NEW: aggregate dispatch across all N STATCOM units. Reports, per phase:
# net P/Q (signed sum -- can cancel across units) and gross P/Q (sum of
# |value| -- reflects how hard the fleet is actually working). Also
# reports the worst individual unit's utilization (S_i/rating) to confirm
# the normalized capability constraint is honored per-unit even in
# aggregate reporting.
# -----------------------------------------------------------------------
function aggregate_statcom_dispatch(gen_ids, sbase_kva, rating_kvar_each)
    isempty(gen_ids) && return nothing

    net_pg  = zeros(3); net_qg  = zeros(3)
    gross_pg = zeros(3); gross_qg = zeros(3)
    worst_util = 0.0
    worst_util_gen = nothing

    for gid in gen_ids
        pg_kw   = [pu_to_kw(JuMP.value(pg[p,gid]), sbase_kva) for p in 1:3]
        qg_kvar = [pu_to_kvar(JuMP.value(qg[p,gid]), sbase_kva) for p in 1:3]
        for p in 1:3
            net_pg[p]  += pg_kw[p];    net_qg[p]  += qg_kvar[p]
            gross_pg[p] += abs(pg_kw[p]); gross_qg[p] += abs(qg_kvar[p])
            s_i = sqrt(pg_kw[p]^2 + qg_kvar[p]^2)
            util = s_i / rating_kvar_each
            if util > worst_util
                worst_util = util
                worst_util_gen = gid
            end
        end
    end

    return (net_pg=net_pg, net_qg=net_qg, gross_pg=gross_pg, gross_qg=gross_qg,
            worst_util=worst_util, worst_util_gen=worst_util_gen, n_units=length(gen_ids))
end

function format_agg_dispatch(d)
    isnothing(d) && return "n/a"
    parts = String[]
    for p in 1:3
        push!(parts, @sprintf("ph%d: net P=%.3f Q=%.3f kVAr | gross P=%.3f Q=%.3f kVAr",
              p, d.net_pg[p], d.net_qg[p], d.gross_pg[p], d.gross_qg[p]))
    end
    header = @sprintf("[%d units, worst per-unit utilization = %.1f%% at gen %s]",
              d.n_units, 100*d.worst_util, string(d.worst_util_gen))
    return header * "  " * join(parts, "  |  ")
end

# -----------------------------------------------------------------------
# Sweep
# -----------------------------------------------------------------------
LOAD_MULTIPLIERS = [2.0, 3, 4.5, 6.0, 7.0]
N_STATCOMS          = 10
RATING_KVAR_EACH     = 20.0   # per-unit nameplate rating, real kVAr

scenarios = [
    ("A: No STATCOM",             :none),
    ("B: 10x STATCOM Q-only",     :qonly),
    ("C: 10x STATCOM P-exchange", :pexchange),
]

results = Dict{Symbol, Vector{NamedTuple}}(
    :none => [], :qonly => [], :pexchange => []
)

still_feasible = Dict(:none => true, :qonly => true, :pexchange => true)

_dm_probe   = load_base_network(data_path; load_multiplier=1.0)
SBASE_KVA   = get_sbase_kva(_dm_probe)
base_pd_kw, base_qd_kvar = total_load_kw_kvar(_dm_probe, SBASE_KVA)

println("="^90)
println(" LOAD-MULTIPLIER SWEEP (v5 -- 10 STATCOMs)")
println(" sbase = $SBASE_KVA kVA")
println(" Feeder total demand @ load_multiplier=1.0: $(round(base_pd_kw,digits=2)) kW / $(round(base_qd_kvar,digits=2)) kVAr")
println(" STATCOM rating: $N_STATCOMS units x $RATING_KVAR_EACH kVAr each = $(N_STATCOMS*RATING_KVAR_EACH) kVAr total fleet capacity")
println("="^90)

for lm in LOAD_MULTIPLIERS
    pd_kw, qd_kvar = base_pd_kw * lm, base_qd_kvar * lm
    println("\n── load_multiplier = $lm   (feeder demand: $(round(pd_kw,digits=1)) kW / $(round(qd_kvar,digits=1)) kVAr) ──")

    for (label, kind) in scenarios
        if !still_feasible[kind]
            println("  $label: skipped (already infeasible at a lower load level)")
            push!(results[kind], (lm=lm, status="SKIPPED", vuf=NaN, vuf_bus=nothing,
                                   vmag=nothing, dispatch=nothing))
            continue
        end

        dm = load_base_network(data_path; load_multiplier=lm, enforce_bounds=true)

        global STATCOM_GEN_IDS = Int[]  # reset before every scenario, then overwrite below if kind != :none
        statcom_gen_ids = Int[]
        if kind != :none
            statcom_gen_ids = add_statcoms!(dm, SBASE_KVA; n_units=N_STATCOMS,
                                            rating_kvar_each=RATING_KVAR_EACH,
                                            p_exchange=(kind == :pexchange))
        end
        global STATCOM_GEN_IDS = statcom_gen_ids

        status = try
            build_and_solve(dm)
        catch e
            println("  $label: ERROR during build/solve -- $e")
            :ERROR
        end

        if status in [JuMP.LOCALLY_SOLVED, JuMP.OPTIMAL, JuMP.ALMOST_LOCALLY_SOLVED]
            vuf, vuf_bus = worst_case_vuf()
            vmag = worst_case_vmag()
            dispatch = aggregate_statcom_dispatch(statcom_gen_ids, SBASE_KVA, RATING_KVAR_EACH)

            println("  $label: SOLVED")
            println("    VUF   : worst = $(round(100*vuf, digits=3))%  (bus $vuf_bus)")
            println("    |V|   : min = $(round(vmag.min_v, digits=4)) pu (bus $(vmag.min_bus), ph$(vmag.min_phase))" *
                    "   max = $(round(vmag.max_v, digits=4)) pu (bus $(vmag.max_bus), ph$(vmag.max_phase))" *
                    "   [bounds: 0.90 - 1.10]")
            if !isnothing(dispatch)
                println("    STATCOM dispatch: $(format_agg_dispatch(dispatch))")
            end

            push!(results[kind], (lm=lm, status="SOLVED", vuf=vuf, vuf_bus=vuf_bus,
                                   vmag=vmag, dispatch=dispatch))
        else
            println("  $label: $status  -- treating as infeasible, will skip further load levels")
            push!(results[kind], (lm=lm, status=string(status), vuf=NaN, vuf_bus=nothing,
                                   vmag=nothing, dispatch=nothing))
            still_feasible[kind] = false
        end
    end
end

# -----------------------------------------------------------------------
# Summary tables
# -----------------------------------------------------------------------
println("\n" * "="^90)
println(" SUMMARY -- VUF and voltage magnitude")
println("="^90)
@printf("  %-6s  %-3s  %-10s  %-22s  %-22s\n", "LM", "Scn", "Status", "VUF (worst)", "|V| min / max (pu)")
println("  " * "-"^80)
for kind in [:none, :qonly, :pexchange]
    label = Dict(:none=>"A", :qonly=>"B", :pexchange=>"C")[kind]
    for r in results[kind]
        if r.status == "SOLVED"
            @printf("  %-6s  %-3s  %-10s  %-22s  %-22s\n",
                r.lm, label, r.status,
                "$(round(100*r.vuf,digits=3))%",
                "$(round(r.vmag.min_v,digits=4)) / $(round(r.vmag.max_v,digits=4))")
        else
            @printf("  %-6s  %-3s  %-10s  %-22s  %-22s\n", r.lm, label, r.status, "-", "-")
        end
    end
end

println("\n" * "="^90)
println(" SUMMARY -- Aggregate STATCOM dispatch (scenarios B and C only)")
println("="^90)
for kind in [:qonly, :pexchange]
    label = Dict(:qonly=>"B (Q-only)", :pexchange=>"C (P-exchange)")[kind]
    println("\n  -- $label --  ($N_STATCOMS units x $RATING_KVAR_EACH kVAr each)")
    for r in results[kind]
        if r.status == "SOLVED" && !isnothing(r.dispatch)
            println("  lm=$(r.lm):  $(format_agg_dispatch(r.dispatch))")
        end
    end
end

println("\nFirst infeasible load_multiplier per scenario:")
for (label, kind) in scenarios
    first_bad = findfirst(r -> r.status != "SOLVED" && r.status != "SKIPPED", results[kind])
    if isnothing(first_bad)
        println("  $label: never went infeasible across the tested range")
    else
        println("  $label: first infeasible at load_multiplier = $(results[kind][first_bad].lm)")
    end
end