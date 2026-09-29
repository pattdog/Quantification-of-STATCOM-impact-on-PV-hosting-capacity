#=
==============================================================================
PV STRESS SWEEP: single-phase, phase-matched PV, with/without STATCOM (v1)
==============================================================================
Built directly on the validated load-multiplier sweep pipeline (raw includes
of variables.jl / constraints_PBalance.jl / objectives_FIXED.jl -- NOT
PMD.instantiate_mc_model/build_mc_opf, since that path has no knowledge of
the statcom_p_exchange patch). Every fix validated in that sweep is carried
over unchanged:
  - sbase set via make_pu=false + make_per_unit!(sbase=...), NOT
    data_eng["settings"]["sbase"] (confirmed not to propagate)
  - add_start_vrvi! for flat positive-sequence voltage start
  - explicit crg/cig = 0 starts for every STATCOM unit
  - warm_start_init_point = "no"
  - objectives_FIXED.jl's vm_012-degeneracy fix (only relevant if you later
    switch objective away from "cost")

NEW IN THIS SCRIPT: single-phase, phase-matched PV.
  add_pv! places one PV unit per (non-source) load bus, connected ONLY to
  the phase(s) the colocated load itself is connected to -- e.g. a load on
  phase B gets a phase-B-only PV gen, not a balanced 3-phase injection. This
  matches real single-phase residential rooftop PV, which is physically
  wired to the same phase as the household's supply.

  PV is P-only, unity power factor (no smart-inverter Q capability):
  qmin = qmax = 0. This is a deliberate simplification for this study --
  the point is to see how much unbalance PV ALONE introduces, before any
  Q-support device is added.

Network stays at load_multiplier = 1.0 throughout -- this sweep stresses the
network via increasing PV penetration (pv_kw per unit), not load growth.
Three scenarios per PV level, exactly mirroring the STATCOM sweep's A/B/C
convention:
    A) PV only, no STATCOM
    B) PV + STATCOM, Q-only (statcom_p_exchange = false)
    C) PV + STATCOM, P-exchange (statcom_p_exchange = true)

Run with:
    julia --project=. pv_statcom_sweep.jl
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

const PMD  = PowerModelsDistribution
const RPMD = rosetta_distribution_opf
const IM   = InfrastructureModels
PMD.silence!()

ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0,"sb"=> "yes",)

#=
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
    #"nlp_scaling_method" => "gradient-based", 
    "nlp_scaling_method" => "none", 
)
=#
data_path = "./rosetta_distribution_opf.jl/data/ENWL_4w_Network1_Feeder1/Master.dss"

# -----------------------------------------------------------------------
# kW/kVAr <-> pu conversions
# -----------------------------------------------------------------------
kw_to_pu(p_kw, sbase_kva)     = p_kw / sbase_kva
kvar_to_pu(q_kvar, sbase_kva) = q_kvar / sbase_kva
pu_to_kw(p_pu, sbase_kva)     = p_pu * sbase_kva
pu_to_kvar(q_pu, sbase_kva)   = q_pu * sbase_kva


# -----------------------------------------------------------------------
# MEN (Multiple Earthed Neutral) -- AS/NZS 3000 earthing arrangement.
#
# Built directly as ENG shunt dicts rather than via DSS `reactor` objects,
# because PMD's single-terminal reactor path ignores R/X and falls back to
# the kvar/kv defaults (produces bs=0.00772 S capacitive, gs=0). The
# two-terminal form (bus1=X.4 bus2=X.0) parses R correctly but materialises
# an explicit earth node, widening conductor_ids to [1,2,3,4,5].
#
# Schema below copied from what PMD itself emits for an ENG shunt.
# ENG is SI: gs in siemens. make_per_unit! converts (zbase = vbase^2/sbase).
# -----------------------------------------------------------------------
function add_men!(data_eng; r_load_ohm=10.0, verbose=true)
    haskey(data_eng, "shunt") || (data_eng["shunt"] = Dict{String,Any}())

    men_buses = String[]
    for (_, load) in data_eng["load"]
        b = load["bus"]
        4 in load["connections"] || continue      # needs a neutral to bond
        b in men_buses || push!(men_buses, b)
    end
    sort!(men_buses)

    g = 1.0 / r_load_ohm
    for b in men_buses
        data_eng["shunt"]["men_$b"] = Dict{String,Any}(
            "source_id"     => "reactor.men_$b",
            "status"        => PMD.ENABLED,
            "model"         => PMD.REACTOR,
            "connections"   => [4],
            "dispatchable"  => PMD.NO,
            "gs"            => fill(g, 1, 1),
            "bs"            => zeros(1, 1),
            "bus"           => b,
            "configuration" => PMD.WYE,
        )
    end

    verbose && println("  MEN: $(length(men_buses)) bonds @ $(r_load_ohm) Ω (g = $(round(g,digits=4)) S)")
    return men_buses
end

# -----------------------------------------------------------------------
# Network loader -- sbase fix via make_pu=false + make_per_unit!(sbase=...),
# validated in the STATCOM load-multiplier sweep.
# -----------------------------------------------------------------------
function load_base_network(data_path; load_multiplier=1.0, enforce_bounds=true, sbase_kva=100, men_ohm=nothing)
    data_eng  = PMD.parse_file(data_path, transformations=[PMD.transform_loops!])
    # --- MEN bonds, added in ENG so make_per_unit! handles the conversion
    isnothing(men_ohm) || add_men!(data_eng; r_load_ohm=men_ohm, verbose=false)
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

function total_load_kw_kvar(data_math, sbase_kva)
    total_pd_pu = sum(sum(l["pd"]) for (i,l) in data_math["load"])
    total_qd_pu = sum(sum(l["qd"]) for (i,l) in data_math["load"])
    return pu_to_kw(total_pd_pu, sbase_kva), pu_to_kvar(total_qd_pu, sbase_kva)
end

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
    println("  Load phase distribution: Phase A=$(tally[1])  Phase B=$(tally[2])  Phase C=$(tally[3])  Three-phase=$threephase")
    return tally
end

# -----------------------------------------------------------------------
# PV PLACEMENT -- single-phase, phase-matched to the colocated load.
# P-only, unity power factor (qmin=qmax=0) -- no smart-inverter Q capability,
# deliberately, so this sweep isolates PV's own unbalance impact before any
# Q-support device (STATCOM) is added on top.
#
# BUGFIX vs. earlier draft: the "trim any other length-3 template field"
# loop now indexes by the ACTUAL phase numbers (phase_conns), not by
# position 1:n_ph. Indexing by position silently mislabels inherited fields
# for any PV unit not on phase 1 (e.g. a phase-B PV would have picked up
# phase-A's values for any leftover length-3 field from the deepcopy'd
# template). Indexing by phase_conns is correct regardless of which phase
# the unit lands on.
# -----------------------------------------------------------------------
function add_pv!(data_math, sbase_kva; pv_kw, spacing=1, pv_cost=-1000.0)
    pv_pu = kw_to_pu(pv_kw, sbase_kva)

    source_buses = Set([i for (i, bus) in data_math["bus"] if bus["bus_type"] == 3])
    load_ids = sort(
        collect(keys(data_math["load"])),
        by = x -> tryparse(Int, x) === nothing ? 0 : parse(Int, x)
    )

    pv_ids = Int[]
    n_ph1  = 0
    n_skip = 0

    for i in 1:spacing:length(load_ids)
        load       = data_math["load"][load_ids[i]]
        target_bus = string(load["load_bus"])
        target_bus ∈ source_buses && continue

        conns       = load["connections"]
        phase_conns = filter(c -> c != 4, conns)
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
        gen["connections"] = deepcopy(conns)

        gen["pmax"] = pv_pu * ones(n_ph)
        gen["pmin"] = zeros(n_ph)
        gen["qmax"] = zeros(n_ph)   # unity PF -- no Q capability
        gen["qmin"] = zeros(n_ph)
        gen["cost"] = [pv_cost, 0.0]

        # BUGFIX: index by actual phase numbers, not position.
        for (k, v) in gen
            if v isa Vector{<:Real} && length(v) == 3 && n_ph != 3
                gen[k] = v[phase_conns]
            end
        end

        push!(pv_ids, parse(Int, gen_id))
        n_ph1 += 1
    end

    println("  PV: $(length(pv_ids)) units  pv_kw=$(pv_kw) kW/unit  spacing=$(spacing)  cost=$(pv_cost)  " *
            "(single-phase: $n_ph1, skipped no-phase-match: $n_skip)")
    return pv_ids
end

# -----------------------------------------------------------------------
# STATCOM PLACEMENT -- merges the two prior versions: kept 3-phase, sized
# directly in real kVAr as a TOTAL nameplate rating split evenly across the
# 3 legs (matching the phase-realistic script's convention), but retains the
# n_units/p_exchange structure from the validated load-multiplier sweep,
# since that's what constraints_PBalance.jl's statcom_p_exchange patch
# actually reads (s_rated, p_loss, statcom_p_exchange fields).
# -----------------------------------------------------------------------
function add_statcoms!(data_math, sbase_kva; n_units, rating_kvar_total, p_exchange::Bool)
    load_ids = sort(
        collect(keys(data_math["load"])),
        by = x -> tryparse(Int, x) === nothing ? 0 : parse(Int, x)
    )
    n_loads = length(load_ids)
    spacing = max(1, div(n_loads, n_units))

    q_total_pu = kvar_to_pu(rating_kvar_total, sbase_kva)
    q_leg_pu   = q_total_pu / 3   # even split across the 3 converter legs

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
            gen["pmax"] =  fill(q_leg_pu, 3)
            gen["pmin"] = -fill(q_leg_pu, 3)
            gen["qmax"] =  fill(q_leg_pu, 3)
            gen["qmin"] = -fill(q_leg_pu, 3)
            gen["statcom_p_exchange"] = true
            gen["s_rated"]            = fill(q_leg_pu, 3)
            gen["p_loss"]             = 0.0
        else
            gen["pmax"] = zeros(3)
            gen["pmin"] = zeros(3)
            gen["qmax"] =  fill(q_leg_pu, 3)
            gen["qmin"] = -fill(q_leg_pu, 3)
            gen["statcom_p_exchange"] = false
        end

        push!(gen_ids, parse(Int, gen_id))
        placed += 1
        idx += spacing
    end

    println("  STATCOM: $(length(gen_ids)) units  rating=$(rating_kvar_total) kVAr/unit total " *
            "($(round(q_leg_pu*sbase_kva,digits=2)) kVAr/leg)  p_exchange=$(p_exchange)")
    return gen_ids
end

# -----------------------------------------------------------------------
# Build + solve -- raw includes, no instantiate_mc_model/build_mc_opf.
# crg/cig starts for STATCOM_GEN_IDS carried over from the validated
# load-multiplier sweep pipeline.
# -----------------------------------------------------------------------
function build_and_solve(data_math)
    #=
    global ref = IM.build_ref(data_math, PMD.ref_add_core!, PMD._pmd_global_keys, PMD.pmd_it_name)[:it][:pmd][:nw][0]
    global model = JuMP.Model(ipopt_solver)
    =#
    pm = PMD.instantiate_mc_model(
    data_math,
    PMD.IVRENPowerModel,        # or whichever type your fork's formulation uses
    (pm) -> nothing;            # no-op build method
    jump_model = JuMP.Model(ipopt_solver),
    ref_extensions = Function[]
)

    global ref   = pm.ref[:it][:pmd][:nw][0]
    global model = pm.model
    include("../core/variables.jl")

    #include("./core/constraints_PBalance.jl")
    include("../core/constraints_unified.jl")
    global objective = "cost"
    println("    [objective = \"$objective\"]")
    global objective_target_gens = [i for (i, g) in ref[:gen] if get(g, "type", "") == "STATCOM"]
    include("../core/objectives_FIXED.jl")

    JuMP.optimize!(model)
    return JuMP.termination_status(model)
end


# -----------------------------------------------------------------------
# Reporting helpers
# -----------------------------------------------------------------------
function worst_case_vuf(; bus_ids=nothing)
    alpha = exp(im*2/3*pi)
    T = 1/3 * [1 1 1 ; 1 alpha alpha^2 ; 1 alpha^2 alpha]
    worst_vuf, worst_bus = 0.0, nothing
    for (i, bus) in ref[:bus]
        !isnothing(bus_ids) && !(i in bus_ids) && continue
        terms = bus["terminals"]
        !(1 in terms && 2 in terms && 3 in terms) && continue
        vr_val = [JuMP.value(vr[p,i]) for p in 1:3]
        vi_val = [JuMP.value(vi[p,i]) for p in 1:3]
        vph = vr_val .+ im .* vi_val
        v012 = T * vph
        vpos, vneg = abs(v012[2]), abs(v012[3])
        if vpos > 1e-6 && vneg/vpos > worst_vuf
            worst_vuf, worst_bus = vneg/vpos, i
        end
    end
    return worst_vuf, worst_bus
end

function worst_case_vmag(; bus_ids=nothing)
    min_v, min_bus, min_phase = Inf, nothing, nothing
    max_v, max_bus, max_phase = -Inf, nothing, nothing
    for (i, bus) in ref[:bus]
        !isnothing(bus_ids) && !(i in bus_ids) && continue
        get(bus, "bus_type", 1) == 3 && continue
        terms = bus["terminals"]
        for p in 1:3
            p in terms || continue
            vm = sqrt(JuMP.value(vr[p,i])^2 + JuMP.value(vi[p,i])^2)
            if vm < min_v; min_v, min_bus, min_phase = vm, i, p; end
            if vm > max_v; max_v, max_bus, max_phase = vm, i, p; end
        end
    end
    return (min_v=min_v, min_bus=min_bus, min_phase=min_phase, max_v=max_v, max_bus=max_bus, max_phase=max_phase)
end

function pv_dispatch(pv_gen_ids, sbase_kva)
    isempty(pv_gen_ids) && return nothing
    output_kw = 0.0
    capacity_kw = 0.0
    for gid in pv_gen_ids
        gen = ref[:gen][gid]
        n_ph = length(gen["connections"]) - 1   # exclude neutral
        phases = gen["connections"][1:end-1]
        output_kw   += sum(pu_to_kw(JuMP.value(pg[p,gid]), sbase_kva) for p in phases)
        capacity_kw += pu_to_kw(gen["pmax"][1] * n_ph, sbase_kva)
    end
    util = capacity_kw > 0 ? output_kw / capacity_kw * 100 : NaN
    return (output_kw=output_kw, capacity_kw=capacity_kw, util=util, n_units=length(pv_gen_ids))
end

function aggregate_statcom_dispatch(gen_ids, sbase_kva, rating_kvar_total)
    isempty(gen_ids) && return nothing
    net_pg  = zeros(3); net_qg  = zeros(3)
    gross_pg = zeros(3); gross_qg = zeros(3)
    worst_util = 0.0
    worst_util_gen = nothing
    for gid in gen_ids
        pg_kw   = [pu_to_kw(JuMP.value(pg[p,gid]), sbase_kva) for p in 1:3]
        qg_kvar = [pu_to_kvar(JuMP.value(qg[p,gid]), sbase_kva) for p in 1:3]
        leg_rating_kvar = rating_kvar_total / 3
        for p in 1:3
            net_pg[p]  += pg_kw[p];    net_qg[p]  += qg_kvar[p]
            gross_pg[p] += abs(pg_kw[p]); gross_qg[p] += abs(qg_kvar[p])
            s_i = sqrt(pg_kw[p]^2 + qg_kvar[p]^2)
            util = s_i / leg_rating_kvar
            if util > worst_util
                worst_util, worst_util_gen = util, gid
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
        push!(parts, @sprintf("ph%d: net P=%.3f kW Q=%.3f kVAr | gross P=%.3f kW Q=%.3f kVAr",
              p, d.net_pg[p], d.net_qg[p], d.gross_pg[p], d.gross_qg[p]))
    end
    header = @sprintf("[%d units, worst per-unit utilization = %.1f%% at gen %s]",
              d.n_units, 100*d.worst_util, string(d.worst_util_gen))
    return header * "  " * join(parts, "  |  ")
end

# -----------------------------------------------------------------------
# Sweep -- load_multiplier fixed at 1.0; PV size is the stress axis.
# -----------------------------------------------------------------------
PV_KW_LEVELS  = [1.0, 3.0, 5.0, 7.0]#, 10.0, 25, 40]
N_STATCOMS    = 1
STATCOM_KVAR  = 60   # total nameplate per unit, kVAr
scenarios = [
    ("A: PV only",                 :none),
    ("B: PV + STATCOM Q-only",     :qonly),
    ("C: PV + STATCOM P-exchange", :pexchange),
]

results = Dict{Symbol, Vector{NamedTuple}}(:none => [], :qonly => [], :pexchange => [])
still_feasible = Dict(:none => true, :qonly => true, :pexchange => true)

_dm_probe = load_base_network(data_path; load_multiplier=1.0)
SBASE_KVA = get_sbase_kva(_dm_probe)
base_pd_kw, base_qd_kvar = total_load_kw_kvar(_dm_probe, SBASE_KVA)

println("="^90)
println(" PV STRESS SWEEP (v1 -- single-phase, phase-matched PV; load_multiplier = 1.0)")
println(" sbase = $SBASE_KVA kVA")
println(" Feeder total demand: $(round(base_pd_kw,digits=2)) kW / $(round(base_qd_kvar,digits=2)) kVAr")
println(" STATCOM: $N_STATCOMS units x $STATCOM_KVAR kVAr total each = $(N_STATCOMS*STATCOM_KVAR) kVAr total fleet capacity")
report_load_phase_distribution(_dm_probe)
println("="^90)

for pv_kw in PV_KW_LEVELS
    println("\n── pv_kw = $pv_kw kW/unit ──")

    for (label, kind) in scenarios
        if !still_feasible[kind]
            println("  $label: skipped (already infeasible at a lower PV level)")
            push!(results[kind], (pv_kw=pv_kw, status="SKIPPED", vuf=NaN, vuf_bus=nothing,
                                   vmag=nothing, pv=nothing, dispatch=nothing))
            continue
        end

        dm = load_base_network(data_path; load_multiplier=1.0, enforce_bounds=true)
        pv_gen_ids = add_pv!(dm, SBASE_KVA; pv_kw=pv_kw)

        statcom_gen_ids = Int[]
        if kind != :none
            statcom_gen_ids = add_statcoms!(dm, SBASE_KVA; n_units=N_STATCOMS,
                                             rating_kvar_total=STATCOM_KVAR,
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

            # STATCOM PCC (point of common coupling) view -- what the
            # STATCOM's own (local-only) objective actually sees, as
            # opposed to the network-wide worst case above.
            statcom_bus_ids = isempty(STATCOM_GEN_IDS) ? Int[] : [ref[:gen][gid]["gen_bus"] for gid in STATCOM_GEN_IDS]
            vuf_local = NaN
            vmag_local = nothing
            if !isempty(statcom_bus_ids)
                vuf_local, _ = worst_case_vuf(bus_ids=statcom_bus_ids)
                vmag_local = worst_case_vmag(bus_ids=statcom_bus_ids)
            end

            pv   = pv_dispatch(pv_gen_ids, SBASE_KVA)
            dispatch = aggregate_statcom_dispatch(statcom_gen_ids, SBASE_KVA, STATCOM_KVAR)

            println("  $label: SOLVED")
            println("    VUF   : worst = $(round(100*vuf, digits=3))%  (bus $vuf_bus)")
            if !isempty(statcom_bus_ids)
                println("    VUF (STATCOM PCC) : $(round(100*vuf_local, digits=3))%  (bus $(statcom_bus_ids[1]))")
            end
            println("    |V|   : min = $(round(vmag.min_v, digits=4)) pu (bus $(vmag.min_bus), ph$(vmag.min_phase))" *
                    "   max = $(round(vmag.max_v, digits=4)) pu (bus $(vmag.max_bus), ph$(vmag.max_phase))" *
                    "   [bounds: 0.90 - 1.10]")
            if !isnothing(vmag_local)
                println("    |V| (STATCOM PCC) : min = $(round(vmag_local.min_v, digits=4)) pu (ph$(vmag_local.min_phase))" *
                        "   max = $(round(vmag_local.max_v, digits=4)) pu (ph$(vmag_local.max_phase))")
            end
            if !isnothing(pv)
                println("    PV output / capacity : $(round(pv.output_kw,digits=1)) / $(round(pv.capacity_kw,digits=1)) kW" *
                        "   ($(round(pv.util,digits=1))% utilisation, $(pv.n_units) units)")
            end
            if !isnothing(dispatch)
                println("    STATCOM dispatch: $(format_agg_dispatch(dispatch))")
            end

            push!(results[kind], (pv_kw=pv_kw, status="SOLVED", vuf=vuf, vuf_bus=vuf_bus,
            vuf_local=vuf_local,   # ADD
            vmag=vmag, pv=pv, dispatch=dispatch))
        else
            println("  $label: $status  -- treating as infeasible, will skip further PV levels")
            push!(results[kind], (pv_kw=pv_kw, status=string(status), vuf=NaN, vuf_bus=nothing,
            vuf_local=NaN,   # ADD
            vmag=nothing, pv=nothing, dispatch=nothing))
            still_feasible[kind] = false
        end
    end
end

# -----------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------
println("\n" * "="^130)
println(" SUMMARY -- VUF, voltage magnitude, PV utilisation (grouped by PV level)")
println("="^130)
@printf("  %-7s  %-3s  %-10s  %-12s  %-12s  %-20s  %-18s  %-30s\n",
        "PV_kW", "Scn", "Status", "VUF (net)", "VUF (PCC)", "|V| min/max (pu)", "PV util", "STATCOM util")
println("  " * "-"^130)

for pv_kw in PV_KW_LEVELS
    for (label, kind) in scenarios
        scn_label = Dict(:none=>"A", :qonly=>"B", :pexchange=>"C")[kind]
        r = findfirst(x -> x.pv_kw == pv_kw, results[kind])
        isnothing(r) && continue
        row = results[kind][r]

        if row.status == "SOLVED"
            vuf_pcc_str = kind == :none ? "-" : "$(round(100*row.vuf_local, digits=3))%"

            pv_util_str = isnothing(row.pv) ? "-" :
                "$(round(row.pv.util,digits=1))% ($(round(row.pv.output_kw,digits=1))kW)"

            statcom_util_str = "-"
            if !isnothing(row.dispatch)
                statcom_kw   = sum(row.dispatch.gross_pg)
                statcom_kvar = sum(row.dispatch.gross_qg)
                statcom_util_str = "$(round(100*row.dispatch.worst_util,digits=1))% " *
                                    "($(round(statcom_kvar,digits=1))kVAr, $(round(statcom_kw,digits=1))kW)"
            end

            @printf("  %-7s  %-3s  %-10s  %-12s  %-12s  %-20s  %-18s  %-30s\n",
                pv_kw, scn_label, row.status,
                "$(round(100*row.vuf,digits=3))%",
                vuf_pcc_str,
                "$(round(row.vmag.min_v,digits=4)) / $(round(row.vmag.max_v,digits=4))",
                pv_util_str,
                statcom_util_str)
        else
            @printf("  %-7s  %-3s  %-10s  %-12s  %-12s  %-20s  %-18s  %-30s\n",
                pv_kw, scn_label, row.status, "-", "-", "-", "-", "-")
        end
    end
    println("  " * "-"^130)   # visual separator between PV levels
end