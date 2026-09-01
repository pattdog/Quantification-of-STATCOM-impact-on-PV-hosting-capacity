#=
# ==============================================================================
# SCRIPT: STATCOM & PV Hosting Capacity Assessment — VVC-Enabled PV
# PROJECT: Undergraduate Thesis - Power Systems Engineering
# AUTHOR: Pat
# DATE: Aug 2026
# ==============================================================================
#
# THIS VERSION ON TOP OF single_phase_pv.jl:
#
#  1. GENERATOR BOUNDS FIXED (PV only): the independent pmax/qmax box has
#     been replaced with the standard joint apparent-power capability curve,
#     P^2 + Q^2 <= S_rating^2, matching Eq. (1p) in Quiertant et al. 2023 and
#     the same constraint idiom already used for branch thermal limits in
#     constraints.jl. The old box bounds were tight in one direction (the
#     P=full-output corner sat exactly on the true S-circle) but understated
#     the inverter's real capability everywhere else, and — more importantly
#     — P and Q were governed by two independent constraints, so absorbing Q
#     never actually cost the inverter any P headroom. The new joint
#     constraint makes that trade-off real. add_pv! now stores a single
#     s_rating_pu per unit (apparent power) instead of deriving a tight
#     q_lim; pmax/qmax become a loose outer box (never binding) and the real
#     shaping happens via a NEW CONSTRAINT that must be added to
#     constraints.jl (see the accompanying patch file / notes below — this
#     part cannot be added from the top-level script, since it lives inside
#     the generator loop in your rosetta_distribution_opf.jl package).
#
#  2. VOLT-VAR DROOP (PV only, STATCOM untouched): each single-phase PV
#     unit's Q is no longer a free OPF decision variable — it's pinned by a
#     per-phase, local, autonomous Volt-Var droop law following the generic
#     VVC curve in Fig. 1 of Quiertant et al. 2023 / Eq. (4) of the Mhanna
#     et al. VVWO paper, using AS/NZS 4777.2-style breakpoints:
#         V1=207V  V2=220V  V3=240V  V4=253V
#     reacting only to that PV unit's OWN phase-to-neutral voltage
#     magnitude — no visibility into any other bus or unit. This again
#     requires a new constraint block in constraints.jl (registered as a
#     JuMP user-defined function so Ipopt can autodiff through the
#     piecewise droop curve — the "nonsmooth" approach from both papers).
#
#  3. STATCOM IS DELIBERATELY UNCHANGED in this script, per discussion with
#     supervisor: STATCOM Q stays a free, centrally-optimized OPF variable
#     for now (Option A). The STATCOM-side realism question (sequence-
#     voltage-minimizing objective, per supervisor guidance) is separate
#     follow-up work, not part of this script.
#
#  4. gen["vvc"] = true/false toggles the droop per PV unit, so you can run
#     UPF (Q=0 baseline), free-Q (old OPF-coordinated ceiling), and droop
#     (realistic autonomous VVC) side by side — see Case 2 below.
#
# IMPORTANT: item 1 and item 2 both require a corresponding block inside
# constraints.jl's generator loop. This script alone is NOT sufficient —
# see the accompanying constraints_patch_vvc.jl file for the exact code to
# insert into your rosetta_distribution_opf.jl package, and where.
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
#  [CORRECTED] pmax stays at the TRUE real-power nameplate (pv_pu) as a hard
#  cap — only qmax/qmin are loosened to +/-s_rating. An earlier version of
#  this function loosened pmax itself, which let P reach s_rating (25% above
#  nameplate at Q=0) in every case — caught via a flat 125% utilisation
#  figure across every run regardless of voltage bounds.
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

using Logging
Logging.disable_logging(Logging.Warn)
# ^ Suppresses ALL @warn-level output for the rest of this session — soft-
# scope ambiguity warnings, JuMP.Containers "single-element axis" notices,
# anything else at Warn level or below. This is blunt (it would also hide
# any genuinely new warning), but the alternative is chasing each one down
# individually as they surface, which isn't a good use of time while we're
# trying to just see solved results. Comment this out later if you want
# warnings back for a specific debugging session.

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
# PV PLACEMENT — REWRITTEN AGAIN: S_rating-based, VVC-ready
#
#   BOUNDS FIX: gen["pmax"]/gen["qmax"] are now a LOOSE outer box
#   (±s_rating on each axis) that should never itself bind. The
#   real shaping comes from the joint capability constraint
#       pg[p,id]^2 + qg[p,id]^2 <= s_rating^2
#   which must be added in constraints.jl (see constraints_patch_vvc.jl).
#   This script only stores the s_rating value each unit needs;
#   it cannot add a JuMP constraint from outside the model-builder.
#
#   VVC: each unit optionally carries a droop toggle (gen["vvc"])
#   and breakpoints (gen["vvc_v1..v4"], in PER-UNIT voltage relative
#   to this network's own 0.2309 kV phase-to-neutral base) plus a max
#   |Q| fraction of s_rating (gen["vvc_qbar"]). When vvc=false, Q is
#   left as a free OPF variable (bounded only by the loose box above)
#   — i.e. the old Option-A "smart inverter ceiling" behavior.
#   When vvc=true, constraints.jl's new block pins Q to the droop law
#   instead — realistic, local, autonomous AS/NZS 4777.2-style control.
# ═══════════════════════════════════════════════════════════════
function add_pv!(data_math;
        pv_kw      = 5.0,       # real active rating per installed unit, kW
        s_scale    = 1.25,      # apparent-power headroom multiplier: s_rating = s_scale*pv_kw
        spacing    = 1,
        pv_cost    = -1000.0,
        enable_vvc = true,      # true: local Volt-Var droop (needs constraints.jl patch)
                                 # false: Q free in OPF (old Option-A ceiling)
        vvc_v1_v   = 207.0,     # VVC breakpoints, REAL VOLTS, phase-to-neutral
        vvc_v2_v   = 220.0,
        vvc_v3_v   = 240.0,
        vvc_v4_v   = 253.0,
        vvc_qbar   = 0.44,      # max |Q| as a fraction of s_rating at V<=V1 / V>=V4
        vnom_kv    = 0.2309,    # network's nominal phase-to-neutral base, from Loads_v2.txt
        sbase_kva  = SBASE_KVA)

    pv_pu = kw_to_pu(pv_kw, sbase_kva)
    s_pu  = s_scale * pv_pu   # apparent power rating, per-unit — this IS s_rating

    # VVC breakpoints: real volts -> per-unit, relative to this network's own base
    vnom_v = vnom_kv * 1000
    v1_pu  = vvc_v1_v / vnom_v
    v2_pu  = vvc_v2_v / vnom_v
    v3_pu  = vvc_v3_v / vnom_v
    v4_pu  = vvc_v4_v / vnom_v

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
        gen["connections"] = deepcopy(conns)       # matches the load's phase(s)

        # BOUNDS FIX: pmax stays at the TRUE real-power nameplate (pv_pu) —
        # a 5kW inverter cannot legitimately output more than 5kW of real
        # power just because it has apparent-power headroom for Q. Only
        # qmax/qmin get loosened to +/-s_pu; the joint S^2 constraint added
        # in constraints.jl then curtails P below pv_pu ONLY once |Q| grows
        # large enough that the S-circle, not this P cap, becomes binding.
        # (Loosening pmax itself was a bug in the previous version of this
        # function — it let P reach s_pu = s_scale*pv_pu at Q=0, i.e. up to
        # 25% more real power than the stated nameplate, in every case.)
        gen["pmax"] =  pv_pu * ones(n_ph)
        gen["pmin"] =  zeros(n_ph)
        gen["qmax"] =  s_pu * ones(n_ph)
        gen["qmin"] = -s_pu * ones(n_ph)
        gen["cost"] = [pv_cost 0.0]

        # New: apparent-power rating for the joint capability constraint,
        # AND the true active-power nameplate rating for reporting — pmax
        # is now a loose s_pu ceiling, so it can no longer double as "the
        # kW capacity" the way it did in the old box-bound version.
        gen["s_rating"] = s_pu
        gen["p_rating"] = pv_pu

        # New: VVC droop configuration (used by constraints.jl only if
        # enable_vvc=true; harmless to carry even when false).
        gen["vvc"]      = enable_vvc
        gen["vvc_v1"]   = v1_pu
        gen["vvc_v2"]   = v2_pu
        gen["vvc_v3"]   = v3_pu
        gen["vvc_v4"]   = v4_pu
        gen["vvc_qbar"] = vvc_qbar

        # Defensive: trim any other length-3 template field (e.g. vg/pg/qg
        # starting points) down to n_ph so the gen record stays internally
        # consistent. PMD versions vary in which fields exist here.
        for (k, v) in gen
            if v isa Vector{<:Real} && length(v) == 3 && n_ph != 3
                gen[k] = v[phase_conns]
            end
        end

        push!(pv_ids, gen_id)
        n_ph == 1 ? (n_ph1 += 1) : (n_ph3 += 1)
    end

    vvc_str = enable_vvc ? "VVC droop ON (V1-V4=$(vvc_v1_v)/$(vvc_v2_v)/$(vvc_v3_v)/$(vvc_v4_v) V, qbar=$(vvc_qbar))" : "VVC OFF (Q free in OPF)"
    println("  PV: $(length(pv_ids)) units  pv_kw=$(pv_kw) kW/unit  s_scale=$(s_scale)  spacing=$(spacing)  cost=$(pv_cost)")
    println("      single-phase: $n_ph1   three-phase: $n_ph3   skipped (no phase match): $n_skip   $(vvc_str)")
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
function solve_and_report(data_math, label; objective_choice="cost", sbase_kva=SBASE_KVA, vbase_kv=VBASE_KV)
    # ══════════════════════════════════════════════════════════════
    # REWRITTEN: this now matches the REAL working pattern (confirmed
    # against a known-good script) — ref/model built manually, custom
    # physics pulled in via include(), NOT PMD.instantiate_mc_model.
    # That earlier path never touched core/constraints.jl at all, which
    # is why the VVC droop code silently had zero effect for this whole
    # debugging thread — it was never in the call graph.
    # ══════════════════════════════════════════════════════════════
    # include() always runs in GLOBAL scope, never the calling function's
    # local scope — variables.jl's very first line references `ref`, so it
    # must already be a global by the time include() runs, not just a
    # local inside this function. Same for `model`.
    global objective = objective_choice   # renamed the keyword arg specifically
                                    # to avoid this: `global objective = objective`
                                    # (same name on both sides) is a real Julia
                                    # scoping trap — once `global objective`
                                    # appears anywhere in a function, EVERY bare
                                    # occurrence of that name in the function is
                                    # treated as the global, including the RHS
                                    # of that same line. That made the RHS read
                                    # the (unassigned) global instead of the
                                    # keyword argument — exactly the error seen.
    global ref = IM.build_ref(data_math, PMD.ref_add_core!, PMD._pmd_global_keys, PMD.pmd_it_name)[:it][:pmd][:nw][0]
    global model = JuMP.Model(ipopt_solver)
    include("./core/variables.jl")
    include("./core/constraints_VVC.jl")
    include("./core/objectives_FIXED.jl")

    println("    solving...")
    t0 = time()
    JuMP.optimize!(model)
    println("    solve took $(round(time()-t0, digits=1))s")
    status = JuMP.termination_status(model)
    println("\n  [$label]  →  $(status)")

    v_min = NaN;  v_mean = NaN;  v_max = NaN
    n_over = 0;   n_under = 0
    pv_util = NaN;  pv_curtail = NaN
    pv_output = NaN;  pv_capacity = NaN
    st_q_total = NaN;  st_q_cap = NaN
    st_util_net = NaN; st_util_gross = NaN
    vm_per_bus = Dict{String, Float64}()

    if status in [JuMP.LOCALLY_SOLVED, JuMP.ALMOST_LOCALLY_SOLVED, JuMP.OPTIMAL]

        # ── Voltage profile ─────────────────────────────────────
        # NOTE: ref[:bus] keys are Int (IM.build_ref convention), not the
        # String keys data_math uses. bus["bus_type"] is read straight off
        # ref[:bus][b] — IM's ref-building generally carries this field
        # through from data_math unchanged; flag if this errors, as it
        # would mean bus_type isn't propagated and needs pulling from
        # data_math["bus"][string(b)] instead.
        vm_all     = Float64[]
        bus_labels = String[]

        for (b, bus) in ref[:bus]
            bus["bus_type"] == 3 && continue
            bus_max = 0.0
            for p in 1:3
                vr_val = JuMP.value(vr[p, b])
                vi_val = JuMP.value(vi[p, b])
                vm_val = abs(vr_val + im * vi_val)
                push!(vm_all,     vm_val)
                push!(bus_labels, string(b))
                bus_max = max(bus_max, vm_val)
            end
            vm_per_bus[string(b)] = bus_max
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
        # Metadata (type, p_rating, s_rating, vvc flags) is pulled from
        # data_math["gen"] (String keys, definitely carries our custom
        # fields since add_pv! set them there directly) rather than
        # ref[:gen], since it's unverified whether IM.build_ref preserves
        # arbitrary custom keys through to ref. Solved VALUES (pg/qg) come
        # from the JuMP variable containers via ref's Int keys, since
        # that's what pg/qg are actually indexed by.
        pv_gen_ids = [parse(Int, i) for (i, g) in data_math["gen"] if get(g, "type", "") == "PV"]
        if !isempty(pv_gen_ids)
            pg_vals = Float64[]
            pg_rated = Float64[]
            for gid in pv_gen_ids
                g = data_math["gen"][string(gid)]
                # pg is zero-filled (not a real variable) on unused phase
                # slots per variables.jl's construction, so summing 1:3
                # unconditionally is safe regardless of this unit's n_ph.
                push!(pg_vals, sum(JuMP.value(pg[p, gid]) for p in 1:3))
                push!(pg_rated, get(g, "p_rating", sum(g["pmax"])))
            end
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
        st_gen_ids = [parse(Int, i) for (i, g) in data_math["gen"] if get(g, "type", "") == "STATCOM"]
        if !isempty(st_gen_ids)
            qg_vecs = Vector{Float64}[]
            qg_rated = Float64[]
            for gid in st_gen_ids
                g = data_math["gen"][string(gid)]
                push!(qg_vecs, [JuMP.value(qg[p, gid]) for p in 1:3])
                push!(qg_rated, sum(g["qmax"]))
            end

            qg_signed_sum = sum(sum(q) for q in qg_vecs)          # NET: signs can cancel
            qg_abs_sum    = sum(sum(abs.(q)) for q in qg_vecs)    # GROSS: no cancellation

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
# DEV_MODE = true: fast sanity-check run — fewer units (wider spacing),
# one PV kW level, Case 3 skipped entirely. Should finish in well under a
# minute. Flip to false once you've confirmed everything looks sane, for
# the real full sweep (52 units, all kW levels, full STATCOM rating sweep).
DEV_MODE = false
DEV_SPACING   = DEV_MODE ? 10 : 1
DEV_KW_LEVELS = DEV_MODE ? [5.0] : [1.0, 3.0, 5.0, 7.0, 10.0]

do_case1 = true
do_case2 = true
do_case3 = !DEV_MODE

# NOTE: not `const` — these get re-run in the same Julia session often while
# iterating, and `const` reassignment (even to an equal value) triggers a
# noisy redefinition warning for array-valued constants specifically (Julia
# allows silent redefinition of scalar consts to equal values, but not
# arrays, since each literal is a fresh, distinct object).
PV_COST      = -1000.0
STATCOM_COST =  1.0

# Case 2: sweep of realistic single-phase inverter sizes (kW)
PV_KW_LEVELS = [1.0, 3.0, 5.0, 7.0, 10.0]

# Case 3: fixed PV size for hosting-capacity stress test, chosen to sit
# inside typical UK single-phase residential inverter range (≤5kW is most
# common; higher values push into stress-test territory deliberately)
CASE3_PV_KW = 20

# Case 3: STATCOM rating sweep, kVAr TOTAL nameplate per unit
STATCOM_RATINGS_KVAR = [1.0, 5.0, 10.0, 20.0, 30.0, 50.0, 70.0, 100.0, 200.0]

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
## Case 2: PV penetration sweep — three-way comparison at each kW level:
##   (a) UPF        — Q forced to 0 (no VVC, no free dispatch)
##   (b) free-Q      — Q a free OPF variable (old Option-A "ceiling")
##   (c) VVC droop   — Q pinned by local per-phase Volt-Var droop (realistic)
##   NOTE: (b) and (c) both require the constraints.jl patch (S^2 capability
##   constraint) to behave correctly; (c) additionally requires the VVC
##   droop block. Without the patch, (b)/(c) will run using only the loose
##   box bounds and will NOT show the intended P/Q competition or droop
##   behavior — see constraints_patch_vvc.jl.
# ───────────────────────────────────────────────────────────────
if do_case2
    println("\n" * "="^55)
    println(" CASE 2: PV Penetration Sweep (no bounds)")
    println(" Single-phase PV, phase-matched to colocated load")
    println(" UPF vs free-Q (Option A ceiling) vs VVC droop (realistic)")
    println("="^55)

    for pv_kw in DEV_KW_LEVELS
        println("\n  ── pv_kw = $(pv_kw) kW/unit ──")

        println("    (a) UPF — Q forced to 0")
        dm_upf = load_base_network(data_path, enforce_bounds=false)
        pv_ids_upf = add_pv!(dm_upf; pv_kw=pv_kw, s_scale=1.25, spacing=DEV_SPACING,
                              pv_cost=PV_COST, enable_vvc=false)
        for id in pv_ids_upf
            dm_upf["gen"][id]["qmax"] = zeros(length(dm_upf["gen"][id]["qmax"]))
            dm_upf["gen"][id]["qmin"] = zeros(length(dm_upf["gen"][id]["qmin"]))
        end
        solve_and_report(dm_upf, "PV $(pv_kw) kW/unit  UPF [Q=0]")

        println("    (b) free-Q — Option A ceiling (centrally-optimized Q)")
        dm_free = load_base_network(data_path, enforce_bounds=false)
        add_pv!(dm_free; pv_kw=pv_kw, s_scale=1.25, spacing=DEV_SPACING,
                pv_cost=PV_COST, enable_vvc=false)
        solve_and_report(dm_free, "PV $(pv_kw) kW/unit  free-Q [Option A]")

        println("    (c) VVC droop — local, autonomous, per-phase")
        dm_vvc = load_base_network(data_path, enforce_bounds=false)
        add_pv!(dm_vvc; pv_kw=pv_kw, s_scale=1.25, spacing=DEV_SPACING,
                pv_cost=PV_COST, enable_vvc=true)
        solve_and_report(dm_vvc, "PV $(pv_kw) kW/unit  VVC droop")
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
    println(" Builder: ref/JuMP.Model + core/{variables,constraints,objectives}.jl  |  Objective: cost")
    println("="^55)

    println("\n  ── 3a: PV only (VVC droop) — hosting capacity baseline ──")
    dm3a = load_base_network(data_path, enforce_bounds=true)
    add_pv!(dm3a; pv_kw=CASE3_PV_KW, s_scale=1.25, spacing=1, pv_cost=PV_COST, enable_vvc=true)
    r3a  = solve_and_report(dm3a, "PV only (VVC droop)  [no STATCOM]")

    println("\n  ── 3b: STATCOM rating sweep (52 STATCOMs, spacing=1) ──")
    rating_results = NamedTuple[]
    r3b_last = nothing
    for rating in STATCOM_RATINGS_KVAR
        println("\n  ── STATCOM rating = $(rating) kVAr/unit  (52 units) ──")
        local dm3b = load_base_network(data_path, enforce_bounds=true)
        add_pv!(dm3b;       pv_kw=CASE3_PV_KW, s_scale=1.25, spacing=1, pv_cost=PV_COST, enable_vvc=true)
        add_statcoms!(dm3b; rating_kvar=rating, spacing=1,   statcom_cost=STATCOM_COST)
        r = solve_and_report(dm3b, "PV+STATCOM $(rating) kVAr/unit  52 units")
        push!(rating_results, (rating_kvar=rating, util=r.pv_util, baseline=r3a.pv_util))
        # NOTE: this loop runs at top-level script scope, where Julia's soft-scope
        # rules treat a bare `r3b_last = r` as creating a NEW local shadowing the
        # outer variable (silently — you only see it as a one-line warning at
        # runtime), leaving the outer r3b_last stuck at `nothing`. `global` forces
        # it to update the actual outer variable instead.
        global r3b_last = r
    end
    plot_rating_sweep(rating_results)

    if !isnothing(r3b_last)
        plot_voltage_profile(r3a, r3b_last, CASE3_PV_KW, STATCOM_RATINGS_KVAR[end])
    end

    table_rows = vcat(
        [(label="No STATCOM (baseline)", util=r3a.pv_util, curtail=r3a.pv_curtail)],
        [(label="STATCOM $(Int(r.rating_kvar)) kVAr (52 units)", util=r.util, curtail=100.0 - r.util)
         for r in rating_results]
    )
    print_delta_hc_table(table_rows)
end

println("\n  Plots saved to ./plots/")
println("  Done.")