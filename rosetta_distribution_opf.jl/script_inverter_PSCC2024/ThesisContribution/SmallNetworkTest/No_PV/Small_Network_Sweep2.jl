#=
==============================================================================
SMALL NETWORK SWEEP 2 -- 4-bus, 4-wire, MEN, 30 kVA STATCOM
==============================================================================
Successor to Small_Network_Sweep.jl (3-bus). Same validated code path as the
908-bus feeder -- load_base_network / add_men! / build_and_solve are carried
over so that what is validated here transfers there.

WHAT CHANGED, AND WHY
---------------------------------------------------------------------------
1. NETWORK: 4 buses (sourcebus -> mid1 -> mid2 -> endbus), 21 lots at the
   Energex LV design ADMD of 4.5 kVA, 185 mm^2 Al mains, 550 m total. Every
   parameter is sourced rather than chosen for convenience -- see the header
   of 4_Bus_Small.txt. The device rating (30 kVA) is set from the class of
   commercially available LV phase-balancing units, and the network is sized
   so that a device of that class is genuinely challenged rather than either
   idle or hopeless.

2. STRESS AXIS: PHASE ALLOCATION of the 9 lots at endbus, from 3/3/3 to
   7/1/1, with total demand held constant at 89.8 kW. This is how real LV
   feeders become unbalanced -- single-phase lots connected to whichever
   phase was convenient. It replaces the previous ad-hoc "unbalance index",
   which was a quantity of the author's own definition with no external
   referent. If more stress is needed, LENGTHEN THE LINES; do not inflate
   loads.

3. SCENARIOS are now named for what the converter can actually do, and
   include the two cases that were missing:
       A   no STATCOM
       B   Q-only, bidirectional         (reactive support, both signs)
       Bc  Q-only, CAPACITIVE ONLY       (single sign -- see note below)
       C   P-only                        (inter-phase real power, Q pinned 0)
       D   P+Q                           (full capability)
   C isolates how much of the benefit comes from real power alone. Bc exists
   because bidirectional reactive power lets devices circulate current that
   substation protection never sees: a downstream conductor can be loaded
   without any upstream relay registering it. Restricting every device to
   one sign removes that failure mode. Comparing Bc against B prices it.

4. STATCOM PLACEMENT IS EXPLICIT. The previous add_statcoms! picked a bus
   from the load dictionary's iteration order. With load on three buses that
   would silently place the unit wherever the dictionary happened to sort.
   add_statcom! now takes a bus NAME. Block 4 uses this to compare siting.

5. gen["connections"] IS SET EXPLICITLY. The previous version relied on the
   deepcopy from gen 1 to supply it -- the same latent defect that add_pv!
   carries a BUGFIX comment about. It never bit because every unit was
   3-phase, but it was never correct.

6. BOTH LIMIT CASES ARE RUN. Block 1 relaxes the statutory limits and draws
   them as reference lines: with the limits enforced, a breach cannot be
   RETURNED by the solver at all, only an infeasibility, so compliance would
   be assumed rather than demonstrated. Block 3 then enforces them, which is
   what an operator bound by those limits actually gets. Both are true; they
   answer different questions.

7. VUF CAP IS A RATIO AND IS GUARDED. constraints_unified.jl previously
   carried an unconditional, hard-coded |V2| <= 0.02 on every bus -- an
   absolute negative-sequence bound, not the ratio |V2|/|V1| that the
   standards define, and invisible from any calling script. It is now
   guarded on VUF_CAP with the original constant as the fallback, so the
   908-bus scripts are unaffected until they opt in. Required patch:

       if @isdefined(VUF_CAP)
           if isfinite(VUF_CAP)
               JuMP.@constraint(model, v_neg_seq_real^2 + v_neg_seq_imag^2
                            <= VUF_CAP^2 * (v_pos_seq_real^2 + v_pos_seq_imag^2))
           end
       else
           JuMP.@constraint(model, v_neg_seq_real^2 + v_neg_seq_imag^2 <= 0.02^2)
       end

   (v_pos_seq_real/imag are built with Tre[2,:]/Tim[2,:], exactly as the
   negative-sequence pair are built with Tre[3,:]/Tim[3,:].)

8. VOLTAGE BOUNDS ARE REFERENCED TO NEUTRAL. constraints_unified.jl enforces
   vr[t]^2 + vi[t]^2 >= vmin[t]^2, which is phase-to-GROUND. On a MEN network
   the quantity a customer's appliance sees is phase-to-NEUTRAL, and the two
   differ by exactly the neutral displacement. Note the contrast, which is
   worth stating explicitly in the methodology: the reference is DECISIVE for
   voltage magnitudes and PROVABLY IRRELEVANT for VUF, because a common-mode
   shift is purely zero-sequence and 1 + a + a^2 = 0 cancels it out of
   |V2|/|V1|.

9. LOSSES are reported everywhere, including the neutral-conductor share.
   Neutral I^2R is invisible in a 3-wire Kron-reduced model and is created
   entirely by unbalance, which makes it the sharpest performance metric
   available here -- and the only one denominated in kW.

A NOTE ON "cost": with the source generator carrying a positive cost and the
STATCOM's own cost term vanishing (P is pinned in B/Bc, and sum(P) = 0 in
C/D), minimising cost reduces to minimising source real power, i.e. load plus
losses. Load is fixed. The "cost" rows are therefore a LOSS MINIMISATION
case, not a degenerate one, and should be labelled as such.

A NOTE ON "loss": objectives_FIXED.jl's "loss" branch indexes br_r[c] with a
single index, which reads a 4x4 matrix column-major and so picks r11, r21,
r31, r41 -- the self resistance of conductor 1 and then three MUTUALS. It
also weights by |Z| rather than r, and ignores the off-diagonal cross terms.
It is not used here; "cost" is the better loss objective until that is fixed.

Run with:
    julia --project=<thesis root> Small_Network_Sweep2.jl
==============================================================================
=#

using Logging
Logging.disable_logging(Logging.Warn)

using Pkg

# --- Path anchoring (depth-independent) ---------------------------------
# Locations are found by SEARCHING UPWARD for marker contents rather than by
# counting ".." segments, so moving this script to a different depth needs no
# edits. THESIS_ROOT deliberately takes the OUTERMOST Project.toml: both the
# package and the thesis directory have one, and a nearest-match search would
# activate the wrong environment.
const SCRIPT_DIR = @__DIR__

function _find_up(start::AbstractString, pred::Function; outermost::Bool=false, maxdepth::Int=12)
    d, hit = start, nothing
    for _ in 1:maxdepth
        if pred(d)
            hit = d
            outermost || return hit
        end
        parent = dirname(d)
        parent == d && break
        d = parent
    end
    return hit
end

const PSCC_DIR    = _find_up(SCRIPT_DIR, d -> isfile(joinpath(d, "core", "variables.jl")))
const RPMD_ROOT   = _find_up(SCRIPT_DIR, d -> isdir(joinpath(d, "data")) && isdir(joinpath(d, "src")))
const THESIS_ROOT = _find_up(SCRIPT_DIR, d -> isfile(joinpath(d, "Project.toml")); outermost=true)

isnothing(PSCC_DIR)    && error("could not find core/variables.jl above $SCRIPT_DIR")
isnothing(RPMD_ROOT)   && error("could not find a directory with both data/ and src/ above $SCRIPT_DIR")
isnothing(THESIS_ROOT) && error("could not find a Project.toml above $SCRIPT_DIR")

const CORE_DIR = joinpath(PSCC_DIR, "core")
println("  core   : $CORE_DIR")
println("  package: $RPMD_ROOT")
println("  project: $THESIS_ROOT")

Pkg.activate(THESIS_ROOT)
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

ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "sb" => "yes")

# <-- set to wherever 4_Bus_Small.txt lives. Both options given; pick one.
data_path = joinpath(RPMD_ROOT, "data", "4_Bus_Small.dss")
# data_path = joinpath(RPMD_ROOT, "data", "4_Bus_Small.txt")
outdir    = joinpath(THESIS_ROOT, "results", "small_network_2")

isfile(data_path) || error("network file not found: $data_path")

# =======================================================================
# NETWORK AND DEVICE CONFIGURATION
# =======================================================================

const ADMD_KVA   = 4.5           # Energex LV design ADMD per residential lot
const LOAD_PF    = 0.95
const KW_PER_LOT = ADMD_KVA * LOAD_PF        # 4.275 kW

# endbus phase allocation, in LOTS. Total is 9 in every case, so total demand
# is constant at 89.8 kW and only the allocation changes.
const ALLOCATIONS = [(3,3,3), (4,3,2), (5,3,1), (6,2,1), (7,1,1)]
const WORST_ALLOC = ALLOCATIONS[end]

const LOAD_BUS_NAME = "endbus"   # bus whose allocation is swept
const STATCOM_BUS   = "endbus"   # default siting; block 4 varies it

# Device rating: TOTAL apparent power, split evenly across three legs.
# 30 kVA = 10 kVA/leg, the class of commercially available LV phase
# balancers. For the Q-only scenarios this is effectively a kvar rating.
const STATCOM_KVA = 30.0
const N_STATCOMS  = 1

# MEN bonding. 10 ohm is the per-installation assumption; real electrode
# resistance spans roughly 1 ohm (old metallic water-pipe bonding) to 100
# ohm+ (single rod, dry sandy soil), so this is the parameter the results
# are most sensitive to and the obvious candidate for a sensitivity sweep.
# MEN_PER_LOT divides it by the number of lots at each bus, since each
# installation carries its own bond in parallel.
const MEN_OHM     = 10.0
const MEN_PER_LOT = true
const LOTS_AT_BUS = Dict("mid1" => 6, "mid2" => 6, "endbus" => 9)

# Statutory limits. Block 1 relaxes them and draws them as references; block
# 3 enforces them. VUF_LIMIT_REF is never enforced anywhere -- it only marks
# figures and sets the breach flag.
const VMIN_PU       = 0.90
const VMAX_PU       = 1.10
const VUF_LIMIT_REF = 0.02
const NEUTRAL_REF   = true

# Rating sweep (block 2). Full correction of the 7/1/1 allocation needs
# roughly 50 kVA of leg capacity, so the knee sits inside this range.
const RATING_KVA_LEVELS = [5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 60.0, 90.0]

# Scenario table: (code, label, mode, objective)
const SCENARIOS = [
    ("A",  "no STATCOM",          :none,  "cost"),
    ("B",  "Q-only",              :qonly, "VUF"),
    ("Bc", "Q-only capacitive",  :qcap,  "VUF"),
    ("C",  "P-only",              :ponly, "VUF"),
    ("D",  "P+Q",                 :pq,    "VUF"),
]
# The same devices under the cost objective, which is loss minimisation here
# (see the header note). Scenario A is excluded because it already runs under
# "cost" above and, having no controllable device, would give an identical row.
const SCENARIOS_COST = [(c, l, m, "cost") for (c, l, m, _) in SCENARIOS if m != :none]

# -----------------------------------------------------------------------
# kW/kVAr <-> pu
# -----------------------------------------------------------------------
kw_to_pu(p_kw, sbase_kva)     = p_kw / sbase_kva
kvar_to_pu(q_kvar, sbase_kva) = q_kvar / sbase_kva
pu_to_kw(p_pu, sbase_kva)     = p_pu * sbase_kva
pu_to_kvar(q_pu, sbase_kva)   = q_pu * sbase_kva

# -----------------------------------------------------------------------
# MEN -- unchanged in mechanism from the validated version, with one
# addition: lots_at_bus divides the per-installation resistance by the
# number of installations at that bus, since the bonds are in parallel.
# -----------------------------------------------------------------------
function add_men!(data_eng; r_lot_ohm=10.0, lots_at_bus=nothing, verbose=false)
    haskey(data_eng, "shunt") || (data_eng["shunt"] = Dict{String,Any}())

    men_buses = String[]
    for (_, load) in data_eng["load"]
        b = load["bus"]
        4 in load["connections"] || continue
        b in men_buses || push!(men_buses, b)
    end
    sort!(men_buses)

    for b in men_buses
        n = isnothing(lots_at_bus) ? 1 : get(lots_at_bus, b, 1)
        r = r_lot_ohm / n
        data_eng["shunt"]["men_$b"] = Dict{String,Any}(
            "source_id"     => "reactor.men_$b",
            "status"        => PMD.ENABLED,
            "model"         => PMD.REACTOR,
            "connections"   => [4],
            "dispatchable"  => PMD.NO,
            "gs"            => fill(1.0 / r, 1, 1),
            "bs"            => zeros(1, 1),
            "bus"           => b,
            "configuration" => PMD.WYE,
        )
        verbose && println("    MEN $b: $n lots -> $(round(r, digits=3)) Ω")
    end
    return men_buses
end

# -----------------------------------------------------------------------
# Network loader
# -----------------------------------------------------------------------
function load_base_network(data_path; enforce_bounds=true, sbase_kva=100,
                           men_ohm=nothing, neutral_ref=true, lots_at_bus=nothing)
    data_eng = PMD.parse_file(data_path, transformations=[PMD.transform_loops!])
    isnothing(men_ohm) || add_men!(data_eng; r_lot_ohm=men_ohm, lots_at_bus=lots_at_bus)

    data_math = PMD.transform_data_model(
        data_eng, multinetwork=false, kron_reduce=false, phase_project=false, make_pu=false)
    PMD.make_per_unit!(data_math; sbase=sbase_kva)
    PMD.add_start_vrvi!(data_math)

    for (i, bus) in data_math["bus"]
        if bus["bus_type"] == 3
            bus["vmin"] = zeros(4)
            bus["vmax"] = 10.0 * ones(4)
        elseif neutral_ref
            # Phase bounds switched OFF on the built-in phase-to-ground path:
            # vmin=0 and vmax=Inf are read as "skip" by that file's own
            # nonzero_vmin_terminals / nonInf_vmax_terminals filters.
            # add_neutral_referenced_bounds! then imposes the correct ones.
            # Terminal 4 keeps its built-in vmax, because neutral-to-GROUND
            # is genuinely what a neutral displacement cap should measure.
            bus["vmin"] = [0.00, 0.00, 0.00, 0.00]
            bus["vmax"] = [Inf,  Inf,  Inf,  enforce_bounds ? 0.20 : 1.50]
        elseif enforce_bounds
            bus["vmin"] = [VMIN_PU, VMIN_PU, VMIN_PU, 0.00]
            bus["vmax"] = [VMAX_PU, VMAX_PU, VMAX_PU, 0.20]
        else
            bus["vmin"] = [0.50, 0.50, 0.50, 0.00]
            bus["vmax"] = [1.50, 1.50, 1.50, 1.50]
        end
    end

    for (_, gen) in data_math["gen"]
        gen["pmax"] =  [1e4, 1e4, 1e4]
        gen["pmin"] = -[1e4, 1e4, 1e4]
        gen["qmax"] =  [1e4, 1e4, 1e4]
        gen["qmin"] = -[1e4, 1e4, 1e4]
    end

    return data_math
end

get_sbase_kva(dm) = dm["settings"]["sbase"] * dm["settings"]["power_scale_factor"] / 1000

function total_load_kw_kvar(dm, sbase_kva)
    p = sum(sum(l["pd"]) for (_, l) in dm["load"])
    q = sum(sum(l["qd"]) for (_, l) in dm["load"])
    return pu_to_kw(p, sbase_kva), pu_to_kvar(q, sbase_kva)
end

# -----------------------------------------------------------------------
# PHASE ALLOCATION -- the stress mechanism.
#
# Overwrites pd/qd on the three single-phase loads at ONE named bus, leaving
# every other bus untouched. Loads are matched by (bus, connection), never by
# name or dictionary order: dict iteration order is not guaranteed, and
# silently assigning phase B's demand to phase C would produce entirely
# plausible wrong figures. Throws rather than guessing.
# -----------------------------------------------------------------------
function set_phase_allocation!(dm, sbase_kva; bus_name, lots_per_phase, pf=LOAD_PF)
    bus_id = get(dm["bus_lookup"], bus_name, nothing)
    isnothing(bus_id) && error("bus '$bus_name' not found; have: $(sort(collect(keys(dm["bus_lookup"]))))")

    tanphi = tan(acos(pf))
    seen = Dict(1 => false, 2 => false, 3 => false)

    for (_, load) in dm["load"]
        load["load_bus"] == bus_id || continue
        ph = filter(c -> c != 4, load["connections"])
        length(ph) == 1 || error("expected single-phase loads at $bus_name, got connections $(load["connections"])")
        p = ph[1]
        haskey(seen, p) || error("load on unexpected phase $p at $bus_name")
        seen[p] && error("more than one load on phase $p at $bus_name -- one load per phase is assumed")

        p_kw = lots_per_phase[p] * KW_PER_LOT
        load["pd"] = [kw_to_pu(p_kw, sbase_kva)]
        load["qd"] = [kw_to_pu(p_kw * tanphi, sbase_kva)]
        seen[p] = true
    end
    all(values(seen)) || error("not every phase at $bus_name received a load: $seen")
    return nothing
end

alloc_label(a) = "$(a[1])/$(a[2])/$(a[3])"
alloc_spread_kw(a) = (maximum(a) - minimum(a)) * KW_PER_LOT

# -----------------------------------------------------------------------
# STATCOM PLACEMENT -- explicit bus, explicit connections, four modes.
#
#   :qonly  reactive only, both signs   -> original box path (no s_rated)
#   :qcap   reactive only, ONE sign     -> box path, qmin = 0
#   :ponly  inter-phase real power only -> S^2 path, Q pinned to 0
#   :pq     full capability             -> S^2 path
#
# The S^2 path in constraints_unified.jl applies, per leg,
#     (pg/s_rated)^2 + (qg/s_rated)^2 <= 1
# plus the DC-bus coupling sum(pg) == p_loss. :ponly keeps that coupling --
# real power is still only MOVED between phases, never imported.
# -----------------------------------------------------------------------
function add_statcom!(dm, sbase_kva; bus_name, rating_kva_total, mode::Symbol)
    bus_id = get(dm["bus_lookup"], bus_name, nothing)
    isnothing(bus_id) && error("STATCOM bus '$bus_name' not found")

    s_leg_pu = kvar_to_pu(rating_kva_total, sbase_kva) / 3

    gen_id = string(length(dm["gen"]) + 1)
    dm["gen"][gen_id] = deepcopy(dm["gen"]["1"])
    gen = dm["gen"][gen_id]

    gen["gen_bus"] = bus_id
    gen["type"]    = "STATCOM"
    gen["name"]    = "statcom_$(bus_name)"
    gen["cost"]    = [0.001, 0.0]
    # Set explicitly rather than inherited from the gen-1 deepcopy.
    gen["connections"]   = [1, 2, 3, 4]
    gen["configuration"] = PMD.WYE

    if mode == :qonly
        gen["pmax"] = zeros(3);            gen["pmin"] = zeros(3)
        gen["qmax"] = fill(s_leg_pu, 3);   gen["qmin"] = -fill(s_leg_pu, 3)
        gen["statcom_p_exchange"] = false
    elseif mode == :qcap
        gen["pmax"] = zeros(3);            gen["pmin"] = zeros(3)
        gen["qmax"] = fill(s_leg_pu, 3);   gen["qmin"] = zeros(3)   # single sign
        gen["statcom_p_exchange"] = false
    elseif mode == :ponly
        gen["pmax"] =  fill(s_leg_pu, 3);  gen["pmin"] = -fill(s_leg_pu, 3)
        gen["qmax"] =  zeros(3);           gen["qmin"] =  zeros(3)
        gen["statcom_p_exchange"] = true
        gen["s_rated"] = fill(s_leg_pu, 3)
        gen["p_loss"]  = 0.0
    elseif mode == :pq
        gen["pmax"] =  fill(s_leg_pu, 3);  gen["pmin"] = -fill(s_leg_pu, 3)
        gen["qmax"] =  fill(s_leg_pu, 3);  gen["qmin"] = -fill(s_leg_pu, 3)
        gen["statcom_p_exchange"] = true
        gen["s_rated"] = fill(s_leg_pu, 3)
        gen["p_loss"]  = 0.0
    else
        error("unknown STATCOM mode :$mode")
    end

    return [parse(Int, gen_id)]
end

# -----------------------------------------------------------------------
# Neutral-referenced phase voltage bounds.
#
# Same nonconvex quadratic form as constraints_unified.jl -- only the
# reference point moves:
#     (vr[p] - vr[4])^2 + (vi[p] - vi[4])^2  in  [vmin^2, vmax^2]
#
# constraints_unified.jl is deliberately NOT modified for this: it is shared
# with the 908-bus runs, and changing it there would silently invalidate
# results already produced. This is additive and opt-in.
# -----------------------------------------------------------------------
function add_neutral_referenced_bounds!(; vmin=VMIN_PU, vmax=VMAX_PU)
    n = 0
    for (i, bus) in ref[:bus]
        is_real_bus(bus) || continue
        terms = bus["terminals"]
        4 in terms || continue
        for p in 1:3
            p in terms || continue
            dvr = vr[p, i] - vr[4, i]
            dvi = vi[p, i] - vi[4, i]
            JuMP.@constraint(model, dvr^2 + dvi^2 >= vmin^2)
            JuMP.@constraint(model, dvr^2 + dvi^2 <= vmax^2)
            n += 1
        end
    end
    return n
end

# -----------------------------------------------------------------------
# Build + solve.
#
# The raw-include pipeline is kept exactly as validated. The only additions
# are the two globals read by the guarded constraint and the opt-in bounds.
# VUF_CAP must be set on EVERY call, not once: it is a global that persists
# across scenarios in one session, so a stale value would carry silently.
# -----------------------------------------------------------------------
function build_and_solve(dm; objective_mode::String="cost", neutral_ref::Bool=NEUTRAL_REF,
                         bounds_active::Bool=false, vmin=VMIN_PU, vmax=VMAX_PU,
                         vuf_cap=Inf)
    pm = PMD.instantiate_mc_model(
        dm, PMD.IVRENPowerModel, (pm) -> nothing;
        jump_model = JuMP.Model(ipopt_solver), ref_extensions = Function[])

    global ref   = pm.ref[:it][:pmd][:nw][0]
    global model = pm.model

    global VUF_CAP = vuf_cap

    include(joinpath(CORE_DIR, "variables.jl"))
    include(joinpath(CORE_DIR, "constraints_unified.jl"))

    neutral_ref && bounds_active && add_neutral_referenced_bounds!(vmin=vmin, vmax=vmax)

    global objective = objective_mode
    global objective_target_gens = [i for (i, g) in ref[:gen] if get(g, "type", "") == "STATCOM"]
    include(joinpath(CORE_DIR, "objectives_FIXED.jl"))

    JuMP.optimize!(model)
    return JuMP.termination_status(model)
end

solved(st) = st in [JuMP.LOCALLY_SOLVED, JuMP.OPTIMAL, JuMP.ALMOST_LOCALLY_SOLVED]

# =======================================================================
# REPORTING
# =======================================================================

# The virtual source bus is bus_type 3 and is excluded by that test alone;
# the name check is belt and braces. sourcebus is an ordinary bus and IS
# included, so vmax ~ 1.00 means "no feeder bus exceeds the substation".
is_virtual(bus)  = startswith(get(bus, "name", ""), "_virtual")
is_real_bus(bus) = !is_virtual(bus) && get(bus, "bus_type", 1) != 3

const ALPHA = exp(im * 2/3 * pi)
const TSEQ  = 1/3 * [1 1 1; 1 ALPHA ALPHA^2; 1 ALPHA^2 ALPHA]

# VUF = |V2|/|V1|. The reference point is PROVABLY IRRELEVANT here: a common
# neutral shift is purely zero-sequence, and 1 + a + a^2 = 0 removes it from
# both V1 and V2. Confirmed empirically on the 3-bus case -- the neutral- and
# ground-referenced values agreed to three decimals in all 30 rows. Only one
# column is therefore kept.
function worst_case_vuf()
    worst, worst_bus = 0.0, nothing
    for (i, bus) in ref[:bus]
        is_virtual(bus) && continue
        t = bus["terminals"]
        (1 in t && 2 in t && 3 in t) || continue
        vph  = [JuMP.value(vr[p, i]) + im * JuMP.value(vi[p, i]) for p in 1:3]
        v012 = TSEQ * vph
        vpos, vneg = abs(v012[2]), abs(v012[3])
        if vpos > 1e-6 && vneg / vpos > worst
            worst, worst_bus = vneg / vpos, i
        end
    end
    return worst, worst_bus
end

# Absolute |V2|, the quantity the ORIGINAL hard-coded cap bounded. Kept as a
# diagnostic so the difference between that cap and a true VUF ratio stays
# visible.
function worst_neg_seq()
    worst, worst_bus = 0.0, nothing
    for (i, bus) in ref[:bus]
        is_virtual(bus) && continue
        t = bus["terminals"]
        (1 in t && 2 in t && 3 in t) || continue
        v2 = abs((TSEQ * [JuMP.value(vr[p, i]) + im * JuMP.value(vi[p, i]) for p in 1:3])[3])
        v2 > worst && ((worst, worst_bus) = (v2, i))
    end
    return worst, worst_bus
end

# Phase-to-NEUTRAL by default: what a customer's appliance actually sees.
function worst_case_vmag(; ref_neutral::Bool=true)
    min_v, min_bus, min_ph = Inf, nothing, nothing
    max_v, max_bus, max_ph = -Inf, nothing, nothing
    for (i, bus) in ref[:bus]
        is_real_bus(bus) || continue
        t = bus["terminals"]
        vr_n = (ref_neutral && 4 in t) ? JuMP.value(vr[4, i]) : 0.0
        vi_n = (ref_neutral && 4 in t) ? JuMP.value(vi[4, i]) : 0.0
        for p in 1:3
            p in t || continue
            vm = sqrt((JuMP.value(vr[p, i]) - vr_n)^2 + (JuMP.value(vi[p, i]) - vi_n)^2)
            vm < min_v && ((min_v, min_bus, min_ph) = (vm, i, p))
            vm > max_v && ((max_v, max_bus, max_ph) = (vm, i, p))
        end
    end
    return (min_v=min_v, min_bus=min_bus, min_phase=min_ph,
            max_v=max_v, max_bus=max_bus, max_phase=max_ph)
end

function worst_case_nev()
    worst, worst_bus = 0.0, nothing
    for (i, bus) in ref[:bus]
        is_real_bus(bus) || continue
        4 in bus["terminals"] || continue
        nev = sqrt(JuMP.value(vr[4, i])^2 + JuMP.value(vi[4, i])^2)
        nev > worst && ((worst, worst_bus) = (nev, i))
    end
    return worst, worst_bus
end

# -----------------------------------------------------------------------
# NETWORK LOSSES.
#
# Series loss per branch from the SERIES current cs (not the terminal current
# cr, which includes shunt charging):
#     Re(S_loss) = csr' r csr + csi' r csi
# The reactance terms cancel exactly because x is symmetric, so only r
# appears -- this is pure I^2R.
#
# The sum runs over all FOUR conductors because the model is not Kron-reduced,
# so neutral-conductor loss is included. That loss is invisible in a 3-wire
# equivalent and is created entirely by unbalance.
#
# neutral_kw is the DIAGONAL term r[4,4]*|cs4|^2 only. With mutual coupling,
# attributing loss to one conductor is not unique -- the off-diagonal terms
# belong to no single conductor -- so it is an indicative share, and it does
# NOT sum with the phase contributions to equal series_kw.
#
# shunt_kw is real power dissipated in the MEN bonds, i.e. current returning
# through earth rather than through the neutral conductor.
# -----------------------------------------------------------------------
function network_losses(sbase_kva)
    p_series, p_neutral = 0.0, 0.0
    for (l, branch) in ref[:branch]
        r = branch["br_r"]
        n = size(r, 1)
        a = [JuMP.value(csr[c, l]) for c in 1:n]
        b = [JuMP.value(csi[c, l]) for c in 1:n]
        p_series += a' * r * a + b' * r * b
        n >= 4 && (p_neutral += r[4, 4] * (a[4]^2 + b[4]^2))
    end

    p_shunt = 0.0
    if haskey(ref, :shunt)
        for (_, sh) in ref[:shunt]
            g = sh["gs"]
            all(iszero, g) && continue
            i, conns = sh["shunt_bus"], sh["connections"]
            a = [JuMP.value(vr[t, i]) for t in conns]
            b = [JuMP.value(vi[t, i]) for t in conns]
            p_shunt += a' * g * a + b' * g * b
        end
    end

    load_pu = sum(sum(l["pd"]) for (_, l) in ref[:load])
    total   = p_series + p_shunt
    return (series_kw  = pu_to_kw(p_series,  sbase_kva),
            neutral_kw = pu_to_kw(p_neutral, sbase_kva),
            shunt_kw   = pu_to_kw(p_shunt,   sbase_kva),
            total_kw   = pu_to_kw(total,     sbase_kva),
            pct_load   = load_pu > 0 ? 100 * total / load_pu : NaN)
end

# -----------------------------------------------------------------------
# Branch currents, in amps.
#
# The four-conductor sum on a real line is NOT zero: current can leave through
# a MEN bond and return through earth, so the residual IS the earth-return
# current on that section. On the virtual source branch the residual is the
# whole return current re-emerging from earth at the substation bond.
# -----------------------------------------------------------------------
function branch_currents(sbase_kva)
    rows = NamedTuple[]
    for (l, br) in ref[:branch]
        startswith(get(br, "name", ""), "_virtual") && continue
        fb = br["f_bus"]
        vb = get(ref[:bus][fb], "vbase", 0.4 / sqrt(3))
        ibase = sbase_kva / vb
        n = length(br["f_connections"])
        I = [JuMP.value(csr[c, l]) + im * JuMP.value(csi[c, l]) for c in 1:n]
        push!(rows, (line = get(br, "name", string(l)),
                     f_bus = get(ref[:bus][fb], "name", string(fb)),
                     t_bus = get(ref[:bus][br["t_bus"]], "name", string(br["t_bus"])),
                     ia = abs(I[1]) * ibase, ib = abs(I[2]) * ibase, ic = abs(I[3]) * ibase,
                     in_ = n >= 4 ? abs(I[4]) * ibase : NaN,
                     iearth = abs(sum(I)) * ibase,
                     rating = (br["c_rating_a"][1] == Inf ? NaN : br["c_rating_a"][1] * ibase)))
    end
    sort!(rows, by = r -> r.line)
    return rows
end

function current_summary(sbase_kva)
    rows = branch_currents(sbase_kva)
    isempty(rows) && return (in_head=NaN, in_max=NaN, iph_max=NaN)
    return (in_head = rows[1].in_,
            in_max  = maximum(r.in_ for r in rows),
            iph_max = maximum(max(r.ia, r.ib, r.ic) for r in rows))
end

function statcom_dispatch(gen_ids, sbase_kva, rating_kva_total)
    isempty(gen_ids) && return nothing
    net_p, net_q = zeros(3), zeros(3)
    worst_util = 0.0
    leg = rating_kva_total / 3
    for gid in gen_ids
        p = [pu_to_kw(JuMP.value(pg[c, gid]), sbase_kva) for c in 1:3]
        q = [pu_to_kvar(JuMP.value(qg[c, gid]), sbase_kva) for c in 1:3]
        for c in 1:3
            net_p[c] += p[c]; net_q[c] += q[c]
            worst_util = max(worst_util, sqrt(p[c]^2 + q[c]^2) / leg)
        end
    end
    return (net_p=net_p, net_q=net_q, worst_util=worst_util, n=length(gen_ids))
end

function bus_profile(sbase_kva)
    rows = NamedTuple[]
    for (i, bus) in ref[:bus]
        is_virtual(bus) && continue
        t = bus["terminals"]
        vb = get(bus, "vbase", 0.4 / sqrt(3))
        vr_n = 4 in t ? JuMP.value(vr[4, i]) : 0.0
        vi_n = 4 in t ? JuMP.value(vi[4, i]) : 0.0
        nev = sqrt(vr_n^2 + vi_n^2)
        for p in 1:3
            p in t || continue
            push!(rows, (bus_id = i, bus_name = get(bus, "name", string(i)), phase = p,
                         vm_ground_pu = sqrt(JuMP.value(vr[p,i])^2 + JuMP.value(vi[p,i])^2),
                         vm_neutral_pu = sqrt((JuMP.value(vr[p,i])-vr_n)^2 + (JuMP.value(vi[p,i])-vi_n)^2),
                         nev_pu = nev, nev_v = nev * vb * 1000))
        end
    end
    return rows
end

# -----------------------------------------------------------------------
# One metric schema shared by every block, so a row can never be the wrong
# width. metric_row returns exactly length(METRIC_HEADER) entries whether the
# case solved or not.
# -----------------------------------------------------------------------
const METRIC_HEADER = [
    "status", "vuf_pct", "breaches_vuf",
    "vmin_pu", "vmin_bus", "vmin_phase", "vmax_pu", "vmax_bus", "vmax_phase",
    "vmin_ground_pu", "breaches_vmin",
    "nev_pu", "nev_bus", "v2_pu",
    "loss_total_kw", "loss_series_kw", "loss_neutral_kw", "loss_shunt_kw", "loss_pct_of_load",
    "i_neutral_head_a", "i_neutral_max_a", "i_phase_max_a",
    "p_ph1_kw", "p_ph2_kw", "p_ph3_kw", "p_sum_kw",
    "q_ph1_kvar", "q_ph2_kvar", "q_ph3_kvar", "util"]

function metric_row(status, gen_ids, sbase_kva, rating_kva)
    if !solved(status)
        return Any[string(status); fill(nothing, length(METRIC_HEADER) - 1)]
    end
    vuf, vuf_bus = worst_case_vuf()
    vm   = worst_case_vmag(ref_neutral=true)
    vmg  = worst_case_vmag(ref_neutral=false)
    nev, nev_bus = worst_case_nev()
    v2, _ = worst_neg_seq()
    ls   = network_losses(sbase_kva)
    cs   = current_summary(sbase_kva)
    d    = statcom_dispatch(gen_ids, sbase_kva, rating_kva)
    bn(i) = isnothing(i) ? nothing : get(ref[:bus][i], "name", string(i))

    return Any[
        "SOLVED", 100 * vuf, vuf > VUF_LIMIT_REF,
        vm.min_v, bn(vm.min_bus), vm.min_phase, vm.max_v, bn(vm.max_bus), vm.max_phase,
        vmg.min_v, vm.min_v < VMIN_PU,
        nev, bn(nev_bus), v2,
        ls.total_kw, ls.series_kw, ls.neutral_kw, ls.shunt_kw, ls.pct_load,
        cs.in_head, cs.in_max, cs.iph_max,
        isnothing(d) ? nothing : d.net_p[1], isnothing(d) ? nothing : d.net_p[2],
        isnothing(d) ? nothing : d.net_p[3], isnothing(d) ? nothing : sum(d.net_p),
        isnothing(d) ? nothing : d.net_q[1], isnothing(d) ? nothing : d.net_q[2],
        isnothing(d) ? nothing : d.net_q[3], isnothing(d) ? nothing : d.worst_util]
end

function print_metrics(tag, row)
    if row[1] != "SOLVED"
        @printf("  %-34s %s\n", tag, row[1]); return
    end
    @printf("  %-34s VUF=%6.3f%%%s |V|n=%.4f-%.4f%s NEV=%.4f In=%6.1fA loss=%6.3fkW (%5.2f%%)%s\n",
            tag, row[2], row[3] ? " BREACH" : "       ",
            row[4], row[7], row[11] ? " LOW" : "    ",
            row[12], row[21], row[15], row[19],
            isnothing(row[30]) ? "" : @sprintf(" util=%.3f", row[30]))
end

# -----------------------------------------------------------------------
# Minimal CSV writer -- no CSV.jl/DataFrames.jl dependency.
# -----------------------------------------------------------------------
# Fields are quoted when they contain a comma, quote or newline, per RFC 4180.
# Without this, any label containing a comma (e.g. a scenario name) silently
# writes an extra column and every reader sees a ragged row.
function csv_escape(x)
    x === nothing && return ""
    x isa AbstractFloat && isnan(x) && return ""
    s = string(x)
    if occursin(',', s) || occursin('"', s) || occursin('\n', s)
        return '"' * replace(s, '"' => "\"\"") * '"'
    end
    return s
end

function write_csv(path, header::Vector{String}, rows::Vector{Vector{Any}})
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(map(csv_escape, header), ","))
        for r in rows
            println(io, join(map(csv_escape, r), ","))
        end
    end
    println("  wrote $path  ($(length(rows)) rows)")
end

# =======================================================================
# ONE CASE = build network + place device + solve + report
# =======================================================================
function run_case(; alloc, mode, objective, rating_kva=STATCOM_KVA,
                  statcom_bus=STATCOM_BUS, bounds_active=false, vuf_cap=Inf)
    dm = load_base_network(data_path; enforce_bounds=bounds_active, men_ohm=MEN_OHM,
                           neutral_ref=NEUTRAL_REF,
                           lots_at_bus=(MEN_PER_LOT ? LOTS_AT_BUS : nothing))
    set_phase_allocation!(dm, SBASE_KVA; bus_name=LOAD_BUS_NAME, lots_per_phase=alloc)

    gids = mode == :none ? Int[] :
        add_statcom!(dm, SBASE_KVA; bus_name=statcom_bus,
                     rating_kva_total=rating_kva, mode=mode)
    global STATCOM_GEN_IDS = gids

    st = try
        build_and_solve(dm; objective_mode=objective, neutral_ref=NEUTRAL_REF,
                        bounds_active=bounds_active, vuf_cap=vuf_cap)
    catch e
        println("      ERROR: $e"); :ERROR
    end
    return st, gids
end

# =======================================================================
# RUN
# =======================================================================

_probe = load_base_network(data_path; men_ohm=MEN_OHM,
                           lots_at_bus=(MEN_PER_LOT ? LOTS_AT_BUS : nothing))
const SBASE_KVA = get_sbase_kva(_probe)
set_phase_allocation!(_probe, SBASE_KVA; bus_name=LOAD_BUS_NAME, lots_per_phase=ALLOCATIONS[1])
probe_kw, probe_kvar = total_load_kw_kvar(_probe, SBASE_KVA)

println("=" ^ 92)
println(" SMALL NETWORK SWEEP 2 -- 4-bus, 4-wire, MEN")
println(" sbase = $SBASE_KVA kVA   |   21 lots @ $ADMD_KVA kVA ADMD   |   total = $(round(probe_kw, digits=1)) kW / $(round(probe_kvar, digits=1)) kVAr")
println(" STATCOM = $N_STATCOMS x $STATCOM_KVA kVA ($(round(STATCOM_KVA/3, digits=1)) kVA/leg) at $STATCOM_BUS")
println(" MEN = $MEN_OHM Ω per lot" * (MEN_PER_LOT ? " (parallelled per bus)" : " (one bond per bus)"))
println(" buses: ", join(sort(collect(keys(_probe["bus_lookup"]))), ", "))
println("=" ^ 92)

summary_rows  = Vector{Any}[]
profile_rows  = Vector{Any}[]
current_rows  = Vector{Any}[]
rating_rows   = Vector{Any}[]
siting_rows   = Vector{Any}[]

# -----------------------------------------------------------------------
# BLOCK 1 -- allocation sweep, limits RELAXED and drawn as references.
#
# Relaxed is the honest setting for a study whose point is that the unaided
# network breaches: with the limits enforced, the solver cannot return a
# breaching operating point at all. It also removes a handicap from the
# comparator, not from the proposed device -- expect the gap to narrow.
# -----------------------------------------------------------------------
println("\n" * "=" ^ 92)
println(" BLOCK 1 -- phase allocation sweep   [limits RELAXED, 2% and 0.90 pu drawn as references]")
println("=" ^ 92)

for alloc in ALLOCATIONS
    lbl = alloc_label(alloc)
    println("\n── endbus allocation $lbl lots   (phase spread $(round(alloc_spread_kw(alloc), digits=1)) kW) ──")
    for (code, label, mode, obj) in vcat(SCENARIOS, SCENARIOS_COST)
        st, gids = run_case(alloc=alloc, mode=mode, objective=obj)
        row = metric_row(st, gids, SBASE_KVA, STATCOM_KVA)
        push!(summary_rows, Any[lbl, alloc_spread_kw(alloc), code, label, obj, row...])
        print_metrics("$code $label [$obj]", row)

        if solved(st)
            for r in bus_profile(SBASE_KVA)
                push!(profile_rows, Any[lbl, code, obj, r.bus_id, r.bus_name, r.phase,
                                        r.vm_ground_pu, r.vm_neutral_pu, r.nev_pu, r.nev_v])
            end
            for r in branch_currents(SBASE_KVA)
                push!(current_rows, Any[lbl, code, obj, r.line, r.f_bus, r.t_bus,
                                        r.ia, r.ib, r.ic, r.in_, r.iearth, r.rating])
            end
        end
    end
end

# -----------------------------------------------------------------------
# BLOCK 2 -- converter rating sweep at the worst allocation.
#
# Answers the first question a reviewer asks about a device that only
# partly fixes something: how big was it? The rating at which utilisation
# falls below 1.0 is the knee between "device too small" and "device
# cannot" -- a device that stays saturated at EVERY rating is the latter.
# -----------------------------------------------------------------------
println("\n" * "=" ^ 92)
println(" BLOCK 2 -- converter rating sweep at $(alloc_label(WORST_ALLOC)) lots   [limits RELAXED]")
println("=" ^ 92)

let st, gids
    st, gids = run_case(alloc=WORST_ALLOC, mode=:none, objective="cost")
    row = metric_row(st, gids, SBASE_KVA, STATCOM_KVA)
    push!(rating_rows, Any[0.0, "A", "no STATCOM", row...])
    print_metrics("   -- kVA  A no STATCOM", row)
end

for rating in RATING_KVA_LEVELS
    for (code, label, mode, _) in SCENARIOS
        mode == :none && continue
        st, gids = run_case(alloc=WORST_ALLOC, mode=mode, objective="VUF", rating_kva=rating)
        row = metric_row(st, gids, SBASE_KVA, rating)
        push!(rating_rows, Any[rating, code, label, row...])
        print_metrics(@sprintf("%6.1f kVA  %-2s %s", rating, code, label), row)
    end
end

# -----------------------------------------------------------------------
# BLOCK 3 -- the same allocation sweep with the statutory limits ENFORCED.
#
# This is the operational view: what a network bound by 0.90/1.10 pu and a
# 2% VUF cap can actually dispatch. Infeasible rows are a legitimate result
# here -- they mean no compliant dispatch exists for that case -- and are
# recorded rather than skipped.
# -----------------------------------------------------------------------
println("\n" * "=" ^ 92)
println(" BLOCK 3 -- same sweep, limits ENFORCED   [0.90/1.10 pu phase-to-neutral, VUF <= $(VUF_LIMIT_REF)]")
println(" Infeasible here means: no compliant dispatch exists for that case.")
println("=" ^ 92)

for alloc in ALLOCATIONS
    lbl = alloc_label(alloc)
    println("\n── endbus allocation $lbl lots ──")
    for (code, label, mode, obj) in SCENARIOS
        st, gids = run_case(alloc=alloc, mode=mode, objective=obj,
                            bounds_active=true, vuf_cap=VUF_LIMIT_REF)
        row = metric_row(st, gids, SBASE_KVA, STATCOM_KVA)
        push!(summary_rows, Any[lbl, alloc_spread_kw(alloc), code, label, obj * "-enforced", row...])
        print_metrics("$code $label [enforced]", row)
    end
end

# -----------------------------------------------------------------------
# BLOCK 4 -- siting.
#
# Two candidate buses rather than an optimisation, per the suggestion to
# SAMPLE locations. The device in every other block sits at endbus, which is
# best-case siting -- colocated with the unbalance it is correcting. This
# block prices that assumption before the same question is asked of the
# 908-bus feeder, where placement is not obvious.
# -----------------------------------------------------------------------
println("\n" * "=" ^ 92)
println(" BLOCK 4 -- siting at $(alloc_label(WORST_ALLOC)) lots, $STATCOM_KVA kVA   [limits RELAXED]")
println("=" ^ 92)

for bus in ["mid2", "endbus"]
    for (code, label, mode, obj) in SCENARIOS
        mode == :none && continue
        st, gids = run_case(alloc=WORST_ALLOC, mode=mode, objective=obj, statcom_bus=bus)
        row = metric_row(st, gids, SBASE_KVA, STATCOM_KVA)
        push!(siting_rows, Any[bus, code, label, row...])
        print_metrics("$bus  $code $label", row)
    end
end

# =======================================================================
# OUTPUT -- CSVs only. Figures are generated separately from these files, so
# a change to a figure never requires a re-solve.
# =======================================================================
println()

write_csv(joinpath(outdir, "sweep2_summary.csv"),
    ["allocation", "spread_kw", "scenario", "scenario_label", "objective", METRIC_HEADER...],
    summary_rows)

write_csv(joinpath(outdir, "sweep2_rating.csv"),
    ["rating_kva", "scenario", "scenario_label", METRIC_HEADER...],
    rating_rows)

write_csv(joinpath(outdir, "sweep2_siting.csv"),
    ["statcom_bus", "scenario", "scenario_label", METRIC_HEADER...],
    siting_rows)

write_csv(joinpath(outdir, "sweep2_profiles.csv"),
    ["allocation", "scenario", "objective", "bus_id", "bus_name", "phase",
     "vm_ground_pu", "vm_neutral_pu", "nev_pu", "nev_v"],
    profile_rows)

write_csv(joinpath(outdir, "sweep2_currents.csv"),
    ["allocation", "scenario", "objective", "line", "f_bus", "t_bus",
     "ia_a", "ib_a", "ic_a", "in_a", "iearth_a", "rating_a"],
    current_rows)

println("\nDone. Figures: run figures2.jl against $outdir")
println()
println("Checks worth making before plotting:")
println("  1. Balanced allocation 3/3/3, scenario A: neutral current, NEV and VUF all ~ 0.")
println("     If not, the MEN bonds or the network are wrong and nothing below means anything.")
println("  2. p_sum_kw ~ 0 for C and D. That is the inter-phase exchange signature: real power")
println("     is MOVED between phases, never imported. Check the PER-PHASE values too -- a sum")
println("     of zero is equally consistent with all three being individually zero.")
println("  3. p_ph*_kw all ~ 0 for B and Bc. If not, the Q-only path is leaking real power.")
println("  4. loss_total_kw at 3/3/3 is the FLOOR. No device can beat it -- it is what it costs")
println("     to deliver 89.8 kW through this feeder. Measure every scenario against that,")
println("     not against zero. The excess above it is the cost of unbalance.")
println("  5. i_neutral_max_a against i_phase_max_a. If the neutral carries more than the worst")
println("     phase, that is a thermal finding a 3-wire model cannot produce, and phase")
println("     overcurrent protection would never see it.")
println("  6. vmin_pu (neutral-referenced) should be BELOW vmin_ground_pu by roughly nev_pu.")
println("     If they are equal, the neutral referencing did not take effect.")
println("  7. Block 2: the rating at which util drops below 1.0 is the knee. A scenario pinned")
println("     at 1.000 across every rating is limited by CAPABILITY, not by size.")
println("  8. Block 3 infeasibility is a result, not a bug -- it means no compliant dispatch")
println("     exists. Read the status column in the CSV, because figures2.jl filters to SOLVED")
println("     and a dropped point looks identical to one that was never run.")