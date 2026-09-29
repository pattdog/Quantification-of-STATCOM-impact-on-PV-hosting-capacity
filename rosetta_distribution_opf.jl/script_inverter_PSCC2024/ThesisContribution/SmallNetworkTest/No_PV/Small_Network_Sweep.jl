#=
==============================================================================
TOY UNBALANCE SWEEP: 3-bus 4-wire minimal case, with/without STATCOM (v1)
==============================================================================
Fork of pv_statcom_sweep.jl. Purpose: demonstrate cause and effect on a small
network before presenting the 908-bus feeder.

KEPT IDENTICAL from the PV sweep, on purpose -- the toy case must exercise the
same code path as the feeder, or it validates nothing:
    load_base_network, add_men!, add_statcoms!, get_sbase_kva,
    total_load_kw_kvar, worst_case_vmag, aggregate_statcom_dispatch
    build_and_solve  (one change only: `objective` is now an argument
                      instead of a hard-coded "cost" -- see note there)

CHANGED:
  1. MEN IS ACTUALLY ON. In pv_statcom_sweep.jl, load_base_network defaults
     men_ohm=nothing and the sweep never passes it, so add_men! is dead code
     and MEN is absent from every PV result. Here men_ohm=10.0 is passed
     explicitly. If that was unintentional in the PV sweep, it is a bigger
     finding than anything in this file.

  2. STRESS AXIS is load unbalance, not PV penetration. Total demand is held
     at 18 kW while the phase split goes from 6/6/6 to 16/1.5/0.5. Unbalance
     severity is the axis on which Australian and ENWL LV networks actually
     differ, so this doubles as a start on the 1:38 validation ask.

  3. OBJECTIVE IS SWEPT. Scenarios B and C run under both "cost" and "VUF".
     Under "cost" the STATCOM only acts when a constraint binds, and on a
     network this lightly loaded nothing binds -- B and C come out identical
     and there is no figure. Under "VUF" the device is asked to minimise
     unbalance and B/C separate continuously. Both are reported because the
     contrast is itself a result.
     Scenario A is forced to "cost": with no controllable device the objective
     cannot change the answer, and it avoids passing an empty
     objective_target_gens to the VUF branch.

  4. VUF REFERENCE: PROVABLY IRRELEVANT, so only one VUF column is kept.
     Subtracting a common neutral voltage from all three phases adds a purely
     ZERO-sequence term; since 1 + a + a^2 = 0, both V1 and V2 are unchanged
     and |V2|/|V1| is invariant. Confirmed empirically: vuf_neutral and
     vuf_ground agreed to 3 decimals in all 30 rows of the first run. The
     earlier claim that this was an omission was wrong.
     (Superseded note) VUF IS REFERENCED TO NEUTRAL. worst_case_vuf in the PV sweep computes
     phase-to-GROUND. With MEN on and neutral displacement being the thing
     under study, phase-to-ground and phase-to-neutral are different
     quantities. Both are computed here so the change is visible rather than
     silent.

  5. VOLTAGE BOUNDS REFERENCED TO NEUTRAL, not ground. constraints_unified.jl
     enforces vr[t]^2+vi[t]^2 >= vmin[t]^2, which is phase-to-GROUND. On a MEN
     network that differs from what a customer sees by the neutral
     displacement -- 0.043 pu at the worst split here. The built-in phase
     bounds are switched off (vmin=0, vmax=Inf, which that file's own filters
     skip) and add_neutral_referenced_bounds! imposes the correct ones.
     constraints_unified.jl itself is untouched, since it is shared with the
     908-bus runs. Set NEUTRAL_REF=false to reproduce the old behaviour.

  6. RATING SWEEP added as block 2. worst_leg_utilisation was pinned at
     0.9999 for every unbalance index above ~1.2, so those results measure
     converter size, not the formulation.

  7. NO TERMINAL OUTPUT AS DELIVERABLE. Two CSVs are written; figures are
     generated separately from them. Console printing is progress only.

Run with:
    julia --project=. toy_unbalance_sweep.jl
==============================================================================
=#

using Logging
Logging.disable_logging(Logging.Warn)

using Pkg

# --- Path anchoring (depth-independent) ---------------------------------
# Locations are found by SEARCHING UPWARD for marker contents, not by
# counting ".." segments. Moving this script to a different depth (e.g. into
# a No_PV/ subfolder) then requires no edits.
#
#   PSCC_DIR    -- nearest ancestor containing core/variables.jl
#   RPMD_ROOT   -- nearest ancestor containing both data/ and src/
#   THESIS_ROOT -- OUTERMOST ancestor with a Project.toml. Both
#                  rosetta_distribution_opf.jl and thesis_STATCOM have one;
#                  the outer is the environment that was originally
#                  activated, so the search deliberately keeps going up.
const SCRIPT_DIR = @__DIR__

function _find_up(start::AbstractString, pred::Function; outermost::Bool=false, maxdepth::Int=10)
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

const PSCC_DIR = _find_up(SCRIPT_DIR, d -> isfile(joinpath(d, "core", "variables.jl")))
const RPMD_ROOT = _find_up(SCRIPT_DIR, d -> isdir(joinpath(d, "data")) && isdir(joinpath(d, "src")))
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

ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "sb" => "yes",)

# <-- data_path: set to wherever Master_3bus_4w.txt actually lives.
#     Two sensible options, pick one:
data_path = joinpath(RPMD_ROOT, "data", "3_Bus_Small.dss")
# data_path = joinpath(SCRIPT_DIR, "3_Bus_Small.dss")   # if kept beside the script

outdir    = joinpath(THESIS_ROOT, "results", "toy_unbalance")

isfile(data_path) || error("network file not found: $data_path")

# -----------------------------------------------------------------------
# kW/kVAr <-> pu conversions
# -----------------------------------------------------------------------
kw_to_pu(p_kw, sbase_kva)     = p_kw / sbase_kva
kvar_to_pu(q_kvar, sbase_kva) = q_kvar / sbase_kva
pu_to_kw(p_pu, sbase_kva)     = p_pu * sbase_kva
pu_to_kvar(q_pu, sbase_kva)   = q_pu * sbase_kva

# -----------------------------------------------------------------------
# MEN -- unchanged from pv_statcom_sweep.jl
# -----------------------------------------------------------------------
function add_men!(data_eng; r_load_ohm=10.0, verbose=true)
    haskey(data_eng, "shunt") || (data_eng["shunt"] = Dict{String,Any}())

    men_buses = String[]
    for (_, load) in data_eng["load"]
        b = load["bus"]
        4 in load["connections"] || continue
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
# Network loader -- unchanged
# -----------------------------------------------------------------------
function load_base_network(data_path; load_multiplier=1.0, enforce_bounds=true, sbase_kva=100, men_ohm=nothing, neutral_ref=true)
    data_eng  = PMD.parse_file(data_path, transformations=[PMD.transform_loops!])
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
        elseif neutral_ref
            # BOUNDS REFERENCE -- the important correction in this version.
            #
            # constraints_unified.jl enforces  vr[t]^2 + vi[t]^2 >= vmin[t]^2,
            # i.e. phase-to-GROUND. On a MEN network the quantity a customer's
            # appliance sees is phase-to-NEUTRAL, and the two differ by exactly
            # the neutral displacement. At the 16/1.5/0.5 split that gap was
            # 0.043 pu: the OPF believed phase A had 5.4% of headroom when the
            # customer had 1.1%.
            #
            # vmin=0 and vmax=Inf on the phase terminals are read as "skip" by
            # the two filters in constraints_unified.jl (nonzero_vmin_terminals
            # / nonInf_vmax_terminals), so the built-in phase bounds go away
            # and add_neutral_referenced_bounds! imposes the correct ones.
            #
            # Terminal 4 keeps its built-in vmax: neutral-to-GROUND is
            # genuinely what the NEV cap should measure.
            bus["vmin"] = [0.00, 0.00, 0.00, 0.00]
            bus["vmax"] = [Inf,  Inf,  Inf,  enforce_bounds ? 0.20 : 1.50]
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

get_sbase_kva(data_math) = data_math["settings"]["sbase"] * data_math["settings"]["power_scale_factor"] / 1000

function total_load_kw_kvar(data_math, sbase_kva)
    total_pd_pu = sum(sum(l["pd"]) for (i,l) in data_math["load"])
    total_qd_pu = sum(sum(l["qd"]) for (i,l) in data_math["load"])
    return pu_to_kw(total_pd_pu, sbase_kva), pu_to_kvar(total_qd_pu, sbase_kva)
end

# -----------------------------------------------------------------------
# UNBALANCE SETTER -- replaces add_pv! as the stress mechanism.
#
# Overwrites pd/qd on the three single-phase loads according to a phase
# split given in kW. PF held at 0.95 lagging throughout, so the only thing
# changing across the sweep is HOW the same total demand is distributed.
#
# Loads are matched to phases by connection, not by name or dict order --
# dict iteration order is not guaranteed and load ids are strings.
# -----------------------------------------------------------------------
function set_unbalance!(data_math, sbase_kva; split_kw, pf=0.95)
    tanphi = tan(acos(pf))
    assigned = Dict(1 => false, 2 => false, 3 => false)

    for (_, load) in data_math["load"]
        phase_conns = filter(c -> c != 4, load["connections"])
        length(phase_conns) == 1 || error("set_unbalance! expects single-phase loads; got connections $(load["connections"])")
        ph = phase_conns[1]
        haskey(assigned, ph) || error("load on unexpected phase $ph")
        assigned[ph] && error("more than one load on phase $ph -- set_unbalance! assumes exactly one per phase")

        p_kw = split_kw[ph]
        load["pd"] = [kw_to_pu(p_kw, sbase_kva)]
        load["qd"] = [kw_to_pu(p_kw * tanphi, sbase_kva)]
        assigned[ph] = true
    end

    all(values(assigned)) || error("not every phase received a load: $assigned")
    return nothing
end

# Load unbalance index used as the x-axis: spread relative to the mean.
# 0 for a perfectly balanced split, 1.0 when max-min equals the mean.
load_unbalance_index(split_kw) = (maximum(split_kw) - minimum(split_kw)) / mean(split_kw)

# -----------------------------------------------------------------------
# STATCOM PLACEMENT -- unchanged from pv_statcom_sweep.jl.
#
# Note carried over, not fixed here: this function never sets
# gen["connections"] and never trims inherited length-3 fields, unlike
# add_pv!. It relies entirely on the deepcopy from gen "1". That is latent
# rather than active while every unit is 3-phase, but it is the same class
# of bug add_pv!'s BUGFIX comment describes. Left untouched so the toy case
# exercises the identical code path as the feeder.
# -----------------------------------------------------------------------
function add_statcoms!(data_math, sbase_kva; n_units, rating_kvar_total, p_exchange::Bool)
    load_ids = sort(
        collect(keys(data_math["load"])),
        by = x -> tryparse(Int, x) === nothing ? 0 : parse(Int, x)
    )
    n_loads = length(load_ids)
    spacing = max(1, div(n_loads, n_units))

    q_total_pu = kvar_to_pu(rating_kvar_total, sbase_kva)
    q_leg_pu   = q_total_pu / 3

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

    return gen_ids
end

# -----------------------------------------------------------------------
# Neutral-referenced phase voltage bounds.
#
# Same nonconvex quadratic form as constraints_unified.jl, so Ipopt sees no
# new class of constraint -- only the reference point moves:
#     (vr[p] - vr[4])^2 + (vi[p] - vi[4])^2   in   [vmin^2, vmax^2]
#
# constraints_unified.jl is deliberately NOT modified: it is shared with the
# 908-bus runs, and changing it would silently invalidate results already
# produced. This function is additive and opt-in instead.
#
# Must run AFTER constraints_unified.jl is included (vr/vi must exist) and
# ONLY when load_base_network was called with neutral_ref=true, or both sets
# of bounds are active at once and the tighter one silently wins.
# -----------------------------------------------------------------------
function add_neutral_referenced_bounds!(; vmin=0.90, vmax=1.10)
    n = 0
    for (i, bus) in ref[:bus]
        is_real_bus(bus) || continue
        terms = bus["terminals"]
        4 in terms || continue
        for p in 1:3
            p in terms || continue
            dvr = vr[p,i] - vr[4,i]
            dvi = vi[p,i] - vi[4,i]
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
# ONLY CHANGE from pv_statcom_sweep.jl: `objective` is passed in rather
# than hard-coded to "cost". Everything else -- instantiate_mc_model with
# the no-op build, the raw includes, the global ref/model handoff -- is
# byte-identical, because that pipeline is what is already validated.
# -----------------------------------------------------------------------
function build_and_solve(data_math; objective_mode::String="cost", neutral_ref::Bool=true,
                         vmin=0.90, vmax=1.10, bounds_active::Bool=true,
                         vuf_cap=0.02)
    pm = PMD.instantiate_mc_model(
        data_math,
        PMD.IVRENPowerModel,
        (pm) -> nothing;
        jump_model = JuMP.Model(ipopt_solver),
        ref_extensions = Function[]
    )

    global ref   = pm.ref[:it][:pmd][:nw][0]
    global model = pm.model

    # Read by the guarded VUF constraint in constraints_unified.jl. Must be
    # set on EVERY call, not once: it is a global that persists across
    # scenarios in the same session, so a stale value silently carries over.
    # Inf disables the cap entirely; 0.02 reproduces the original behaviour.
    global VUF_CAP = vuf_cap
    # Absolute paths built from CORE_DIR, so these survive the script being
    # moved to a different depth. (A literal "../../core/..." is file-relative
    # and therefore depth-dependent.)
    include(joinpath(CORE_DIR, "variables.jl"))
    include(joinpath(CORE_DIR, "constraints_unified.jl"))

    # bounds_active=false relaxes the phase voltage limits entirely. For a
    # scenario with no controllable device this changes nothing about the
    # answer -- the power flow is fully determined -- it only stops the
    # limits acting as a feasibility filter that forbids the very violation
    # being demonstrated. Applied uniformly across scenarios so no case is
    # ORDERED to comply while another is not; limits are drawn as reference
    # lines on the figures instead.
    neutral_ref && bounds_active && add_neutral_referenced_bounds!(vmin=vmin, vmax=vmax)

    global objective = objective_mode
    global objective_target_gens = [i for (i, g) in ref[:gen] if get(g, "type", "") == "STATCOM"]
    include(joinpath(CORE_DIR, "objectives_FIXED.jl"))

    JuMP.optimize!(model)
    return JuMP.termination_status(model)
end

# -----------------------------------------------------------------------
# Reporting helpers
# -----------------------------------------------------------------------

# Both references computed. ref_neutral=true is the physically meaningful
# one on a MEN network: phase-to-ground and phase-to-neutral differ by
# exactly the neutral displacement, which is the quantity under study.
function worst_case_vuf(; bus_ids=nothing, ref_neutral::Bool=true)
    alpha = exp(im*2/3*pi)
    T = 1/3 * [1 1 1 ; 1 alpha alpha^2 ; 1 alpha^2 alpha]
    worst_vuf, worst_bus = 0.0, nothing
    for (i, bus) in ref[:bus]
        !isnothing(bus_ids) && !(i in bus_ids) && continue
        startswith(get(bus, "name", ""), "_virtual") && continue
        terms = bus["terminals"]
        !(1 in terms && 2 in terms && 3 in terms) && continue
        ref_neutral && !(4 in terms) && continue

        vr_n = ref_neutral ? JuMP.value(vr[4,i]) : 0.0
        vi_n = ref_neutral ? JuMP.value(vi[4,i]) : 0.0

        vph = [(JuMP.value(vr[p,i]) - vr_n) + im*(JuMP.value(vi[p,i]) - vi_n) for p in 1:3]
        v012 = T * vph
        vpos, vneg = abs(v012[2]), abs(v012[3])
        if vpos > 1e-6 && vneg/vpos > worst_vuf
            worst_vuf, worst_bus = vneg/vpos, i
        end
    end
    return worst_vuf, worst_bus
end

# PATCH 3: skip PMD's internal virtual source bus. It is not bus_type 3, so
# the original filter missed it, and it pinned every reported Vmax at ~1.0000
# regardless of scenario -- the column was reading the virtual bus, not a
# real one.
is_real_bus(bus) = get(bus, "bus_type", 1) != 3 && !startswith(get(bus, "name", ""), "_virtual")

# PATCH 1 (reporting half): both references returned.
#   ground  = vr[p]^2 + vi[p]^2                     <- what constraints_unified.jl bounds
#   neutral = (vr[p]-vr[4])^2 + (vi[p]-vi[4])^2     <- what a customer's appliance sees
# On a MEN network these differ by the neutral displacement. At the 16/1.5/0.5
# split that gap was 0.0426 pu: the OPF believed phase A had 5.4% of headroom
# when the real figure was 1.1%.
function worst_case_vmag(; bus_ids=nothing, ref_neutral::Bool=false)
    min_v, min_bus, min_phase = Inf, nothing, nothing
    max_v, max_bus, max_phase = -Inf, nothing, nothing
    for (i, bus) in ref[:bus]
        !isnothing(bus_ids) && !(i in bus_ids) && continue
        is_real_bus(bus) || continue
        terms = bus["terminals"]
        ref_neutral && !(4 in terms) && continue
        vr_n = ref_neutral ? JuMP.value(vr[4,i]) : 0.0
        vi_n = ref_neutral ? JuMP.value(vi[4,i]) : 0.0
        for p in 1:3
            p in terms || continue
            vm = sqrt((JuMP.value(vr[p,i]) - vr_n)^2 + (JuMP.value(vi[p,i]) - vi_n)^2)
            if vm < min_v; min_v, min_bus, min_phase = vm, i, p; end
            if vm > max_v; max_v, max_bus, max_phase = vm, i, p; end
        end
    end
    return (min_v=min_v, min_bus=min_bus, min_phase=min_phase, max_v=max_v, max_bus=max_bus, max_phase=max_phase)
end

# -----------------------------------------------------------------------
# Worst |V2| in pu, against the HARD-CODED 0.02 limit in
# constraints_unified.jl:
#
#     JuMP.@constraint(model, v_neg_seq_real^2 + v_neg_seq_imag^2 <= 0.02^2)
#
# That constraint is unconditional -- every bus, every run, no data field
# and no guard. It is invisible from this script, so it is reported here.
# If v2_max_pu approaches 0.02, results are being shaped by that constant
# rather than by the network or the device, and pushing unbalance or
# LOAD_SCALE further will make the NO-STATCOM case infeasible. That would
# look like "the network needs a STATCOM" when it is an artefact.
# -----------------------------------------------------------------------
const V2_HARD_LIMIT = 0.02

function worst_neg_seq()
    alpha = exp(im*2/3*pi)
    T = 1/3 * [1 1 1 ; 1 alpha alpha^2 ; 1 alpha^2 alpha]
    worst, worst_bus = 0.0, nothing
    for (i, bus) in ref[:bus]
        # Only the virtual bus is excluded, NOT the source bus: the 0.02
        # constraint in constraints_unified.jl applies at every bus in
        # ref[:bus], so the reported maximum should cover the same set.
        startswith(get(bus, "name", ""), "_virtual") && continue
        terms = bus["terminals"]
        !(1 in terms && 2 in terms && 3 in terms) && continue
        vph = [JuMP.value(vr[p,i]) + im*JuMP.value(vi[p,i]) for p in 1:3]
        v2  = abs((T * vph)[3])
        v2 > worst && ((worst, worst_bus) = (v2, i))
    end
    return worst, worst_bus
end

# -----------------------------------------------------------------------
# NETWORK LOSSES.
#
# Series loss per branch, from the SERIES current cs (not the terminal
# current cr, which includes shunt charging):
#     S_loss = (Z cs) . conj(cs),  Z = r + jx
#     Re(S_loss) = csr' r csr + csi' r csi
# The reactance terms cancel exactly because x is symmetric, so only r
# appears -- worth knowing, since it means this is pure I^2R.
#
# Because the model is 4-wire and NOT Kron-reduced, the sum runs over all
# four conductors, so neutral-conductor loss is included. That loss is
# invisible in a 3-wire equivalent and is precisely what unbalanced loading
# creates, which makes it the sharpest performance metric here.
#
# p_neutral is the DIAGONAL contribution r[4,4]*|cs4|^2 only. With mutual
# coupling, attributing loss to one conductor is not unique -- the
# off-diagonal terms belong to no single conductor -- so treat it as an
# indicative share, not an exact decomposition. It does NOT sum with the
# phase contributions to p_series.
#
# Shunt loss is the real power dissipated in the MEN bonds, gs*|v|^2, i.e.
# current returning through earth rather than through the neutral.
# -----------------------------------------------------------------------
function network_losses(sbase_kva)
    p_series  = 0.0
    p_neutral = 0.0
    for (l, branch) in ref[:branch]
        r = branch["br_r"]
        n = size(r, 1)
        csr_v = [JuMP.value(csr[c, l]) for c in 1:n]
        csi_v = [JuMP.value(csi[c, l]) for c in 1:n]
        p_series += csr_v' * r * csr_v + csi_v' * r * csi_v
        n >= 4 && (p_neutral += r[4,4] * (csr_v[4]^2 + csi_v[4]^2))
    end

    p_shunt = 0.0
    if haskey(ref, :shunt)
        for (_, sh) in ref[:shunt]
            g = sh["gs"]
            all(iszero, g) && continue
            i     = sh["shunt_bus"]
            conns = sh["connections"]
            v_r = [JuMP.value(vr[t, i]) for t in conns]
            v_i = [JuMP.value(vi[t, i]) for t in conns]
            p_shunt += v_r' * g * v_r + v_i' * g * v_i
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

# Worst neutral-to-earth voltage, in pu and volts.
function worst_case_nev()
    worst, worst_bus = 0.0, nothing
    for (i, bus) in ref[:bus]
        is_real_bus(bus) || continue
        4 in bus["terminals"] || continue
        nev = sqrt(JuMP.value(vr[4,i])^2 + JuMP.value(vi[4,i])^2)
        if nev > worst
            worst, worst_bus = nev, i
        end
    end
    return worst, worst_bus
end

# Unchanged from pv_statcom_sweep.jl
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

# Per-bus per-phase profile, for the voltage-profile figure.
function bus_profile(sbase_kva)
    rows = NamedTuple[]
    for (i, bus) in ref[:bus]
        # PATCH 3: virtual bus has no distance and would plot on top of sourcebus
        startswith(get(bus, "name", ""), "_virtual") && continue
        terms = bus["terminals"]
        vbase_kv = get(bus, "vbase", NaN)
        vr_n = 4 in terms ? JuMP.value(vr[4,i]) : 0.0
        vi_n = 4 in terms ? JuMP.value(vi[4,i]) : 0.0
        nev  = sqrt(vr_n^2 + vi_n^2)
        for p in 1:3
            p in terms || continue
            vm_g = sqrt(JuMP.value(vr[p,i])^2 + JuMP.value(vi[p,i])^2)
            vm_n = sqrt((JuMP.value(vr[p,i])-vr_n)^2 + (JuMP.value(vi[p,i])-vi_n)^2)
            push!(rows, (bus_id=i,
                         bus_name=get(bus, "name", string(i)),
                         phase=p,
                         vm_ground_pu=vm_g,
                         vm_neutral_pu=vm_n,
                         nev_pu=nev,
                         nev_v=nev * vbase_kv * 1000))
        end
    end
    return rows
end

# -----------------------------------------------------------------------
# Minimal CSV writer -- no new package dependency.
# -----------------------------------------------------------------------
function write_csv(path, header::Vector{String}, rows::Vector{Vector{Any}})
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(header, ","))
        for r in rows
            println(io, join(map(x -> x === nothing ? "" : (x isa AbstractFloat && isnan(x) ? "" : string(x)), r), ","))
        end
    end
    println("  wrote $path  ($(length(rows)) rows)")
end

# =======================================================================
# SWEEP CONFIGURATION
# =======================================================================

# Total demand held constant at 18 kW; only the split changes.
SPLITS = [
    Dict(1 =>  6.0, 2 => 6.0,  3 => 6.0),
    Dict(1 =>  8.0, 2 => 5.5,  3 => 4.5),
    Dict(1 => 10.0, 2 => 5.0,  3 => 3.0),
    Dict(1 => 12.0, 2 => 4.0,  3 => 2.0),
    Dict(1 => 14.0, 2 => 3.0,  3 => 1.0),
    Dict(1 => 16.0, 2 => 1.5,  3 => 0.5),
]

# Raise to ~2.0 to drive phase A below 0.90 pu and make the "cost" objective
# actually bind. At 1.0 nothing binds, which is the point being demonstrated.
LOAD_SCALE = 1.0

# Phase voltage bounds referenced to neutral rather than ground. See
# add_neutral_referenced_bounds!. Set false to reproduce the old behaviour.
# Statutory limits: ENFORCED as constraints, or RELAXED and drawn as
# reference lines on the figures. Relaxed is the honest setting for a study
# whose point is that the unaided network breaches them: with the cap on,
# a breach cannot be returned by the solver at all, only an infeasibility.
# Note this also unhandicaps the Q-only comparator, which currently sits
# pinned at |V|n = 0.9000 -- expect the gap to NARROW, not widen.
BOUNDS_ACTIVE = false          # phase voltage limits 0.90/1.10
VUF_CAP_PU    = Inf            # Inf disables; 0.02 restores the original cap
VUF_LIMIT_REF = 0.02           # drawn on figures, never enforced

NEUTRAL_REF   = true
VMIN_PU       = 0.90
VMAX_PU       = 1.10

MEN_OHM       = 10.0
N_STATCOMS    = 1
STATCOM_KVAR  = 18     # ~30% of the 18 kW total. Put this number in the
                        # figure caption: a device sized to the network is a
                        # result, one sized to 7x the network is an artefact.

# Rating sweep (block 2): the unbalance sweep showed worst_leg_utilisation
# pinned at 0.9999 from an unbalance index of ~1.2 onward, so every result
# past that point is a converter-size limit, not a formulation limit. This
# axis separates the two.
RATING_KVAR_LEVELS = [1.5, 3.0, 4.5, 6.0, 9.0, 12.0, 18.0, 27.0]
RATING_SPLIT       = SPLITS[end]        # most unbalanced case

# Block 3: load-scale stress. Holds the split and drives total demand up
# until the UNAIDED network exceeds the 2% VUF limit, then asks whether each
# device can bring it back. This is the result the solver was never told to
# produce, and it only exists with VUF_CAP_PU = Inf.
# Stops at 2.0 deliberately. At 2.4x the unaided network sits at |V|n = 0.734
# (169 V) with losses at 31.6% of load -- a numerical solution rather than a
# physical operating state, since constant-power loads at that voltage are a
# modelling fiction. The story is complete by 2.0: A and B both breach 2% VUF
# at 1.8 and 2.0 while C complies, which is the result.
LOAD_SCALE_LEVELS = [1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
STRESS_KVAR       = 18.0       # rating at which P-exchange left its S² circle

scenarios = [
    ("A: no STATCOM",           :none,      "cost"),
    ("B: STATCOM Q-only",       :qonly,     "cost"),
    ("C: STATCOM P-exchange",   :pexchange, "cost"),
    ("B: STATCOM Q-only",       :qonly,     "VUF"),
    ("C: STATCOM P-exchange",   :pexchange, "VUF"),
]

# =======================================================================
# RUN
# =======================================================================

_dm_probe = load_base_network(data_path; men_ohm=MEN_OHM)
SBASE_KVA = get_sbase_kva(_dm_probe)
set_unbalance!(_dm_probe, SBASE_KVA; split_kw=SPLITS[1])
base_pd_kw, base_qd_kvar = total_load_kw_kvar(_dm_probe, SBASE_KVA)

println("="^80)
println(" TOY UNBALANCE SWEEP -- 3-bus 4-wire")
println(" sbase = $SBASE_KVA kVA   |   total demand = $(round(base_pd_kw,digits=2)) kW / $(round(base_qd_kvar,digits=2)) kVAr")
println(" MEN = $MEN_OHM Ω   |   STATCOM = $N_STATCOMS x $STATCOM_KVAR kVAr   |   load scale = $LOAD_SCALE")
println(" Bus names present: ", sort([get(b,"name",string(i)) for (i,b) in _dm_probe["bus"]]))
println("="^80)

summary_rows = Vector{Any}[]
profile_rows = Vector{Any}[]

for split in SPLITS
    ubi = load_unbalance_index([split[1], split[2], split[3]])
    split_label = "$(split[1])/$(split[2])/$(split[3])"
    println("\n── split $split_label kW   (unbalance index $(round(ubi,digits=3))) ──")

    for (label, kind, obj) in scenarios
        scn = Dict(:none=>"A", :qonly=>"B", :pexchange=>"C")[kind]

        dm = load_base_network(data_path; load_multiplier=LOAD_SCALE,
                               enforce_bounds=true, men_ohm=MEN_OHM,
                               neutral_ref=NEUTRAL_REF)
        set_unbalance!(dm, SBASE_KVA; split_kw=split)
        LOAD_SCALE == 1.0 || (for (_,l) in dm["load"]; l["pd"] *= LOAD_SCALE; l["qd"] *= LOAD_SCALE; end)

        statcom_gen_ids = kind == :none ? Int[] :
            add_statcoms!(dm, SBASE_KVA; n_units=N_STATCOMS,
                          rating_kvar_total=STATCOM_KVAR,
                          p_exchange=(kind == :pexchange))
        global STATCOM_GEN_IDS = statcom_gen_ids

        status = try
            build_and_solve(dm; objective_mode=obj, neutral_ref=NEUTRAL_REF,
                            vmin=VMIN_PU, vmax=VMAX_PU,
                            bounds_active=BOUNDS_ACTIVE, vuf_cap=VUF_CAP_PU)
        catch e
            println("  $label [$obj]: ERROR -- $e")
            :ERROR
        end

        solved = status in [JuMP.LOCALLY_SOLVED, JuMP.OPTIMAL, JuMP.ALMOST_LOCALLY_SOLVED]

        if solved
            vuf_n, vuf_n_bus = worst_case_vuf(ref_neutral=true)
            vmag             = worst_case_vmag(ref_neutral=true)   # what a customer sees
            vmag_g           = worst_case_vmag(ref_neutral=false)  # what the old bounds measured
            nev, nev_bus     = worst_case_nev()
            v2, _            = worst_neg_seq()
            loss             = network_losses(SBASE_KVA)
            disp             = aggregate_statcom_dispatch(statcom_gen_ids, SBASE_KVA, STATCOM_KVAR)

            push!(summary_rows, Any[
                split_label, ubi, scn, obj, "SOLVED",
                vuf_n, vuf_n_bus,
                vmag.min_v, vmag.max_v, vmag.min_phase, vmag.max_phase,
                vmag_g.min_v, vmag_g.max_v,
                nev, nev_bus, v2,
                loss.total_kw, loss.series_kw, loss.neutral_kw, loss.shunt_kw, loss.pct_load,
                isnothing(disp) ? nothing : disp.net_pg[1],
                isnothing(disp) ? nothing : disp.net_pg[2],
                isnothing(disp) ? nothing : disp.net_pg[3],
                isnothing(disp) ? nothing : sum(disp.net_pg),
                isnothing(disp) ? nothing : disp.net_qg[1],
                isnothing(disp) ? nothing : disp.net_qg[2],
                isnothing(disp) ? nothing : disp.net_qg[3],
                isnothing(disp) ? nothing : disp.worst_util,
            ])

            for r in bus_profile(SBASE_KVA)
                push!(profile_rows, Any[split_label, ubi, scn, obj,
                                        r.bus_id, r.bus_name, r.phase,
                                        r.vm_ground_pu, r.vm_neutral_pu, r.nev_pu, r.nev_v])
            end

            @printf("  %-24s [%-4s] VUF=%6.3f%%  |V|n=%.4f-%.4f  |V|g=%.4f  NEV=%.4f  loss=%.3f kW (%.2f%%)%s\n",
                    label, obj, 100*vuf_n, vmag.min_v, vmag.max_v, vmag_g.min_v, nev,
                    loss.total_kw, loss.pct_load,
                    isnothing(disp) ? "" : @sprintf("  ΣP=%+.3f kW", sum(disp.net_pg)))
        else
            push!(summary_rows, Any[split_label, ubi, scn, obj, string(status),
                                    NaN, nothing, NaN, NaN, nothing, nothing, NaN, NaN, NaN, nothing, NaN,
                                    NaN, NaN, NaN, NaN, NaN,
                                    nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing])
            println("  $label [$obj]: $status")
        end
    end
end

# =======================================================================
# BLOCK 2 -- RATING SWEEP at the most unbalanced split.
#
# Answers the question a reviewer asks first about a device that only halves
# VUF: how big was it? Holds the load split fixed and sweeps converter size,
# so the knee between "rating-limited" and "formulation-limited" is visible.
# =======================================================================

rating_rows = Vector{Any}[]
rating_split_label = "$(RATING_SPLIT[1])/$(RATING_SPLIT[2])/$(RATING_SPLIT[3])"

println("\n" * "="^80)
println(" BLOCK 2 -- rating sweep at $rating_split_label kW split")
println("="^80)

# Baseline: no device, for the horizontal reference line on the figure.
let
    dm = load_base_network(data_path; load_multiplier=LOAD_SCALE, enforce_bounds=true,
                           men_ohm=MEN_OHM, neutral_ref=NEUTRAL_REF)
    set_unbalance!(dm, SBASE_KVA; split_kw=RATING_SPLIT)
    global STATCOM_GEN_IDS = Int[]
    st = build_and_solve(dm; objective_mode="cost", neutral_ref=NEUTRAL_REF,
                         vmin=VMIN_PU, vmax=VMAX_PU,
                         bounds_active=BOUNDS_ACTIVE, vuf_cap=VUF_CAP_PU)
    if st in [JuMP.LOCALLY_SOLVED, JuMP.OPTIMAL, JuMP.ALMOST_LOCALLY_SOLVED]
        vuf, _ = worst_case_vuf(ref_neutral=true)
        vm     = worst_case_vmag(ref_neutral=true)
        nev, _ = worst_case_nev()
        ls     = network_losses(SBASE_KVA)
        push!(rating_rows, Any[0.0, "A", "SOLVED", vuf, vm.min_v, vm.max_v, nev,
                               nothing, nothing, nothing, nothing,
                               nothing, nothing, nothing, nothing,
                               ls.total_kw, ls.neutral_kw, ls.pct_load])
        @printf("  %6s kVAr  A  VUF=%6.3f%%  |V|n=%.4f-%.4f  NEV=%.4f  loss=%.3f kW\n",
                "--", 100*vuf, vm.min_v, vm.max_v, nev, ls.total_kw)
    end
end

for rating in RATING_KVAR_LEVELS
    for (kind, scn) in [(:qonly, "B"), (:pexchange, "C")]
        dm = load_base_network(data_path; load_multiplier=LOAD_SCALE, enforce_bounds=true,
                               men_ohm=MEN_OHM, neutral_ref=NEUTRAL_REF)
        set_unbalance!(dm, SBASE_KVA; split_kw=RATING_SPLIT)

        gids = add_statcoms!(dm, SBASE_KVA; n_units=N_STATCOMS,
                             rating_kvar_total=rating, p_exchange=(kind == :pexchange))
        global STATCOM_GEN_IDS = gids

        st = try
            build_and_solve(dm; objective_mode="VUF", neutral_ref=NEUTRAL_REF,
                            vmin=VMIN_PU, vmax=VMAX_PU,
                            bounds_active=BOUNDS_ACTIVE, vuf_cap=VUF_CAP_PU)
        catch e
            println("  $rating kVAr $scn: ERROR -- $e"); :ERROR
        end

        if st in [JuMP.LOCALLY_SOLVED, JuMP.OPTIMAL, JuMP.ALMOST_LOCALLY_SOLVED]
            vuf, _ = worst_case_vuf(ref_neutral=true)
            vm     = worst_case_vmag(ref_neutral=true)
            nev, _ = worst_case_nev()
            d      = aggregate_statcom_dispatch(gids, SBASE_KVA, rating)
            ls     = network_losses(SBASE_KVA)
            push!(rating_rows, Any[rating, scn, "SOLVED", vuf, vm.min_v, vm.max_v, nev,
                                   d.net_pg[1], d.net_pg[2], d.net_pg[3], sum(d.net_pg),
                                   d.net_qg[1], d.net_qg[2], d.net_qg[3], d.worst_util,
                                   ls.total_kw, ls.neutral_kw, ls.pct_load])
            @printf("  %6.1f kVAr  %s  VUF=%6.3f%%  |V|n=%.4f-%.4f  NEV=%.4f  util=%.3f  loss=%.3f kW\n",
                    rating, scn, 100*vuf, vm.min_v, vm.max_v, nev, d.worst_util, ls.total_kw)
        else
            push!(rating_rows, Any[rating, scn, string(st), NaN, NaN, NaN, NaN,
                                   nothing, nothing, nothing, nothing,
                                   nothing, nothing, nothing, nothing,
                                   NaN, NaN, NaN])
            println("  $rating kVAr $scn: $st")
        end
    end
end

# =======================================================================
# BLOCK 3 -- LOAD-SCALE STRESS. Does the unaided network breach 2% VUF, and
# can each device bring it back?
#
# Requires VUF_CAP_PU = Inf. With the cap active the solver cannot return a
# breaching point at all, so this block would report infeasibility instead
# of the phenomenon.
# =======================================================================

scale_rows = Vector{Any}[]
scale_split_label = "$(RATING_SPLIT[1])/$(RATING_SPLIT[2])/$(RATING_SPLIT[3])"

println("\n" * "="^80)
println(" BLOCK 3 -- load-scale stress at $scale_split_label kW split, STATCOM $STRESS_KVAR kVAr")
println(" limits: bounds_active=$BOUNDS_ACTIVE  vuf_cap=$VUF_CAP_PU  (2% drawn as reference only)")
println("="^80)

for scale in LOAD_SCALE_LEVELS
    for (kind, scn, obj) in [(:none, "A", "cost"), (:qonly, "B", "VUF"), (:pexchange, "C", "VUF")]
        dm = load_base_network(data_path; enforce_bounds=BOUNDS_ACTIVE,
                               men_ohm=MEN_OHM, neutral_ref=NEUTRAL_REF)
        set_unbalance!(dm, SBASE_KVA; split_kw=RATING_SPLIT)
        # Scaling AFTER set_unbalance!, which overwrites pd/qd outright.
        for (_, l) in dm["load"]; l["pd"] *= scale; l["qd"] *= scale; end

        gids = kind == :none ? Int[] :
            add_statcoms!(dm, SBASE_KVA; n_units=N_STATCOMS,
                          rating_kvar_total=STRESS_KVAR, p_exchange=(kind == :pexchange))
        global STATCOM_GEN_IDS = gids

        st = try
            build_and_solve(dm; objective_mode=obj, neutral_ref=NEUTRAL_REF,
                            vmin=VMIN_PU, vmax=VMAX_PU,
                            bounds_active=BOUNDS_ACTIVE, vuf_cap=VUF_CAP_PU)
        catch e
            println("  scale $scale $scn: ERROR -- $e"); :ERROR
        end

        if st in [JuMP.LOCALLY_SOLVED, JuMP.OPTIMAL, JuMP.ALMOST_LOCALLY_SOLVED]
            vuf, _ = worst_case_vuf(ref_neutral=true)
            vm     = worst_case_vmag(ref_neutral=true)
            nev, _ = worst_case_nev()
            v2, _  = worst_neg_seq()
            d      = aggregate_statcom_dispatch(gids, SBASE_KVA, STRESS_KVAR)
            tot_kw, _ = total_load_kw_kvar(dm, SBASE_KVA)
            ls        = network_losses(SBASE_KVA)
            push!(scale_rows, Any[scale, tot_kw, scn, "SOLVED", vuf, vm.min_v, vm.max_v, nev, v2,
                                  vuf > VUF_LIMIT_REF,
                                  isnothing(d) ? nothing : sum(d.net_pg),
                                  isnothing(d) ? nothing : d.worst_util,
                                  ls.total_kw, ls.neutral_kw, ls.shunt_kw, ls.pct_load])
            @printf("  x%.1f (%5.1f kW)  %s  VUF=%6.3f%%%s  |V|n=%.4f-%.4f  NEV=%.4f  loss=%.3f kW (%.2f%%)%s\n",
                    scale, tot_kw, scn, 100*vuf,
                    vuf > VUF_LIMIT_REF ? " BREACH" : "      ",
                    vm.min_v, vm.max_v, nev, ls.total_kw, ls.pct_load,
                    isnothing(d) ? "" : @sprintf("  util=%.3f", d.worst_util))
        else
            push!(scale_rows, Any[scale, NaN, scn, string(st), NaN, NaN, NaN, NaN, NaN,
                                  nothing, nothing, nothing, NaN, NaN, NaN, NaN])
            println("  x$scale $scn: $st")
        end
    end
end

# =======================================================================
# OUTPUT -- CSVs only. Figures are generated separately from these files,
# so a change to a figure never requires a re-solve.
# =======================================================================
println()

write_csv(joinpath(outdir, "toy_summary.csv"),
    ["split_kw","unbalance_index","scenario","objective","status",
     "vuf_neutral","vuf_bus",
     "vmin_pu","vmax_pu","vmin_phase","vmax_phase",
     "vmin_ground_pu","vmax_ground_pu",
     "nev_pu","nev_bus","v2_max_pu",
     "loss_total_kw","loss_series_kw","loss_neutral_kw","loss_shunt_kw","loss_pct_of_load",
     "statcom_p_ph1_kw","statcom_p_ph2_kw","statcom_p_ph3_kw","statcom_p_sum_kw",
     "statcom_q_ph1_kvar","statcom_q_ph2_kvar","statcom_q_ph3_kvar",
     "worst_leg_utilisation"],
    summary_rows)

write_csv(joinpath(outdir, "toy_profiles.csv"),
    ["split_kw","unbalance_index","scenario","objective",
     "bus_id","bus_name","phase","vm_ground_pu","vm_neutral_pu","nev_pu","nev_v"],
    profile_rows)

write_csv(joinpath(outdir, "toy_rating.csv"),
    ["rating_kvar","scenario","status","vuf_neutral","vmin_pu","vmax_pu","nev_pu",
     "statcom_p_ph1_kw","statcom_p_ph2_kw","statcom_p_ph3_kw","statcom_p_sum_kw",
     "statcom_q_ph1_kvar","statcom_q_ph2_kvar","statcom_q_ph3_kvar",
     "worst_leg_utilisation",
     "loss_total_kw","loss_neutral_kw","loss_pct_of_load"],
    rating_rows)

write_csv(joinpath(outdir, "toy_loadscale.csv"),
    ["load_scale","total_kw","scenario","status","vuf_neutral","vmin_pu","vmax_pu",
     "nev_pu","v2_max_pu","breaches_2pct","statcom_p_sum_kw","worst_leg_utilisation",
     "loss_total_kw","loss_neutral_kw","loss_shunt_kw","loss_pct_of_load"],
    scale_rows)

println("\nDone. Figures: run toy_figures.jl against $outdir")
println("Sanity checks before plotting:")
println("  1. Balanced split (6/6/6), scenario A: NEV ≈ 0 and VUF ≈ 0. If not, MEN or the network is wrong.")
println("  2. statcom_p_sum_kw ≈ 0 for scenario C. That is the power-exchange signature.")
println("  3. statcom_p_ph*_kw all ≈ 0 for scenario B. If not, the p_exchange=false branch is leaking.")
println("  4. Under objective=cost, B and C differ even though neither has a Q term in the")
println("     objective. Both return arbitrary feasible points, but C's P-exchange changes the")
println("     FEASIBLE SET itself, so the arbitrary points are not the same point. The cost rows")
println("     are degenerate and should not be read as behaviour.")
println("     (Note: at the BALANCED split both cut loss 0.294 -> 0.264 kW. That is power-factor")
println("      correction on 0.95 PF loads, and it disappears under the VUF objective -- which")
println("      has no reason to correct PF. An argument for a combined objective.)")
println("  5. vmin_pu (neutral-ref) should now be LOWER than vmin_ground_pu by roughly nev_pu.")
println("     If they are equal, NEUTRAL_REF did not take effect.")
println("  6. v2_max_pu against the HARD-CODED $(V2_HARD_LIMIT) pu limit in constraints_unified.jl.")
println("     That constraint is unconditional and invisible from this script. If v2_max_pu")
println("     nears it, the NO-STATCOM case will go infeasible at higher unbalance or LOAD_SCALE,")
println("     which would look like a result but is an artefact of a hard-coded constant.")
println("  7. loss_total_kw at the BALANCED split is the floor -- no device can beat it.")
println("     Compare each scenario's loss against that floor, not against zero.")
println("  8. In toy_rating.csv, worst_leg_utilisation should fall below 1.0 as rating grows.")
println("     The rating where it unpins is the knee between size-limited and formulation-limited.")