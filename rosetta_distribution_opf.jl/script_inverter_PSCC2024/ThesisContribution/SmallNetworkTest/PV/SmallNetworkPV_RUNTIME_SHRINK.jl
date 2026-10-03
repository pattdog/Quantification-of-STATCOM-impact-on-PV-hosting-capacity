#=
==============================================================================
SMALL NETWORK PV -- hosting capacity with AS/NZS 4777.2 volt-var / volt-watt
==============================================================================
Successor to Small_Network_Sweep2.jl. Same network, same MEN treatment, same
raw-include pipeline (variables.jl -> constraints_unified.jl -> objectives),
so what Sweep2 validated carries over unchanged. What is new is the PV fleet
and its mandated control, added AFTER constraints_unified.jl by
add_pv_control!(). The shared file is not edited: it still serves the 908-bus
runs, and its own PV patch only fires for type "PV", which this script never
uses.

WHAT IT RUNS (methodology doc, sections 3 and 5)
---------------------------------------------------------------------------
  Block 0  sanity         3/3/3, unity PF: neutral current, NEV, VUF ~ 0
  Block 1  E1/E2/E3       HC0 (unity PF), HC1 (volt-var), HC1 (volt-var +
                          volt-watt), across the allocation axis. HC0 is
                          checked against the independent-solver calibration
                          (11.56 / 6.01 / 4.68 kW/lot).
  Block 2  E4/E5          HC2 = HC1 + STATCOM, one row per converter mode, and
                          the fraction of attainable (HC2-HC1)/(HC_bal-HC1).
                          Bc vs B under export is E5.
  Block 2b                The same, with the device on plain minimum-grid-
                          import control (no knowledge of the limits).
  Block 3  E7             Substitution: PV fleet reactive output with and
                          without the STATCOM at the same operating point
                          (x = HC1 of that allocation).
  Block 4                 Uniform size sweep with endbus per-phase phasors,
                          and an R/X test: why volt-var lowers voltage but
                          raises unbalance.

DESIGN DECISIONS, AND WHY
---------------------------------------------------------------------------
1. PV ACTIVE POWER IS AN EQUALITY, NOT A CAP. The BMOPFTools VVWO tutorial
   writes volt-watt as P <= f_VW(U) and relies on a zero PV price to push P
   up to the cap. That only holds under an import-cost objective. Under the
   VUF objective the optimiser would curtail heavy-phase PV to cut unbalance,
   and E4 would credit that curtailment to the STATCOM. A real inverter has
   no such discretion, so here
         P = smin( p_avail, S*f_VW(U), S*sqrt(1 - f_VV(U)^2) )     [var priority]
   and PV has NO degrees of freedom. The STATCOM is the only decision
   variable in any case, which is also why enforcing bounds in block 2 is
   legitimate: PV cannot curtail its way to feasibility.

2. SMOOTH DROOP. Each piecewise-linear curve is a sum of ReLU hinges, each
   hinge replaced by softplus: relu_e(z) = e*log(1+exp(z/e)). Error is at
   most e*log(2) per hinge, at the kink only; the telescoping hinge pairs
   cancel the bias exactly away from the kinks. This replaces the min/max
   clipped ramps in constraints_unified.jl, whose Hessian jumps at every
   breakpoint (Ipopt assumes C2). Registered as one JuMP operator with
   analytic, overflow-safe derivatives. Checked numerically before writing:
   at EPS_V_VOLTS = 0.1 V the curve error is <= 0.24% of S (volt-var) and
   <= 0.8% of S (volt-watt), and the smooth volt-var output never leaves
   [-0.60, +0.44], so sqrt(1 - f_VV^2) >= 0.8 and is always safe.
   Reference: Mhanna, Geth, Quiertant & Mancarella, IEEE TPWRS 2026;
   BMOPFTools "Smooth droop encoding" tech note.

3. CONTROLLING VOLTAGE is phase-to-NEUTRAL at the inverter's own bus and
   phase (what the inverter measures), not the phase-to-ground quantity the
   built-in bounds in constraints_unified.jl use.

4. APPARENT-POWER PRIORITY: REACTIVE POWER, PER AS/NZS 4777.2 (methodology
   1.3, resolved). The standard: "Where the inverter apparent power rating is
   reached, active power output level shall be reduced to meet the inverter
   apparent power rating while meeting the reactive power output level
   required by the power quality response mode." That is PRIORITY = :var:
   Q follows the volt-var curve exactly and P takes the remainder of the
   capability circle, sqrt(S^2 - Q^2). Consequence: with p_avail = S (a
   5 kVA inverter behind a 6.6 kW array at noon), volt-var curtails active
   power BELOW 253 V, before volt-watt engages. Curtailment = volt-watt +
   circle; counting volt-watt alone under-reports it.
   PRIORITY = :watt (P first, Q clipped to the remaining headroom) is kept
   only as a sensitivity case -- it is not the standard's behaviour.

5. AGGREGATION. The network carries one load per (bus, phase). PV follows
   the same aggregation: one single-phase unit per (bus, phase) rated
   PENETRATION x lots x S_lot. Identical lots at one bus see one voltage and
   the droop is in percent of rating, so this is exact for uniform
   penetration. Random placement draws (methodology 3.2) need per-lot buses
   and are NOT supported here.

6. HOSTING-CAPACITY TEST, per point:
     no device : bounds OFF, solve (a power flow -- PV is fully determined),
                 then check the four criteria post hoc. Thermal ratings are
                 also relaxed in the model and checked post hoc, otherwise a
                 thermal breach shows up as an opaque infeasibility.
     STATCOM   : route 1 -- bounds off, "cost" objective (header 9), check post
                 hoc. If that fails, route 2 -- bounds ENFORCED (phase-to-
                 neutral 216.2-253.0 V, VUF 2%, thermal): a solved point means
                 a compliant dispatch exists. Pass if either route passes.
                 Route 2 alone would be the definition; route 1 guards
                 against Ipopt's local infeasibility verdicts on a nonconvex
                 problem.
   Bisection on per-lot inverter rating, bracket expanded if the top passes.

7. HC IS INSTALLED kVA PER LOT (= kW/lot at PAVAIL_FRAC = 1, the clipped-
   inverter solar-noon case). Delivered active power is reported alongside,
   because under var priority volt-var itself costs active power.

8. VOLTAGE LIMITS ARE IN VOLTS (253.0 / 216.2 V, AS 60038), converted with
   each bus's own vbase. Sweep2's VMAX_PU = 1.10 is 1.10 x 230.9 = 254.0 V,
   not 253.0 V, so the two scripts' "1.10 pu" are not the same limit.

9. THE "cost" OBJECTIVE IS MINIMUM GRID IMPORT, NOT LOSS MINIMISATION.
   objectives_FIXED.jl's "cost" minimises the priced real power drawn at the
   source. Loads and PV are fixed by equalities and must be delivered; the
   only free variable is the STATCOM dispatch. So
         source P = load + losses - PV output.
   Without PV (Sweep2) load is fixed and this IS loss minimisation for a
   fixed delivery. With PV under var priority, PV output depends on voltage
   through the capability circle, so the objective also rewards pulling
   voltage down to free PV active power -- with the feeder exporting it is
   maximum net export. It is labelled "minimum grid import" here for that
   reason. (No "loss" objective is used: that branch of objectives_FIXED.jl
   indexes br_r incorrectly -- see small_network_state.md, section 4.)

10. PV CONTROL BY BLOCK. Block 1 runs every allocation three ways: :unity
   (no inverter response), :vv (volt-var only) and :vvw (volt-var + volt-
   watt). Blocks 2, 2b and 3 use HC1_MODE = :vvw -- both curves, because
   both are mandated and that is the status quo. The first run showed :vv and
   :vvw give identical HC1 (volt-watt only acts above 253 V, the HC limit),
   so this changes no HC value; it matters in block 3, whose relaxed runs
   can exceed 253 V, where a real inverter would curtail.

CHANGES AFTER THE FIRST FULL RUN (2026-09-28)
---------------------------------------------------------------------------
- Parse once: the .dss is parsed once per bounds setting and deep-copied.
- Block 2 bisection starts at HC1 (HC2 >= HC1 by construction; confirmed).
- HC_TOL 0.02 -> 0.05 kVA/lot.
- 'active' column: limits near-active at the HC point. The old 'binding'
  label for STATCOM rows came from the minimum-grid-import dispatch at the
  first failing point, not from the dispatch that defines HC2 -- misleading.
- Block 2b added: HC2 under plain cost dispatch (no route 2).
- HC1_MODE :vv -> :vvw (header 10). "Loss minimisation" relabelled
  "minimum grid import" throughout (header 9).
- E7: substitution ratio for Q-only modes only; VUF objective skipped at
  balanced allocation (flat objective, arbitrary dispatch).
- Per-stage timing printed at the end.

CHANGES 2026-10-01
---------------------------------------------------------------------------
- Core files compiled once and reused instead of include()d per solve
  (COMPILE-ONCE section). Profiling put ~63 s of every ~63 s solve in the
  three include()s; Ipopt itself takes 0.06-0.08 s. Validated against
  include() automatically before each run (VALIDATE_CORE).

- Block 4 added: uniform size sweep with endbus per-phase phasors, and an
  R/X test (line reactance scaled). Writes pv_sweep.csv.
- Substitution ratio suppressed when the device supplies < 0.1 kvar.

KNOWN LIMITS OF THIS SCRIPT
---------------------------------------------------------------------------
- Not executed by its author: written against Sweep2 and
  constraints_unified.jl without access to variables.jl or
  objectives_FIXED.jl. First run: expect to fix small interface issues
  (see the checks printed at the end).
- Uses JuMP's nonlinear interface (JuMP >= 1.15), which the existing
  min/max volt-var patch already requires.
- Snapshot only. "Curtailment" in the output is snapshot kW, not energy.

Run with:
    julia --project=<thesis root> Small_Network_PV.jl
==============================================================================
=#

# Wall clock for the whole run, package loading included (Base only, no extra packages).
RUN_T0 = time()      # not const: the value changes on every include() in one REPL session
println("Run started:  ", Libc.strftime("%Y-%m-%d %H:%M:%S", RUN_T0))

using Logging
Logging.disable_logging(Logging.Warn)

using Pkg

# --- Path anchoring: identical to Small_Network_Sweep2.jl ---------------------
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

# max_iter raised from Ipopt's default: the droop adds nonconvex equalities.
ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer,
    "print_level" => 0, "sb" => "yes", "max_iter" => 3000)

data_path = joinpath(RPMD_ROOT, "data", "4_Bus_Small.dss")
# data_path = joinpath(RPMD_ROOT, "data", "4_Bus_Small.txt")
outdir    = joinpath(THESIS_ROOT, "results", "small_network_pv")

isfile(data_path) || error("network file not found: $data_path")

# =======================================================================
# CONFIGURATION
# =======================================================================

# --- Network: unchanged from Sweep2 -------------------------------------
const ADMD_KVA   = 4.5
const LOAD_PF    = 0.95
const KW_PER_LOT = ADMD_KVA * LOAD_PF        # 4.275 kW

const ALLOCATIONS = [(3,3,3), (4,3,2), (5,3,1), (6,2,1), (7,1,1)]
const BALANCED    = (3,3,3)
const LOAD_BUS_NAME = "endbus"
const STATCOM_BUS   = "endbus"   # best-case siting; say so in every caption
const STATCOM_KVA   = 30.0

const MEN_OHM     = 10.0
const MEN_PER_LOT = true
const LOTS_AT_BUS = Dict("mid1" => 6, "mid2" => 6, "endbus" => 9)

# --- Operating point: solar noon (methodology section 2) ----------------
# Calibration used 25% of ADMD. Load and PV are one joint state: never sweep
# one of them here as if it were time of day.
const LOAD_FRAC_NOON = 0.25

# Source voltage. `nothing` keeps the value in the .dss so HC0 stays
# comparable with the calibration. HC is very sensitive to this (the VVWO
# tutorial runs its head at 1.05 pu); set it deliberately, never by accident.
const SOURCE_VM_PU = nothing

# --- PV fleet -----------------------------------------------------------
const PENETRATION = 1.0      # fraction of lots with PV, uniform (see header 5)
const PAVAIL_FRAC = 1.0      # p_avail / S at solar noon. 1.0 = clipped inverter
                             # (5 kVA behind 6.6 kW). Lower it for an unclipped
                             # array; HC is then still reported per kVA installed.

# AS/NZS 4777.2:2020 Australia A. (volts, fraction of rated S; +ve = supplying)
const VV_POINTS = [(207.0, 0.44), (220.0, 0.0), (240.0, 0.0), (258.0, -0.60)]
const VW_POINTS = [(253.0, 1.0), (260.0, 0.2)]

# Apparent-power priority -- methodology 1.3, UNRESOLVED. See header 4.
const PRIORITY = :var        # :var = AS/NZS 4777.2 (header 4); :watt = sensitivity only

# Softplus smoothing. Voltage hinges in volts; power hinges as a fraction of S.
const EPS_V_VOLTS = 0.1      # curve error <= 0.24% S (VV), <= 0.8% S (VW)
const EPS_P_FRAC  = 1e-3     # smooth-min error <= 0.07% S

# --- Hosting-capacity criteria (methodology 3.1) -------------------------
const V_MAX_V   = 253.0      # phase-to-neutral, AS 60038
const V_MIN_V   = 216.2
const VUF_LIMIT = 0.02
const V_TOL_V   = 0.05       # post-hoc tolerance on voltage checks
const VUF_TOL   = 1e-4
const THERM_TOL = 1e-3

# Bisection on per-lot inverter rating (kVA/lot).
const HC_LO  = 0.0
const HC_HI  = 20.0          # expanded x1.5 while the top still passes
const HC_CAP = 60.0          # give up expanding beyond this; reported as ">"
const HC_TOL = 0.05          # 50 W/lot: well inside the model's own accuracy

# Calibration from the independent nodal solver, unity PF, 25% load.
const CALIB_HC0 = Dict((3,3,3) => 11.56, (5,3,1) => 6.01, (7,1,1) => 4.68)

# Which blocks to run. Block 3 needs block 1's HC1 values.
const RUN_BLOCK0 = true
const RUN_BLOCK1 = true
const RUN_BLOCK2 = true
const RUN_BLOCK2_COST = true # 2b: HC2 when the STATCOM just minimises grid import
const RUN_BLOCK3 = true

# Block 4: uniform size sweep + R/X test (the mechanism behind "volt-var
# lowers voltage but raises unbalance"). No bisection: every size is solved
# and the endbus per-phase phasors are recorded.
const RUN_BLOCK4     = true
const SWEEP_ALLOCS   = [(3,3,3), (5,3,1), (7,1,1)]
const SWEEP_MODES    = [:unity, :vvw]
const SWEEP_SIZES    = collect(0.5:0.5:12.0)     # kVA per customer
const SWEEP_X_SCALES = [1.0, 3.0]                # line reactance multiplier (1.0 = as built)
const SWEEP_BUS      = "endbus"
const SWEEP_PRINT_X  = 5.0                       # size at which the summary line is printed

# Compile-once loading of the core files (see the COMPILE-ONCE section). The
# global is deliberately non-const so validate_compiled_core() can toggle it.
USE_COMPILED_CORE = true
const VALIDATE_CORE = true   # before the run: compare compiled vs include on 2 cases

const PV_MODES       = [:unity, :vv, :vvw]          # HC0, HC1, HC1 + volt-watt
const HC1_MODE       = :vvw                          # status quo: both mandated curves (header 10)
const STATCOM_MODES  = [(:qonly, "B", "Q-only"), (:qcap, "Bc", "Q-only capacitive"),
                        (:ponly, "C", "P-only"), (:pq,   "D", "P+Q")]
const STATCOM_HC_OBJECTIVE = "cost"                  # route 1: minimum grid import (header 9)
const E7_OBJECTIVES  = ["cost", "VUF"]

# -----------------------------------------------------------------------
# kW/kVAr <-> pu  (verbatim from Sweep2)
# -----------------------------------------------------------------------
kw_to_pu(p_kw, sbase_kva)     = p_kw / sbase_kva
kvar_to_pu(q_kvar, sbase_kva) = q_kvar / sbase_kva
pu_to_kw(p_pu, sbase_kva)     = p_pu * sbase_kva
pu_to_kvar(q_pu, sbase_kva)   = q_pu * sbase_kva

# =======================================================================
# CARRIED OVER FROM Small_Network_Sweep2.jl -- keep in sync with that file
# =======================================================================

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

    for (_, gen) in data_math["gen"]
        gen["pmax"] =  [1e4, 1e4, 1e4]
        gen["pmin"] = -[1e4, 1e4, 1e4]
        gen["qmax"] =  [1e4, 1e4, 1e4]
        gen["qmin"] = -[1e4, 1e4, 1e4]
    end
    return data_math
end

get_sbase_kva(dm) = dm["settings"]["sbase"] * dm["settings"]["power_scale_factor"] / 1000

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
    gen["connections"]   = [1, 2, 3, 4]
    gen["configuration"] = PMD.WYE
    if mode == :qonly
        gen["pmax"] = zeros(3);            gen["pmin"] = zeros(3)
        gen["qmax"] = fill(s_leg_pu, 3);   gen["qmin"] = -fill(s_leg_pu, 3)
        gen["statcom_p_exchange"] = false
    elseif mode == :qcap
        gen["pmax"] = zeros(3);            gen["pmin"] = zeros(3)
        gen["qmax"] = fill(s_leg_pu, 3);   gen["qmin"] = zeros(3)
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

is_virtual(bus)  = startswith(get(bus, "name", ""), "_virtual")
is_real_bus(bus) = !is_virtual(bus) && get(bus, "bus_type", 1) != 3

# Same as Sweep2 except vmin/vmax are passed in pu of each bus's own vbase.
function add_neutral_referenced_bounds!(; vmin_V=V_MIN_V, vmax_V=V_MAX_V)
    n = 0
    for (i, bus) in ref[:bus]
        is_real_bus(bus) || continue
        terms = bus["terminals"]
        4 in terms || continue
        vb_V = get(bus, "vbase", 0.4 / sqrt(3)) * 1000
        vmin, vmax = vmin_V / vb_V, vmax_V / vb_V
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

solved(st) = st in [JuMP.LOCALLY_SOLVED, JuMP.OPTIMAL, JuMP.ALMOST_LOCALLY_SOLVED]

const ALPHA = exp(im * 2/3 * pi)
const TSEQ  = 1/3 * [1 1 1; 1 ALPHA ALPHA^2; 1 ALPHA^2 ALPHA]

function worst_case_vuf()
    worst, worst_bus = 0.0, nothing
    for (i, bus) in ref[:bus]
        is_virtual(bus) && continue
        t = bus["terminals"]
        (1 in t && 2 in t && 3 in t) || continue
        v012 = TSEQ * [JuMP.value(vr[p, i]) + im * JuMP.value(vi[p, i]) for p in 1:3]
        vpos, vneg = abs(v012[2]), abs(v012[3])
        if vpos > 1e-6 && vneg / vpos > worst
            worst, worst_bus = vneg / vpos, i
        end
    end
    return worst, worst_bus
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

function network_losses(sbase_kva)
    p_series = 0.0
    for (l, branch) in ref[:branch]
        r = branch["br_r"]
        n = size(r, 1)
        a = [JuMP.value(csr[c, l]) for c in 1:n]
        b = [JuMP.value(csi[c, l]) for c in 1:n]
        p_series += a' * r * a + b' * r * b
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
    return pu_to_kw(p_series + p_shunt, sbase_kva)
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
    return (net_p=net_p, net_q=net_q, worst_util=worst_util)
end

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
# NEW -- operating point, thermal bookkeeping, source voltage
# =======================================================================

# Lots per (bus_id, phase), read from the allocated demand BEFORE scaling to
# noon. Throws on a non-integer count rather than rounding it away: a
# fractional lot count means the allocation or the .dss is not what this
# script assumes, and PV sized from it would be silently wrong.
function lots_by_bus_phase(dm, sbase_kva)
    lots = Dict{Tuple{Int,Int},Float64}()
    for (lid, load) in dm["load"]
        ph = filter(c -> c != 4, load["connections"])
        length(ph) == 1 || error("load $lid: expected single-phase, got $(load["connections"])")
        n = pu_to_kw(sum(load["pd"]), sbase_kva) / KW_PER_LOT
        abs(n - round(n)) < 1e-6 || error("load $lid: $(round(n, digits=4)) lots is not an integer")
        key = (load["load_bus"], ph[1])
        lots[key] = get(lots, key, 0.0) + round(n)
    end
    return lots
end

function scale_loads!(dm, frac)
    for (_, load) in dm["load"]
        load["pd"] = load["pd"] .* frac
        load["qd"] = load["qd"] .* frac
    end
end

# R/X test (block 4): multiply every real line's reactance matrix by k, self
# and mutual alike, so the phase-plus-neutral loop reactance scales by k and
# its resistance is untouched. The virtual source branch is left alone.
function scale_line_reactance!(dm, k)
    for (_, br) in dm["branch"]
        startswith(get(br, "name", ""), "_virtual") && continue
        br["br_x"] = br["br_x"] .* k
    end
end

function set_source_voltage!(dm, vm_pu)
    isnothing(vm_pu) && return
    for (_, bus) in dm["bus"]
        bus["bus_type"] == 3 || continue
        bus["vm"][1:3] .= vm_pu
    end
end

# Ratings recorded before relaxation so they can be checked post hoc.
const THERMAL_RATINGS = Dict{Int,Vector{Float64}}()

function record_and_maybe_relax_thermal!(dm; relax::Bool)
    empty!(THERMAL_RATINGS)
    for (k, br) in dm["branch"]
        THERMAL_RATINGS[parse(Int, k)] = collect(Float64, br["c_rating_a"])
        if relax
            br["c_rating_a"] = fill(Inf, length(br["c_rating_a"]))
            haskey(br, "rate_a") && (br["rate_a"] = fill(Inf, length(br["rate_a"])))
        end
    end
end

# =======================================================================
# NEW -- PV fleet
# =======================================================================

# One single-phase unit per (bus, phase), sized PENETRATION x lots x S_lot.
# type "PV_VVW" deliberately does NOT match the "PV" guard in
# constraints_unified.jl, so that file's clipped-ramp patch never fires and
# the only PV control in the model is add_pv_control! below.
#
# qmin/qmax must be FINITE: constraints_unified.jl only links qg to the
# terminal currents for phases with a finite q bound. They are set loose
# (1.2 S) so they never bind -- the droop equality is the real constraint.
function add_pv_fleet!(dm, sbase_kva; s_lot_kva, lots, mode::Symbol)
    ids = Int[]
    for ((bus_id, ph), n) in sort(collect(lots))
        n_eff = PENETRATION * n
        n_eff > 0 || continue
        S   = kw_to_pu(n_eff * s_lot_kva, sbase_kva)
        pav = S * PAVAIL_FRAC

        gid = maximum(parse.(Int, collect(keys(dm["gen"])))) + 1
        g = deepcopy(dm["gen"]["1"])
        g["index"]         = gid
        g["gen_bus"]       = bus_id
        g["type"]          = "PV_VVW"
        g["name"]          = "pv_$(bus_id)_$(ph)"
        g["cost"]          = [0.0, 0.0]
        g["connections"]   = [ph, 4]
        g["configuration"] = PMD.WYE
        g["pmin"] = [0.0];       g["pmax"] = [2.0 * pav]
        g["qmin"] = [-1.2 * S];  g["qmax"] = [1.2 * S]
        g["pv_s"]      = S
        g["pv_pavail"] = pav
        g["pv_mode"]   = String(mode)
        g["pv_phase"]  = ph
        g["pv_lots"]   = n_eff
        dm["gen"][string(gid)] = g
        push!(ids, gid)
    end
    return ids
end

# --- smooth hinge: overflow-safe softplus and its analytic derivatives -----
_sp(t::Real)   = t > 0 ? t + log1p(exp(-t)) : log1p(exp(t))
_dsp(t::Real)  = t >= 0 ? 1 / (1 + exp(-t)) : exp(t) / (1 + exp(t))
_d2sp(t::Real) = (s = _dsp(t); s * (1 - s))

# Exact piecewise-linear curve, clamped flat outside the breakpoints. Used
# only for post-solve droop verification.
function pwl_exact(U, pts)
    U <= pts[1][1]   && return pts[1][2]
    U >= pts[end][1] && return pts[end][2]
    for k in 1:length(pts)-1
        (x0, y0), (x1, y1) = pts[k], pts[k+1]
        U <= x1 && return y0 + (y1 - y0) * (U - x0) / (x1 - x0)
    end
end

# Smooth curve as a JuMP expression: y1 + sum over sloped segments of
# slope * (relu(U - x0) - relu(U - x1)). The pair telescopes, so the
# softplus bias cancels exactly away from the kinks.
function curve_expr(U, pts_pu, relu)
    f = pts_pu[1][2]
    for k in 1:length(pts_pu)-1
        (x0, y0), (x1, y1) = pts_pu[k], pts_pu[k+1]
        s = (y1 - y0) / (x1 - x0)
        s == 0 && continue
        f = f + s * (relu(U - x0) - relu(U - x1))
    end
    return f
end

# Monitored voltage per PV unit, kept for post-solve verification.
const PV_VM = Dict{Int,Any}()

# Call AFTER include(constraints_unified.jl), so pg/qg are already tied to
# the terminal currents. Adds, per PV unit:
#   :unity  Q = 0,            P = p_avail
#   :vv     Q = S f_VV(U),    P = smin(p_avail, S sqrt(1 - f_VV^2))            [var]
#   :vvw    Q = S f_VV(U),    P = smin(p_avail, S f_VW(U), S sqrt(1 - f_VV^2)) [var]
# with the :watt variants in the else branch.
function add_pv_control!()
    empty!(PV_VM)
    pv = [(id, g) for (id, g) in ref[:gen] if get(g, "type", "") == "PV_VVW"]
    isempty(pv) && return 0

    op = JuMP.add_nonlinear_operator(model, 1, _sp, _dsp, _d2sp; name = :pv_softplus)

    for (id, g) in pv
        ph, bus_id = g["pv_phase"], g["gen_bus"]
        S, pav, mode = g["pv_s"], g["pv_pavail"], g["pv_mode"]

        if mode == "unity"
            JuMP.@constraint(model, qg[ph, id] == 0.0)
            JuMP.@constraint(model, pg[ph, id] == pav)
            continue
        end

        bus = ref[:bus][bus_id]
        4 in bus["terminals"] || error("PV at bus $bus_id has no neutral terminal")
        vb_V = get(bus, "vbase", 0.4 / sqrt(3)) * 1000

        # Phase-to-neutral magnitude. Start at 1.0 so the defining equality
        # never starts at the origin, where its gradient vanishes.
        U = JuMP.@variable(model, base_name = "U_pv_$(id)", lower_bound = 0.0, start = 1.0)
        JuMP.@constraint(model, U^2 == (vr[ph, bus_id] - vr[4, bus_id])^2 +
                                       (vi[ph, bus_id] - vi[4, bus_id])^2)

        eV = EPS_V_VOLTS / vb_V
        eP = EPS_P_FRAC * S
        reluV(z) = eV * op(z / eV)
        reluP(z) = eP * op(z / eP)
        smin(a, b) = a - reluP(a - b)

        vv_pu = [(v / vb_V, y) for (v, y) in VV_POINTS]
        vw_pu = [(v / vb_V, y) for (v, y) in VW_POINTS]
        fvv = curve_expr(U, vv_pu, reluV)

        # Active-power ceiling before the capability circle.
        p_head = mode == "vvw" ? smin(pav, S * curve_expr(U, vw_pu, reluV)) : min(pav, S)

        if PRIORITY == :var
            Qexpr = S * fvv
            Pexpr = smin(p_head, S * sqrt(1 - fvv^2))   # 1 - fvv^2 >= 0.64 always
        elseif PRIORITY == :watt
            Pexpr = p_head
            qh    = sqrt(S^2 - Pexpr^2 + eP^2)          # eP^2 keeps sqrt off zero
            qc    = S * fvv
            Qexpr = qc - reluP(qc - qh) + reluP(-qh - qc)
        else
            error("PRIORITY must be :var or :watt")
        end

        JuMP.@constraint(model, qg[ph, id] == Qexpr)
        JuMP.@constraint(model, pg[ph, id] == Pexpr)
        PV_VM[id] = U
    end
    return length(pv)
end

# Post-solve: fleet totals, per-phase Q, and the droop verification demanded
# by methodology section 6. droop_err is the largest deviation of any unit's
# (P, Q) from the EXACT piecewise curve at its solved voltage, as a fraction
# of its rating. It should be of the order of the smoothing error, not more.
function pv_dispatch(sbase_kva)
    P = Q = PAV = 0.0
    q_ph = zeros(3)
    err = 0.0
    for (id, g) in ref[:gen]
        get(g, "type", "") == "PV_VVW" || continue
        ph, S, pav = g["pv_phase"], g["pv_s"], g["pv_pavail"]
        p = JuMP.value(pg[ph, id]); q = JuMP.value(qg[ph, id])
        P += p; Q += q; PAV += pav; q_ph[ph] += q

        if g["pv_mode"] == "unity"
            p_ex, q_ex = pav, 0.0
        else
            bus = ref[:bus][g["gen_bus"]]
            U_V = JuMP.value(PV_VM[id]) * get(bus, "vbase", 0.4 / sqrt(3)) * 1000
            fvv = pwl_exact(U_V, VV_POINTS)
            p_head = g["pv_mode"] == "vvw" ? min(pav, S * pwl_exact(U_V, VW_POINTS)) : min(pav, S)
            if PRIORITY == :var
                q_ex = S * fvv
                p_ex = min(p_head, S * sqrt(1 - fvv^2))
            else
                p_ex = p_head
                qh   = sqrt(max(S^2 - p_ex^2, 0.0))
                q_ex = clamp(S * fvv, -qh, qh)
            end
        end
        S > 0 && (err = max(err, abs(p - p_ex) / S, abs(q - q_ex) / S))
    end
    return (p_kw = pu_to_kw(P, sbase_kva), q_kvar = pu_to_kvar(Q, sbase_kva),
            pavail_kw = pu_to_kw(PAV, sbase_kva), curt_kw = pu_to_kw(PAV - P, sbase_kva),
            q_ph_kvar = pu_to_kvar.(q_ph, sbase_kva), droop_err = err)
end

# =======================================================================
# NEW -- hosting-capacity criteria
# =======================================================================

bus_name(i) = isnothing(i) ? "" : get(ref[:bus][i], "name", string(i))

function vpn_extremes_V()
    vmax, vmax_at = -Inf, ""
    vmin, vmin_at = Inf, ""
    for (i, bus) in ref[:bus]
        is_real_bus(bus) || continue
        t = bus["terminals"]
        vb_V = get(bus, "vbase", 0.4 / sqrt(3)) * 1000
        vrn = 4 in t ? JuMP.value(vr[4, i]) : 0.0
        vin = 4 in t ? JuMP.value(vi[4, i]) : 0.0
        for p in 1:3
            p in t || continue
            v = sqrt((JuMP.value(vr[p, i]) - vrn)^2 + (JuMP.value(vi[p, i]) - vin)^2) * vb_V
            v > vmax && ((vmax, vmax_at) = (v, "$(bus_name(i)) ph$p"))
            v < vmin && ((vmin, vmin_at) = (v, "$(bus_name(i)) ph$p"))
        end
    end
    return (vmax=vmax, vmax_at=vmax_at, vmin=vmin, vmin_at=vmin_at)
end

# Worst conductor utilisation |I| / rating, against the ratings recorded
# before relaxation. Series current, as in Sweep2's branch_currents.
function thermal_worst(sbase_kva)
    worst, at = 0.0, ""
    in_max_a = 0.0
    for (l, br) in ref[:branch]
        startswith(get(br, "name", ""), "_virtual") && continue
        rating = get(THERMAL_RATINGS, l, Float64[])
        ibase = sbase_kva / get(ref[:bus][br["f_bus"]], "vbase", 0.4 / sqrt(3))   # A
        n = length(br["f_connections"])
        for c in 1:n
            I = hypot(JuMP.value(csr[c, l]), JuMP.value(csi[c, l]))
            c == 4 && (in_max_a = max(in_max_a, I * ibase))
            c <= length(rating) && isfinite(rating[c]) && rating[c] > 0 || continue
            u = I / rating[c]
            u > worst && ((worst, at) = (u, "$(get(br, "name", string(l))) c$c"))
        end
    end
    return (util=worst, at=at, in_max_a=in_max_a)
end

function evaluate_point(gen_ids, sbase_kva, statcom_kva)
    v = vpn_extremes_V()
    vuf, vuf_bus = worst_case_vuf()
    nev, _ = worst_case_nev()
    th = thermal_worst(sbase_kva)
    pvd = pv_dispatch(sbase_kva)
    sd = statcom_dispatch(gen_ids, sbase_kva, statcom_kva)
    return (vmax_V=v.vmax, vmax_at=v.vmax_at, vmin_V=v.vmin, vmin_at=v.vmin_at,
            vuf=vuf, vuf_at=bus_name(vuf_bus), nev=nev,
            therm=th.util, therm_at=th.at, in_max_a=th.in_max_a,
            loss_kw=network_losses(sbase_kva),
            pv=pvd, statcom=sd)
end

# Violations as (criterion, relative excess), largest first. Empty = passes.
function violations(m)
    v = Tuple{String,Float64}[]
    m.vmax_V > V_MAX_V + V_TOL_V     && push!(v, ("vmax "  * m.vmax_at, (m.vmax_V - V_MAX_V) / V_MAX_V))
    m.vmin_V < V_MIN_V - V_TOL_V     && push!(v, ("vmin "  * m.vmin_at, (V_MIN_V - m.vmin_V) / V_MIN_V))
    m.vuf    > VUF_LIMIT + VUF_TOL   && push!(v, ("vuf "   * m.vuf_at,  (m.vuf - VUF_LIMIT) / VUF_LIMIT))
    m.therm  > 1 + THERM_TOL         && push!(v, ("therm " * m.therm_at, m.therm - 1))
    sort!(v, by = x -> -x[2])
    return v
end

# -----------------------------------------------------------------------
# Per-phase phasors at one bus, for block 4.
#
# Angles are reported as the SHIFT from each phase's nominal position
# (0, -120, +120 degrees), so a perfectly balanced bus reads 0/0/0 whatever
# its magnitude. Two references, because the two limits see different things:
#   vpn : phase-to-NEUTRAL -- what the 253 V limit and the inverter measure
#   vpg : phase-to-GROUND  -- what VUF is computed from (a neutral shift is
#         zero-sequence and cannot appear in |V2|/|V1|)
# v1/v2 are the positive- and negative-sequence magnitudes in volts.
# -----------------------------------------------------------------------
function bus_id_by_name(name)
    for (i, bus) in ref[:bus]
        get(bus, "name", "") == name && return i
    end
    error("bus '$name' not found in the solved model")
end

_wrap180(a) = mod(a + 180.0, 360.0) - 180.0
const NOMINAL_ANGLE_DEG = (0.0, -120.0, 120.0)

function bus_phasors(name, sbase_kva)
    i   = bus_id_by_name(name)
    bus = ref[:bus][i]
    vb  = get(bus, "vbase", 0.4 / sqrt(3)) * 1000
    vg  = [(JuMP.value(vr[p, i]) + im * JuMP.value(vi[p, i])) * vb for p in 1:3]
    vn  = 4 in bus["terminals"] ? (JuMP.value(vr[4, i]) + im * JuMP.value(vi[4, i])) * vb : 0.0im
    vpn = vg .- vn
    shift(v, p) = _wrap180(rad2deg(angle(v)) - NOMINAL_ANGLE_DEG[p])
    v012 = TSEQ * vg
    pv_p, pv_q = zeros(3), zeros(3)
    for (id, g) in ref[:gen]
        (get(g, "type", "") == "PV_VVW" && g["gen_bus"] == i) || continue
        ph = g["pv_phase"]
        pv_p[ph] += pu_to_kw(JuMP.value(pg[ph, id]), sbase_kva)
        pv_q[ph] += pu_to_kvar(JuMP.value(qg[ph, id]), sbase_kva)
    end
    return (vpn_mag = abs.(vpn), vpn_ang = [shift(vpn[p], p) for p in 1:3],
            vpg_mag = abs.(vg),  vpg_ang = [shift(vg[p], p)  for p in 1:3],
            nev = abs(vn), v1 = abs(v012[2]), v2 = abs(v012[3]), pv_p = pv_p, pv_q = pv_q)
end

# =======================================================================
# BUILD + SOLVE -- Sweep2's pipeline plus add_pv_control!
# =======================================================================
function build_and_solve_pv(dm; objective_mode::String="cost", bounds_active::Bool=false,
                            vuf_cap=Inf)
    t0 = time()
    pm = PMD.instantiate_mc_model(
        dm, PMD.IVRENPowerModel, (pm) -> nothing;
        jump_model = JuMP.Model(ipopt_solver), ref_extensions = Function[])

    global ref   = pm.ref[:it][:pmd][:nw][0]
    global model = pm.model
    global VUF_CAP = vuf_cap

    if USE_COMPILED_CORE
        _core_variables!()                       # compiled once, reused (see COMPILE-ONCE)
        _core_constraints!()
    else
        include(joinpath(CORE_DIR, "variables.jl"))
        include(joinpath(CORE_DIR, "constraints_unified.jl"))
    end

    add_pv_control!()
    bounds_active && add_neutral_referenced_bounds!()

    global objective = objective_mode
    global objective_target_gens = [i for (i, g) in ref[:gen] if get(g, "type", "") == "STATCOM"]
    if USE_COMPILED_CORE
        _core_objectives!()
    else
        include(joinpath(CORE_DIR, "objectives_FIXED.jl"))
    end

    TIMING[:build] += time() - t0
    t1 = time()
    JuMP.optimize!(model)
    TIMING[:optimize] += time() - t1
    TIMING[:ipopt]    += JuMP.solve_time(model)
    TIMING[:n]        += 1
    return JuMP.termination_status(model)
end

# Per-stage wall time, accumulated over the run and printed at the end.
#   parse    : .dss parse + transform (should be ~0 after the first two calls)
#   build    : instantiate + include(variables/constraints/objectives) + droop
#   optimize : optimize!, of which ipopt is Ipopt's own reported time; the
#              difference is JuMP/MOI model passing and nonlinear setup.
const TIMING = Dict{Symbol,Float64}(:parse => 0.0, :build => 0.0, :optimize => 0.0,
                                    :ipopt => 0.0, :n => 0.0)

# The .dss is parsed once per bounds setting and deep-copied for every solve.
# Nothing downstream of load_base_network depends on the case, so this is
# identical to re-parsing, just without the cost.
const _BASE_CACHE = Dict{Bool,Any}()
function base_network(bounds_active::Bool)
    if !haskey(_BASE_CACHE, bounds_active)
        t = time()
        _BASE_CACHE[bounds_active] = load_base_network(data_path; enforce_bounds=bounds_active,
            men_ohm=MEN_OHM, neutral_ref=true,
            lots_at_bus=(MEN_PER_LOT ? LOTS_AT_BUS : nothing))
        TIMING[:parse] += time() - t
    end
    t = time()
    dm = deepcopy(_BASE_CACHE[bounds_active])
    TIMING[:parse] += time() - t
    return dm
end

function run_pv_case(; alloc, s_lot_kva, pv_mode::Symbol, statcom_mode::Symbol=:none,
                     statcom_kva=STATCOM_KVA, objective="cost", bounds_active=false,
                     x_scale=1.0)
    dm = base_network(bounds_active)
    x_scale == 1.0 || scale_line_reactance!(dm, x_scale)
    set_source_voltage!(dm, SOURCE_VM_PU)
    record_and_maybe_relax_thermal!(dm; relax = !bounds_active)
    set_phase_allocation!(dm, SBASE_KVA; bus_name=LOAD_BUS_NAME, lots_per_phase=alloc)
    lots = lots_by_bus_phase(dm, SBASE_KVA)
    scale_loads!(dm, LOAD_FRAC_NOON)
    s_lot_kva > 0 && add_pv_fleet!(dm, SBASE_KVA; s_lot_kva=s_lot_kva, lots=lots, mode=pv_mode)

    gids = statcom_mode == :none ? Int[] :
        add_statcom!(dm, SBASE_KVA; bus_name=STATCOM_BUS,
                     rating_kva_total=statcom_kva, mode=statcom_mode)

    st = try
        build_and_solve_pv(dm; objective_mode=objective, bounds_active=bounds_active,
                           vuf_cap = bounds_active ? VUF_LIMIT : Inf)
    catch e
        println("      ERROR: ", sprint(showerror, e)); :ERROR
    end
    return st, gids
end

# =======================================================================
# HOSTING CAPACITY -- one point, then bisection
# =======================================================================

const POINT_HEADER = [
    "block", "allocation", "pv_mode", "statcom", "x_kva_per_lot", "route", "status", "pass", "binding",
    "vmax_v", "vmax_at", "vmin_v", "vmin_at", "vuf_pct", "vuf_at", "nev_pu", "therm_util", "therm_at",
    "i_neutral_max_a", "loss_kw", "pv_p_kw", "pv_pavail_kw", "pv_curt_kw", "pv_q_kvar",
    "pv_q_ph1_kvar", "pv_q_ph2_kvar", "pv_q_ph3_kvar", "droop_err_frac",
    "sc_p1_kw", "sc_p2_kw", "sc_p3_kw", "sc_q1_kvar", "sc_q2_kvar", "sc_q3_kvar", "sc_util"]

point_rows = Vector{Any}[]

function point_row(block, alloc, pv_mode, sc, x, route, st, pass, binding, m)
    head = Any[block, alloc_label(alloc), String(pv_mode), sc, x, route, string(st), pass, binding]
    isnothing(m) && return Any[head; fill(nothing, length(POINT_HEADER) - length(head))]
    s = m.statcom
    sv(f, k) = isnothing(s) ? nothing : f(s)[k]
    return Any[head...,
        m.vmax_V, m.vmax_at, m.vmin_V, m.vmin_at, 100 * m.vuf, m.vuf_at, m.nev, m.therm, m.therm_at,
        m.in_max_a, m.loss_kw, m.pv.p_kw, m.pv.pavail_kw, m.pv.curt_kw, m.pv.q_kvar,
        m.pv.q_ph_kvar..., m.pv.droop_err,
        sv(x -> x.net_p, 1), sv(x -> x.net_p, 2), sv(x -> x.net_p, 3),
        sv(x -> x.net_q, 1), sv(x -> x.net_q, 2), sv(x -> x.net_q, 3),
        isnothing(s) ? nothing : s.worst_util]
end

# Does installed size x (kVA/lot) pass every criterion? See header 6.
# enforce_fallback=false skips route 2: the answer is then "does the STATCOM,
# running STATCOM_HC_OBJECTIVE with no knowledge of the limits, keep the
# network compliant?" -- the cost-dispatch HC of block 2b.
function hc_point(alloc, x; pv_mode, statcom_mode=:none, block="", enforce_fallback=true)
    sc = statcom_mode == :none ? "none" : String(statcom_mode)
    relaxed_only = statcom_mode == :none || !enforce_fallback

    st, gids = run_pv_case(alloc=alloc, s_lot_kva=x, pv_mode=pv_mode, statcom_mode=statcom_mode,
                           objective=STATCOM_HC_OBJECTIVE, bounds_active=false)
    if solved(st)
        m = evaluate_point(gids, SBASE_KVA, STATCOM_KVA)
        v = violations(m)
        pass = isempty(v)
        push!(point_rows, point_row(block, alloc, pv_mode, sc, x, "relaxed", st, pass,
                                    pass ? "" : v[1][1], m))
        (pass || relaxed_only) && return (pass=pass, binding=pass ? "" : v[1][1],
                                          status=string(st), m=m)
        binding = v[1][1]
    else
        push!(point_rows, point_row(block, alloc, pv_mode, sc, x, "relaxed", st, false,
                                    "not solved", nothing))
        binding = "relaxed solve: $st"
        relaxed_only && return (pass=false, binding=binding, status=string(st), m=nothing)
    end

    # Route 2: a compliant STATCOM dispatch exists iff this solves.
    st2, gids2 = run_pv_case(alloc=alloc, s_lot_kva=x, pv_mode=pv_mode, statcom_mode=statcom_mode,
                             objective=STATCOM_HC_OBJECTIVE, bounds_active=true)
    if solved(st2)
        m2 = evaluate_point(gids2, SBASE_KVA, STATCOM_KVA)
        v2 = violations(m2)
        pass2 = isempty(v2)
        push!(point_rows, point_row(block, alloc, pv_mode, sc, x, "enforced", st2, pass2,
                                    pass2 ? "" : v2[1][1], m2))
        return (pass=pass2, binding=pass2 ? "" : v2[1][1], status=string(st2), m=m2)
    end
    push!(point_rows, point_row(block, alloc, pv_mode, sc, x, "enforced", st2, false,
                                "infeasible: $binding", nothing))
    return (pass=false, binding=binding, status=string(st2), m=nothing)
end

# Constraints within a small margin of their limit at the HC point -- what
# actually stops the search, as opposed to `binding`, which is the worst
# violation at the first FAILING point. For STATCOM rows those differ: the
# failing point's violation comes from the relaxed (minimum-grid-import) dispatch,
# which is not the dispatch that defines HC2. Margins are loose enough to
# cover the HC_TOL gap between the reported HC and the true boundary.
function active_set(m)
    isnothing(m) && return ""
    a = String[]
    m.vmax_V >= V_MAX_V - 0.5          && push!(a, "vmax")
    m.vmin_V <= V_MIN_V + 0.5          && push!(a, "vmin")
    m.vuf    >= VUF_LIMIT - 5e-4       && push!(a, "vuf")
    m.therm  >= 0.98                   && push!(a, "therm")
    !isnothing(m.statcom) && m.statcom.worst_util >= 0.98 && push!(a, "device")
    return isempty(a) ? "none" : join(a, "+")
end

# Largest x that passes. Assumes a single pass/fail crossing; a bracket that
# fails at the bottom is reported rather than bisected.
#
# lo_start: a size already known to pass. Block 2 passes HC1: every STATCOM
# mode can dispatch zero, so HC2 >= HC1 by construction -- and the first run
# confirmed it (no HC2 below HC1; C at 3/3/3 equal to it). If lo_start
# unexpectedly fails, the search falls back to HC_LO rather than trusting it.
function hc_bisect(alloc; pv_mode, statcom_mode=:none, block="", lo_start=HC_LO,
                   enforce_fallback=true)
    n = 0
    pt(x) = (n += 1; hc_point(alloc, x; pv_mode, statcom_mode, block, enforce_fallback))

    r_lo = pt(lo_start)
    if !r_lo.pass && lo_start > HC_LO
        println("      note: lo_start = $(round(lo_start, digits=2)) failed ($(r_lo.binding)); restarting from $(HC_LO)")
        lo_start = HC_LO
        r_lo = pt(lo_start)
    end
    r_lo.pass || return (hc=NaN, binding="fails at x=$(lo_start): $(r_lo.binding)", n=n,
                         capped=false, active="", m=nothing)

    lo, m_lo = lo_start, r_lo.m
    hi = max(HC_HI, lo_start + 2.0)
    r_hi = pt(hi)
    while r_hi.pass
        lo, m_lo = hi, r_hi.m
        hi >= HC_CAP && return (hc=hi, binding="passes at cap", n=n, capped=true,
                                active=active_set(m_lo), m=m_lo)
        hi = min(1.5 * hi, HC_CAP)
        r_hi = pt(hi)
    end
    binding = r_hi.binding
    while hi - lo > HC_TOL
        mid = (lo + hi) / 2
        r = pt(mid)
        if r.pass
            lo, m_lo = mid, r.m
        else
            hi, binding = mid, r.binding
        end
    end
    return (hc=lo, binding=binding, n=n, capped=false, active=active_set(m_lo), m=m_lo)
end

fmt_hc(r) = isnan(r.hc) ? "   --  " : (r.capped ? @sprintf(">%6.2f", r.hc) : @sprintf("%7.2f", r.hc))

# =======================================================================
# COMPILE-ONCE LOADING OF THE CORE FILES
#
# Profiling (2026-09-29) showed ~99.8% of every solve's time in the three
# include()s -- variables.jl 28 s, constraints_unified.jl 22 s,
# objectives_FIXED.jl 13 s -- against 0.06 s for add_pv_control!() and
# 0.08 s for Ipopt. include() re-lowers and recompiles a file's top-level
# code on every call, and the closures/generators in the JuMP macros leave
# new compiled code behind each time, so a long session also slows down.
#
# Fix, WITHOUT editing the core files (they stay the single source of truth,
# and other scripts that include() them are unaffected): read each file once
# at start-up, wrap its parsed contents in a function, and call that function
# per solve. The first call compiles; later calls reuse the compiled code.
#
# Scoping is the only subtlety, handled mechanically:
#   * Names a file assigns at TOP LEVEL -- including inside top-level
#     if/elseif blocks, which do not introduce a scope -- are globals under
#     include(), and the other files and the reporting code read them (vr,
#     pg, crg_bus, ...). They are declared `global` in the wrapper.
#   * Names assigned only inside loops or comprehensions are locals under
#     include() (soft scope) and stay locals in the wrapper.
#   * Top-level function definitions (objectives_FIXED.jl's four helpers) are
#     evaluated once at global scope, which also stops them being redefined
#     on every solve.
#   * One harmless difference: constraints_unified.jl reuses `alpha` as a
#     load exponent inside its load loop. Under include() that copy is
#     local; in the wrapper it overwrites the global sequence operator. T,
#     Tre and Tim are built from alpha before that loop, and
#     objectives_FIXED.jl rebuilds all four, so no result can change.
#
# validate_compiled_core() solves two cases both ways and requires the
# objective and every bus voltage to agree to 1e-7 before the run starts.
# Set USE_COMPILED_CORE = false to fall back to include() exactly as before.
# If you edit a core file mid-REPL-session, re-run load_compiled_core!().
# =======================================================================

_is_funcdef(ex) = ex isa Expr && (ex.head === :function || ex.head === :macro ||
    (ex.head === :(=) && ex.args[1] isa Expr && ex.args[1].head in (:call, :where)))

function _lhs_names!(names::Set{Symbol}, lhs)
    if lhs isa Symbol
        push!(names, lhs)
    elseif lhs isa Expr && lhs.head === :tuple
        foreach(a -> _lhs_names!(names, a), lhs.args)
    elseif lhs isa Expr && lhs.head === :(::)
        _lhs_names!(names, lhs.args[1])
    end                      # x[i] = ..., x.f = ... mutate, they do not bind
    return names
end

# Collect names bound in the same scope as the statement. Descends only into
# constructs that do NOT open a new scope (block, if/elseif); stops at
# for/while/let/function/comprehension/generator, whose bindings are local.
function _collect_assigned!(names::Set{Symbol}, ex)
    ex isa Expr || return names
    h = ex.head
    if h in (:(=), Symbol("+="), Symbol("-="), Symbol("*="), Symbol("/="))
        _is_funcdef(ex) || _lhs_names!(names, ex.args[1])
        _collect_assigned!(names, ex.args[2])         # chained a = b = ...
    elseif h in (:block, :if, :elseif, :const)
        foreach(a -> _collect_assigned!(names, a), ex.args)
    end
    return names
end

function _compile_core_file(path::AbstractString, fname::Symbol)
    top = Meta.parseall(read(path, String); filename = path)
    body = Any[]
    nfun = 0
    for st in top.args
        if _is_funcdef(st)
            Core.eval(Main, st)              # helper functions: defined once, globally
            nfun += 1
        else
            push!(body, st)
        end
    end
    names = Set{Symbol}()
    foreach(st -> _collect_assigned!(names, st), body)
    decl = isempty(names) ? Any[] : Any[Expr(:global, sort!(collect(names))...)]
    Core.eval(Main, Expr(:function, Expr(:call, fname),
                         Expr(:block, decl..., body..., :(return nothing))))
    return (globals = length(names), helpers = nfun)
end

function load_compiled_core!()
    for (file, fname) in (("variables.jl", :_core_variables!),
                          ("constraints_unified.jl", :_core_constraints!),
                          ("objectives_FIXED.jl", :_core_objectives!))
        info = _compile_core_file(joinpath(CORE_DIR, file), fname)
        println("  compiled-once: $(rpad(file, 24)) $(info.globals) globals",
                info.helpers > 0 ? ", $(info.helpers) helper functions hoisted" : "")
    end
end

function _core_snapshot()
    v = Dict{Tuple{Int,Int},ComplexF64}()
    for (i, bus) in ref[:bus], t in bus["terminals"]
        v[(i, t)] = JuMP.value(vr[t, i]) + im * JuMP.value(vi[t, i])
    end
    return JuMP.objective_value(model), v
end

# Same case built both ways; objective and every terminal voltage must agree.
# Two objectives so both branches of objectives_FIXED.jl that this script
# uses ("cost", "VUF") are exercised, with a P+Q STATCOM so the S^2 and DC-bus
# patches in constraints_unified.jl are too.
function validate_compiled_core(; alloc=(7,1,1), x=5.0, tol=1e-7)
    global USE_COMPILED_CORE
    keep = USE_COMPILED_CORE
    try
        for obj in ("cost", "VUF")
            snap = Dict{Bool,Any}()
            for compiled in (false, true)
                USE_COMPILED_CORE = compiled
                st, _ = run_pv_case(alloc=alloc, s_lot_kva=x, pv_mode=:vvw,
                                    statcom_mode=:pq, objective=obj)
                solved(st) || error("validation solve failed ($obj, compiled=$compiled): $st")
                snap[compiled] = _core_snapshot()
            end
            dobj = abs(snap[false][1] - snap[true][1])
            dv   = maximum(abs(snap[false][2][k] - snap[true][2][k]) for k in keys(snap[false][2]))
            @printf("  validate [%-4s]  |Δobjective| = %.2e   max |ΔV| = %.2e pu   %s\n",
                    obj, dobj, dv, (dobj < tol && dv < tol) ? "OK" : "MISMATCH")
            (dobj < tol && dv < tol) ||
                error("compiled core differs from include() -- set USE_COMPILED_CORE = false and report")
        end
    finally
        USE_COMPILED_CORE = keep
    end
end

# =======================================================================
# RUN
# =======================================================================

USE_COMPILED_CORE && load_compiled_core!()

_probe = load_base_network(data_path; men_ohm=MEN_OHM,
                           lots_at_bus=(MEN_PER_LOT ? LOTS_AT_BUS : nothing))
const SBASE_KVA = get_sbase_kva(_probe)
set_source_voltage!(_probe, SOURCE_VM_PU)
set_phase_allocation!(_probe, SBASE_KVA; bus_name=LOAD_BUS_NAME, lots_per_phase=BALANCED)
_lots = lots_by_bus_phase(_probe, SBASE_KVA)
_src_vm = [b["vm"][1:3] for (_, b) in _probe["bus"] if b["bus_type"] == 3]

println("=" ^ 92)
println(" SMALL NETWORK PV -- hosting capacity, AS/NZS 4777.2 Australia A")
println(" sbase = $SBASE_KVA kVA   |   solar noon: load = $(Int(100*LOAD_FRAC_NOON))% of ADMD")
println(" source vm (pu) = $(_src_vm)" * (isnothing(SOURCE_VM_PU) ? "  [from .dss]" : "  [overridden]"))
println(" lots per (bus, phase) at 3/3/3: ", join(["$(get(_probe["bus"][string(b)], "name", b)) ph$p=$(Int(n))"
                                              for ((b, p), n) in sort(collect(_lots))], ", "))
println(" PV: penetration $(PENETRATION), p_avail = $(PAVAIL_FRAC) x S, priority = :$PRIORITY")
println(" smoothing: $(EPS_V_VOLTS) V on voltage hinges, $(EPS_P_FRAC) x S on power hinges")
println(" criteria: $(V_MIN_V)-$(V_MAX_V) V phase-to-neutral, VUF $(100*VUF_LIMIT)%, conductor rating")
println(" STATCOM: $(STATCOM_KVA) kVA at $STATCOM_BUS (best-case siting)")
println(" core files: ", USE_COMPILED_CORE ? "compiled once" : "include() per solve")
println("=" ^ 92)

if USE_COMPILED_CORE && VALIDATE_CORE
    println("\n Validating compiled core against include() (slow once: the include() side)...")
    validate_compiled_core()
    foreach(k -> TIMING[k] = 0.0, collect(keys(TIMING)))   # time the run, not the check
end

hc_rows = Vector{Any}[]
e7_rows = Vector{Any}[]
HC = Dict{Tuple{Any,Symbol,Symbol},Any}()   # (alloc, pv_mode, statcom_mode) => result

# -----------------------------------------------------------------------
# BLOCK 0 -- sanity. Balanced allocation, unity PF, a moderate size.
# Neutral current, NEV and VUF must all be ~0; droop_err must be ~0.
# -----------------------------------------------------------------------
if RUN_BLOCK0
    println("\n BLOCK 0 -- sanity: 3/3/3, unity PF, 5 kVA/lot, and 3/3/3 volt-var at 10 kVA/lot")
    for (mode, x) in ((:unity, 5.0), (:vv, 10.0))
        st, gids = run_pv_case(alloc=BALANCED, s_lot_kva=x, pv_mode=mode)
        if solved(st)
            m = evaluate_point(gids, SBASE_KVA, STATCOM_KVA)
            v = violations(m)
            push!(point_rows, point_row("0", BALANCED, mode, "none", x, "relaxed", st, isempty(v),
                                        isempty(v) ? "" : v[1][1], m))
            @printf("   %-6s x=%5.1f  Vmax=%6.1f V  VUF=%.2e  NEV=%.2e pu  In_max=%.2e A  P=%6.1f kW  Q=%+6.1f kvar  droop_err=%.1e\n",
                    mode, x, m.vmax_V, m.vuf, m.nev, m.in_max_a, m.pv.p_kw, m.pv.q_kvar, m.pv.droop_err)
        else
            println("   $mode x=$x: $st  <-- fix this before trusting anything below")
        end
    end
end

# -----------------------------------------------------------------------
# BLOCK 1 -- E1 / E2 / E3. No device. HC0, HC1(vv), HC1(vv+vw).
# -----------------------------------------------------------------------
if RUN_BLOCK1
    println("\n" * "=" ^ 92)
    println(" BLOCK 1 -- HC ladder, no device   [kVA/customer; 'active' = limits at the HC point]")
    println("=" ^ 92)
    for alloc in ALLOCATIONS
        for mode in PV_MODES
            r = hc_bisect(alloc; pv_mode=mode, block="1")
            HC[(alloc, mode, :none)] = r
            cal = get(CALIB_HC0, alloc, NaN)
            push!(hc_rows, Any["1", alloc_label(alloc), String(mode), "none", "A",
                               r.hc, r.capped, r.active, r.binding, r.n,
                               mode == :unity ? cal : nothing, nothing, nothing])
            @printf("  %-6s %-6s HC = %s   active: %-16s first fail: %-24s (%2d)%s\n",
                    alloc_label(alloc), mode, fmt_hc(r), r.active, r.binding, r.n,
                    (mode == :unity && !isnan(cal)) ? @sprintf("   calib %.2f", cal) : "")
        end
    end
end

# -----------------------------------------------------------------------
# BLOCK 2 -- E4 / E5. HC2 = HC1 + STATCOM, per mode.
# Fraction of attainable = (HC2 - HC1) / (HC_bal - HC1), HC_bal = HC1 at 3/3/3.
# -----------------------------------------------------------------------
if RUN_BLOCK2
    println("\n" * "=" ^ 92)
    println(" BLOCK 2 -- HC2 on top of HC1 (:$HC1_MODE), $(STATCOM_KVA) kVA at $STATCOM_BUS")
    println("=" ^ 92)
    hc_bal = haskey(HC, (BALANCED, HC1_MODE, :none)) ? HC[(BALANCED, HC1_MODE, :none)].hc :
             hc_bisect(BALANCED; pv_mode=HC1_MODE, block="2").hc
    for alloc in ALLOCATIONS
        hc1 = haskey(HC, (alloc, HC1_MODE, :none)) ? HC[(alloc, HC1_MODE, :none)].hc :
              (HC[(alloc, HC1_MODE, :none)] = hc_bisect(alloc; pv_mode=HC1_MODE, block="2")).hc
        @printf("\n  %-6s HC1 = %6.2f   ceiling HC_bal = %6.2f\n", alloc_label(alloc), hc1, hc_bal)
        for (smode, code, label) in STATCOM_MODES
            r = hc_bisect(alloc; pv_mode=HC1_MODE, statcom_mode=smode, block="2", lo_start=hc1)
            HC[(alloc, HC1_MODE, smode)] = r
            gap  = hc_bal - hc1
            frac = gap > HC_TOL ? (r.hc - hc1) / gap : NaN
            push!(hc_rows, Any["2", alloc_label(alloc), String(HC1_MODE), String(smode), code,
                               r.hc, r.capped, r.active, r.binding, r.n, nothing, r.hc - hc1, frac])
            @printf("    %-2s %-18s HC2 = %s   dHC = %+6.2f   frac = %s   active: %-20s (%2d)\n",
                    code, label, fmt_hc(r), r.hc - hc1,
                    isnan(frac) ? "  -- " : @sprintf("%5.2f", frac), r.active, r.n)
        end
    end
end

# -----------------------------------------------------------------------
# BLOCK 2b -- HC2 under COST DISPATCH: relaxed route only, no fallback.
# Block 2 asks whether ANY compliant dispatch exists (ideal, limit-aware
# control). This asks what the device delivers when it simply minimises
# grid import (the "cost" objective, header 9), with no knowledge of the
# voltage or VUF limits. The first run showed that at 6/2/1 and 7/1/1,
# 5 kVA/lot, a Q-only device on this objective pushes VUF over 2% on a
# network that is compliant with no device at all -- so this HC can fall
# BELOW HC1, and the search starts from zero, not from HC1.
# -----------------------------------------------------------------------
if RUN_BLOCK2_COST
    println("\n" * "=" ^ 92)
    println(" BLOCK 2b -- HC2 under cost dispatch ([$STATCOM_HC_OBJECTIVE], limits not known to the device)")
    println("=" ^ 92)
    for alloc in ALLOCATIONS
        hc1 = haskey(HC, (alloc, HC1_MODE, :none)) ? HC[(alloc, HC1_MODE, :none)].hc : NaN
        @printf("\n  %-6s HC1 = %6.2f\n", alloc_label(alloc), hc1)
        for (smode, code, label) in STATCOM_MODES
            r = hc_bisect(alloc; pv_mode=HC1_MODE, statcom_mode=smode, block="2b",
                          enforce_fallback=false)
            push!(hc_rows, Any["2b", alloc_label(alloc), String(HC1_MODE), String(smode), code,
                               r.hc, r.capped, r.active, r.binding, r.n, nothing, r.hc - hc1, nothing])
            @printf("    %-2s %-18s HC2 = %s   dHC = %+6.2f   %s  first fail: %-24s (%2d)\n",
                    code, label, fmt_hc(r), r.hc - hc1, r.hc < hc1 - HC_TOL ? "WORSE THAN NO DEVICE" : "                    ",
                    r.binding, r.n)
        end
    end
end

# -----------------------------------------------------------------------
# BLOCK 3 -- E7, substitution. At x = HC1(alloc), bounds relaxed, compare
# the PV fleet's reactive output with and without the STATCOM.
#   displaced   = |Q_pv, no device| - |Q_pv, with device|        (kvar)
#   subst_ratio = displaced / sum|Q_statcom|
# A ratio near 1 means every kvar the STATCOM supplied was a kvar the PV
# fleet stopped supplying: displacement, not addition.
#
# The ratio is computed for Q-ONLY modes (B, Bc) only. For C and D part of
# the displacement is caused by real-power exchange, so dividing by the
# device's Q alone is meaningless (the first run gave D a ratio of 1.36).
# Their displacement is reported in kvar, which is the number that matters.
#
# The VUF objective is skipped at balanced allocation: VUF is already ~0
# there, the objective is flat, and the dispatch Ipopt returns is arbitrary.
#
# Under var priority, freed Q also frees P through the circle: dP_pv. This
# is why the "cost" objective is minimum grid import rather than loss
# minimisation in PV cases (header 9).
# -----------------------------------------------------------------------
if RUN_BLOCK3
    println("\n" * "=" ^ 92)
    println(" BLOCK 3 -- E7 substitution at x = HC1(:$HC1_MODE)   [bounds relaxed]")
    println("=" ^ 92)
    for alloc in ALLOCATIONS
        haskey(HC, (alloc, HC1_MODE, :none)) || continue
        x = HC[(alloc, HC1_MODE, :none)].hc
        (isnan(x) || x <= 0) && continue

        st0, _ = run_pv_case(alloc=alloc, s_lot_kva=x, pv_mode=HC1_MODE)
        solved(st0) || (println("  $(alloc_label(alloc)) baseline: $st0"); continue)
        m0 = evaluate_point(Int[], SBASE_KVA, STATCOM_KVA)
        push!(point_rows, point_row("3", alloc, HC1_MODE, "none", x, "relaxed", st0, true, "", m0))
        @printf("\n  %-6s x = %5.2f kVA/lot   no device: Q_pv = %+7.2f kvar  P_pv = %6.2f kW  Vmax = %6.1f V\n",
                alloc_label(alloc), x, m0.pv.q_kvar, m0.pv.p_kw, m0.vmax_V)

        for (smode, code, label) in STATCOM_MODES, obj in E7_OBJECTIVES
            (obj == "VUF" && alloc == BALANCED) && continue
            st, gids = run_pv_case(alloc=alloc, s_lot_kva=x, pv_mode=HC1_MODE,
                                   statcom_mode=smode, objective=obj)
            if !solved(st)
                push!(e7_rows, Any[alloc_label(alloc), x, code, obj, string(st), fill(nothing, 9)...])
                println("    $code [$obj]: $st"); continue
            end
            m = evaluate_point(gids, SBASE_KVA, STATCOM_KVA)
            push!(point_rows, point_row("3", alloc, HC1_MODE, String(smode), x, "relaxed-$obj", st, true, "", m))
            sq  = sum(abs, m.statcom.net_q)
            dis = abs(m0.pv.q_kvar) - abs(m.pv.q_kvar)
            rat = (smode in (:qonly, :qcap) && sq > 0.1) ? dis / sq : NaN   # < 0.1 kvar: device idle, ratio is noise
            push!(e7_rows, Any[alloc_label(alloc), x, code, obj, "SOLVED",
                               m0.pv.q_kvar, m.pv.q_kvar, dis, sq, rat,
                               m.pv.p_kw - m0.pv.p_kw, m0.vmax_V, m.vmax_V, 100 * m.vuf])
            @printf("    %-2s [%-4s] Q_pv = %+7.2f  displaced = %+6.2f kvar  sum|Q_sc| = %6.2f  ratio = %s  dP_pv = %+5.2f kW  Vmax = %6.1f V\n",
                    code, obj, m.pv.q_kvar, dis, sq, isnan(rat) ? " -- " : @sprintf("%4.2f", rat),
                    m.pv.p_kw - m0.pv.p_kw, m.vmax_V)
        end
    end
end

# -----------------------------------------------------------------------
# BLOCK 4 -- uniform size sweep and R/X test. No device.
#
# Question: WHY does volt-var lower the worst voltage but raise VUF at the
# same PV size? Proposed answer: on a resistive loop (R/X ~ 2.6 here),
# absorbing Q on the heavy phase reduces that phase's voltage MAGNITUDE only
# by ~X|Q| but shifts its ANGLE by ~R|Q|, and the absorption is single-phase.
# The 253 V limit sees magnitude; VUF sees the whole phasor.
#
# Two tests of that answer, both read off pv_sweep.csv:
#   (a) Phasors. Under volt-var, the heavy phase should show LOWER magnitude
#       but a LARGER angle shift than at unity PF, at the same size.
#   (b) R/X. With the line reactance scaled up (SWEEP_X_SCALES), Q becomes a
#       better magnitude lever and a weaker angle lever, so the VUF penalty
#       of volt-var (dVUF = VUF_vvw - VUF_unity) should shrink or change sign.
# If (b) does not happen, the R/X explanation is wrong and needs replacing.
# -----------------------------------------------------------------------
const SWEEP_HEADER = [
    "allocation", "pv_mode", "x_scale", "x_kva_per_lot", "status",
    "vmax_v", "vmin_v", "vuf_pct", "v1_v", "v2_v", "pv_p_kw", "pv_q_kvar", "pv_curt_kw",
    "vpn_mag_a", "vpn_mag_b", "vpn_mag_c", "vpn_ang_a", "vpn_ang_b", "vpn_ang_c",
    "vpg_mag_a", "vpg_mag_b", "vpg_mag_c", "vpg_ang_a", "vpg_ang_b", "vpg_ang_c",
    "nev_v", "bus_pv_p_a", "bus_pv_p_b", "bus_pv_p_c", "bus_pv_q_a", "bus_pv_q_b", "bus_pv_q_c"]

sweep_rows = Vector{Any}[]

if RUN_BLOCK4
    println("\n" * "=" ^ 92)
    println(" BLOCK 4 -- uniform size sweep + R/X test at $SWEEP_BUS   [no device, limits relaxed]")
    println(" summary line at x = $(SWEEP_PRINT_X) kVA/customer; phase a = the heavy phase; angles are shifts from nominal")
    println("=" ^ 92)
    for xs in SWEEP_X_SCALES, alloc in SWEEP_ALLOCS
        probe = Dict{Symbol,Any}()
        for mode in SWEEP_MODES, x in SWEEP_SIZES
            st, gids = run_pv_case(alloc=alloc, s_lot_kva=x, pv_mode=mode, x_scale=xs)
            if !solved(st)
                push!(sweep_rows, Any[alloc_label(alloc), String(mode), xs, x, string(st),
                                      fill(nothing, length(SWEEP_HEADER) - 5)...])
                continue
            end
            m  = evaluate_point(gids, SBASE_KVA, STATCOM_KVA)
            ph = bus_phasors(SWEEP_BUS, SBASE_KVA)
            push!(sweep_rows, Any[alloc_label(alloc), String(mode), xs, x, "SOLVED",
                m.vmax_V, m.vmin_V, 100 * m.vuf, ph.v1, ph.v2, m.pv.p_kw, m.pv.q_kvar, m.pv.curt_kw,
                ph.vpn_mag..., ph.vpn_ang..., ph.vpg_mag..., ph.vpg_ang...,
                ph.nev, ph.pv_p..., ph.pv_q...])
            x == SWEEP_PRINT_X && (probe[mode] = (m, ph))
        end
        if haskey(probe, :unity) && haskey(probe, :vvw)
            (mu, pu_), (mv, pv_) = probe[:unity], probe[:vvw]
            @printf("  X x%.0f  %-6s unity: Vmax %6.1f V  VUF %5.2f%%  |Va| %6.1f V  ang_a %+5.2f deg | vvw: Vmax %6.1f V  VUF %5.2f%%  |Va| %6.1f V  ang_a %+5.2f deg | dVUF %+5.2f pp\n",
                    xs, alloc_label(alloc),
                    mu.vmax_V, 100 * mu.vuf, pu_.vpn_mag[1], pu_.vpg_ang[1],
                    mv.vmax_V, 100 * mv.vuf, pv_.vpn_mag[1], pv_.vpg_ang[1],
                    100 * (mv.vuf - mu.vuf))
        end
    end
end

# =======================================================================
# OUTPUT
# =======================================================================
println()
write_csv(joinpath(outdir, "pv_hc.csv"),
    ["block", "allocation", "pv_mode", "statcom", "scenario", "hc_kva_per_lot", "capped",
     "active_at_hc", "first_fail", "n_solves", "calib_hc0", "delta_hc", "frac_attainable"],
    hc_rows)
write_csv(joinpath(outdir, "pv_points.csv"), POINT_HEADER, point_rows)
write_csv(joinpath(outdir, "pv_substitution.csv"),
    ["allocation", "x_kva_per_lot", "scenario", "objective", "status",
     "q_pv_none_kvar", "q_pv_with_kvar", "displaced_kvar", "statcom_abs_q_kvar", "subst_ratio",
     "dp_pv_kw", "vmax_none_v", "vmax_with_v", "vuf_with_pct"],
    e7_rows)
RUN_BLOCK4 && write_csv(joinpath(outdir, "pv_sweep.csv"), SWEEP_HEADER, sweep_rows)

println("""

Checks before believing any number above:
  1. Block 0: VUF, NEV ~ 0 at 3/3/3 and droop_err ~ 1e-3 or below. If droop_err is large
     the droop constraints are not the ones being enforced -- check that constraints_unified.jl
     linked qg for the PV gens (finite qmin/qmax) and that no "PV"-type patch fired.
  2. Block 1, unity rows: HC0 against the calibration (11.56 / 6.01 / 4.68). A mismatch is
     either the source voltage, the load fraction, or the 253.0 V vs 254.0 V limit (header 8).
  3. Block 1, vv vs vvw: if volt-watt only acts above 253 V, the two HC1 columns should agree
     to within the smoothing error. A material gap means volt-watt is acting below 253 V
     (softplus leakage -- lower EPS_V_VOLTS) or the worst bus has no inverter.
  4. droop_err_frac in pv_points.csv for every solved row. The equality guarantees it; verify.
  5. Block 2: no P-ONLY (C) HC2 may exceed HC_bal (methodology 3.5 bounds real-power exchange
     only). Q-capable modes CAN exceed it when HC_bal is thermal-limited: they supply the PV
     fleet's volt-var absorption locally and relieve the head line (first run: 14.10 vs 13.91).
  6. Block 2 'active' column: a genuine boundary shows several limits active at once, usually
     with the device saturated (e.g. vmax+vuf+device). A pass with 'none' active next to a
     LOCALLY_INFEASIBLE point is the pattern to re-check by hand.
  7. Block 2b vs 2: the gap is what limit-aware control is worth over minimum grid import.
     'WORSE THAN NO DEVICE' rows are a finding, not a bug -- confirm one in pv_points.csv.
  8. Block 3: per-phase sc_p of C/D sum to ~0. Ratio column is for B/Bc only.
  9. pv_curt_kw in pv_points.csv below 253 V is capability-circle curtailment (reactive
     priority, AS/NZS 4777.2) -- real, mandated, and invisible to a volt-watt-only count.
  10. Block 4: at 3/3/3 every angle shift should be equal across phases and dVUF ~ 0. At
     7/1/1, X x1: volt-var should show LOWER |Va| but a LARGER ang_a than unity. At X x3 the
     dVUF column should shrink or change sign; if it does not, the R/X explanation is wrong.
""")

n = max(TIMING[:n], 1)
@printf("Timing over %d solves (s/solve):  parse+copy %.2f   build %.2f   optimize %.2f (Ipopt %.2f)\n",
        Int(TIMING[:n]), TIMING[:parse]/n, TIMING[:build]/n, TIMING[:optimize]/n, TIMING[:ipopt]/n)

let dt = time() - RUN_T0
    println("Run finished: ", Libc.strftime("%Y-%m-%d %H:%M:%S", time()))
    @printf("Total wall time: %.1f s  (%d min %04.1f s), package loading and compilation included\n",
            dt, floor(Int, dt / 60), dt % 60)
end