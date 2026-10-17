#=
==============================================================================
LV12_248 OPF -- hosting capacity on the real feeder, staged and self-checking
==============================================================================
Same pipeline as Small_Network_PV.jl (variables.jl -> constraints_unified.jl ->
add_pv_control! -> objectives_FIXED.jl, core files compiled once), applied to
LV12_248 (Geth et al. 2025) with measured SGSC loads (draw 45, 13:00).

Every generic function below is copied VERBATIM from Small_Network_PV.jl
(add_men!, add_statcom!, add_neutral_referenced_bounds!, add_pv_control!,
pv_dispatch, evaluate_point, the bisection, the compile-once loader, ...).
What is new is only what the 4-bus layout cannot express: one PV unit per
customer, loads from the SGSC files, and the staged checks.

STAGES -- each prints PASS/FAIL against OpenDSS and stops at the first failure
---------------------------------------------------------------------------
  1  parse Master_pmd.dss       248 buses, 247 four-wire lines, 121 loads,
                                56.829 kW, source neutral grounded
  2  earthing, source           121 MEN bonds = load buses = men_buses.csv;
                                source impedance in the math model = SRC_Z1/Z0
  3  power flow, no PV          source P/Q, every bus V_pn and NEV vs OpenDSS
  4  121 x 3.0 kW PV, unity PF  the same, near the voltage limit
  5  HC0 (uniform, unity PF)    3.025 kW/customer (OpenDSS, same 253.05 V test)
  6  HC1 (volt-var + volt-watt) new result; droop error checked
  7  HC2 per STATCOM mode       new result; 30 kVA at b733
  8  siting sweep               HC2 for B/C/D at 7 buses (lv12_siting.csv)
  9  one vs two devices         60 kVA at b3278 vs 2 x 30 kVA, one per branch (lv12_pairs.csv)
 10  unity-PF baseline          STATCOM on top of HC0 instead of HC1 (lv12_base0.csv)
 11  sequential adoption        55% existing 5 kVA fleet, then random applicants; HC = how many more
                                connect, per configuration (lv12_adoption.csv; quick test with
                                LV12_NINIT=1 LV12_NNEW=2)

SETTINGS (environment variables)
  LV12_PEN=0.55         existing penetration (nested fleets; outputs in results/pen55/)
  LV12_CONFIGS=HC1,...  run only these stage-11 configurations (parallel terminals)
  LV12_THERMAL=0        switch off thermal (400 A/conductor) + 500 kVA transformer limits
  LV12_THERMAL_A=400    assumed conductor rating, amps (feeder cables)
  LV12_LEAD_A=721.7     transformer LV leads l_4140 / switch_1916 / l_514 (default = 500 kVA
                        per-phase current; 400 reproduces the old uniform-400 A runs)
  LV12_NINIT/LV12_NNEW  number of fleets / applicant orders per fleet
  LV12_LOCAL_VSET=240   voltage set point (V) of the local_vbal_vv controller
  LV12_LOCAL_VV_W=1.0   weight of its voltage-regulation term
  STAGE 11 configurations: HC0, HC1, HC2_D_one_coord/_local, HC2_D_pair_coord/_local,
    HC2_D_one/pair_vbal, _vbalvv, HC2_D_head_coord
    (coord = minimise the feeder's highest voltage, upper bound, network limits enforced if
     needed; every local_* = the device sees only its own connection and NO network limit
     is imposed on it -- pass/fail is judged afterwards:
       local   |V2|^2 at the PCC          vbal  |V0pn|^2+|V2|^2 at the PCC
       vbalvv  vbal + V1 regulation above LV12_LOCAL_VSET;  head = device at b800)
    e.g. LV12_CONFIGS=HC1,HC2_D_one_coord,HC2_D_one_vbal,HC2_D_one_vbalvv
SOURCE SENSITIVITY: LV12_SOURCE=isc100k|tx_dataset|tx_4pct julia --project=<root> lv12_opf.jl
  (default stiff; other sources write to results_<source>/ and skip the OpenDSS port checks)

Stages 1-5 are the port check: nothing new is claimed until they pass.

WHERE TO PUT IT
---------------------------------------------------------------------------
Put the whole lv12_sgsc/ folder next to Small_Network_PV.jl (anywhere below
the folder holding core/variables.jl and below your Project.toml). The core
files and the project are found exactly as Small_Network_PV.jl finds them.
Override with ENV["THESIS_CORE_DIR"] / ENV["THESIS_PROJECT"] if needed.

    julia --project=<thesis root> lv12_sgsc/lv12_opf.jl

DECISIONS, AND WHY
---------------------------------------------------------------------------
1. THERMAL LIMITS (14 Oct): ON by default with an ASSUMED 400 A per conductor on
   every branch (the dataset has no ratings), plus the dataset's 500 kVA
   transformer as a source apparent-power limit. Caution: at 55% (67 x 5 kW,
   unity PF) the circuit head cables l_514 / l_4140 already carry ~400-500 A in
   OpenDSS, so 400 A binds at the baseline in most fleets. LV12_THERMAL=0 turns
   both off (voltage + VUF only); LV12_THERMAL_A changes the rating.
2. SOURCE: stiff (OpenDSS default impedance, 0.08 mOhm), neutral solidly
   earthed with the 4_Bus_Small.dss near-zero reactor. The impedance is SET
   in Julia (SRC_Z1/SRC_Z0, set_source_impedance!) and checked in stage 2:
   PMD's parser otherwise computes OpenDSS's default fault level at 115 kV
   (1.60 + j6.42 Ohm) and the loaded feeder has no solution. Your earlier LV12
   folder used ISC3=100000 (2.3 mOhm): that gives HC0 ~2.5% lower, and is the
   whole 2.71-vs-2.78 gap. Keeping the dataset's transformer resistance
   (3 mOhm) gives ~13% lower. Open question for Fred.
3. MEN: 10 ohm per customer; one customer per bus, so no paralleling.
4. PV: one single-phase unit per customer, on the customer's own bus and
   phase, p_avail = S (clipped inverter at noon). Uniform size here;
   count-based HC with placement draws comes after stage 5 passes.
5. CRITERIA as Small_Network_PV.jl: 216.2-253.0 V phase-to-neutral (+0.05 V
   post-hoc tolerance), VUF 2%. STATCOM HC: route 1 (cost dispatch, checked
   post hoc) then route 2 (bounds enforced).

Not executed by its author (no Julia where it was written). Expect small
interface fixes on the first run; the stage checks say where.
==============================================================================
=#

RUN_T0 = time()
println("Run started:  ", Libc.strftime("%Y-%m-%d %H:%M:%S", RUN_T0))

using Logging
Logging.disable_logging(Logging.Warn)
using Pkg

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

const PSCC_DIR = haskey(ENV, "THESIS_CORE_DIR") ? dirname(rstrip(ENV["THESIS_CORE_DIR"], ['/', '\\'])) :
                 _find_up(SCRIPT_DIR, d -> isfile(joinpath(d, "core", "variables.jl")))
const THESIS_ROOT = get(ENV, "THESIS_PROJECT",
                 _find_up(SCRIPT_DIR, d -> isfile(joinpath(d, "Project.toml")); outermost=true))
isnothing(PSCC_DIR)    && error("could not find core/variables.jl above $SCRIPT_DIR (set ENV[\"THESIS_CORE_DIR\"])")
isnothing(THESIS_ROOT) && error("could not find a Project.toml above $SCRIPT_DIR (set ENV[\"THESIS_PROJECT\"])")
const CORE_DIR = joinpath(PSCC_DIR, "core")
println("  core   : $CORE_DIR")
println("  project: $THESIS_ROOT")

Pkg.activate(THESIS_ROOT)
using rosetta_distribution_opf
import PowerModelsDistribution
import InfrastructureModels
using Ipopt
using JuMP
using Printf
using Statistics
using Random

const PMD  = PowerModelsDistribution
const RPMD = rosetta_distribution_opf
const IM   = InfrastructureModels
PMD.silence!()

ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer,
    "print_level" => 0, "sb" => "yes", "max_iter" => 3000)

# --- Files (all relative to this script) -------------------------------
const NET_FILE   = joinpath(SCRIPT_DIR, "reference_rebuild", "Master_pmd.dss")
const REF_CSV    = joinpath(SCRIPT_DIR, "reference_rebuild", "opendss_reference_h13_draw045.csv")
const MEN_CSV    = joinpath(SCRIPT_DIR, "reference_rebuild", "men_buses.csv")
const PROF_CSV   = joinpath(SCRIPT_DIR, "out", "sgsc_profiles_hourly.csv")
const ASSIGN_CSV = joinpath(SCRIPT_DIR, "out", "lv12_assignment.csv")
# Source model (see SOURCE_CASES below). "stiff" is the validated default; any other
# choice writes to results_<source>/ and skips the OpenDSS port checks (stages 3-4,
# stage-5 gate), whose targets were computed with the stiff source.
const SOURCE     = get(ENV, "LV12_SOURCE", "stiff")
const PORT_CHECK = SOURCE == "stiff"
# LV12_CONFIGS = "HC1,HC2_D_one" runs only those stage-11 configurations (parallel runs in
# separate terminals); each such run writes to its own subfolder so outputs do not clash.
const ADOPT_ONLY = filter(!isempty, strip.(split(get(ENV, "LV12_CONFIGS", ""), ",")))
const PEN_TAG    = haskey(ENV, "LV12_PEN") ? "pen" * string(round(Int, 100 * parse(Float64, ENV["LV12_PEN"]))) : ""
const THERM_TAG  = get(ENV, "LV12_THERMAL", "1") == "1" ? "" : "nothermal"
const OUTDIR     = joinpath(SCRIPT_DIR, PORT_CHECK ? "results" : "results_" * SOURCE, PEN_TAG, THERM_TAG,
                            isempty(ADOPT_ONLY) ? "" : "adopt_" * join(ADOPT_ONLY, "_"))
for f in (NET_FILE, REF_CSV, MEN_CSV, PROF_CSV, ASSIGN_CSV)
    isfile(f) || error("missing input: $f")
end

# =======================================================================
# CONFIGURATION
# =======================================================================

# Loads. Master_pmd.dss already carries draw 45 at 13:00; any other
# (draw, hour) is applied from the SGSC files by set_loads!.
const DSS_DRAW = 45
const DSS_HOUR = 13
const LOAD_DRAW = 45
const HOUR      = 13
const LOAD_PF   = 0.95

const MEN_OHM = 10.0                 # per customer; one customer per bus

# Source impedance (ohm, sequence values) -- OpenDSS's default fault level at
# 0.4 kV (MVAsc3 2000, MVAsc1 2100), i.e. a stiff source. Set explicitly by
# set_source_impedance! because PMD's parser, given no impedance in the DSS,
# computes that default at OpenDSS's default 115 kV instead: 1.60 + j6.42 ohm,
# which cannot carry 57 kW (stage 3 came back LOCALLY_INFEASIBLE).
# Sensitivity cases (ohm, Z1 and Z0), selected with ENV["LV12_SOURCE"]:
#   stiff       OpenDSS default at 0.4 kV (validated; stated in Master_pmd.dss)
#   isc100k     ISC3 = ISC1 = 100 kA, as in 4_Bus_Small.dss and Patrick's earlier LV12 folder
#   tx_dataset  the dataset's 500 kVA transformer: %r 0.4+0.4, xhl 0.039% (3.0 + j0.15 mOhm)
#   tx_4pct     the same R with a typical 4% impedance (3.0 + j14.7 mOhm)
# Transformer cases use Z0 = Z1 (Dyn unit, LV side).
const SOURCE_CASES = Dict(
    "stiff"      => (complex(1.94028500029066e-5, 7.76114000116266e-5),
                     complex(2.17289779069744e-5, 6.51869337209233e-5)),
    "isc100k"    => (complex(5.60112033611204e-4, 2.24044813444482e-3),
                     complex(7.31733159547457e-4, 2.19519947864237e-3)),
    "tx_dataset" => (complex(3.0e-3, 1.4696e-4), complex(3.0e-3, 1.4696e-4)),
    "tx_4pct"    => (complex(3.0e-3, 1.4696e-2), complex(3.0e-3, 1.4696e-2)))
haskey(SOURCE_CASES, SOURCE) || error("LV12_SOURCE = $SOURCE; choose one of $(sort(collect(keys(SOURCE_CASES))))")
const SRC_Z1 = SOURCE_CASES[SOURCE][1]
const SRC_Z0 = SOURCE_CASES[SOURCE][2]

const STATCOM_BUS = "b733"           # where HC0 binds most often across draws
const STATCOM_KVA = 30.0

const PAVAIL_FRAC = 1.0              # p_avail / S at solar noon (clipped inverter)

# AS/NZS 4777.2:2020 Australia A -- identical to Small_Network_PV.jl
const VV_POINTS = [(207.0, 0.44), (220.0, 0.0), (240.0, 0.0), (258.0, -0.60)]
const VW_POINTS = [(253.0, 1.0), (260.0, 0.2)]
const PRIORITY  = :var
const EPS_V_VOLTS = 0.1
const EPS_P_FRAC  = 1e-3

# Hosting-capacity criteria -- identical to Small_Network_PV.jl
const V_MAX_V   = 253.0
const V_MIN_V   = 216.2
const VUF_LIMIT = 0.02
const V_TOL_V   = 0.05
const VUF_TOL   = 1e-4
const THERM_TOL = 1e-3
# Thermal limits (header, decision 1). The dataset has no conductor ratings: every
# conductor of every branch (phases and neutral) gets THERMAL_A, an assumption.
# The transformer limit TX_KVA is the dataset's own 500 kVA rating. Both are
# criteria when THERMAL_ON[] (LV12_THERMAL=0 switches them off; LV12_THERMAL_A sets A).
const THERMAL_ON = Ref(get(ENV, "LV12_THERMAL", "1") == "1")
const THERMAL_A  = parse(Float64, get(ENV, "LV12_THERMAL_A", "400"))
const TX_KVA     = 500.0
# Transformer LV leads (15 Oct): l_4140, switch_1916_closed and l_514 (0.25 m, 1 m, 0.24 m,
# b908 -> b1122 -> b238 -> b800) carry BOTH circuits' current, i.e. they are the transformer's
# tails, not feeder cable. They are rated at the transformer's per-phase current,
# 500 kVA / (3 x 230.94 V) = 721.7 A, not at THERMAL_A. LV12_LEAD_A=400 reproduces the old
# (uniform 400 A) behaviour as a sensitivity.
const TX_LEADS   = Set(["l_4140", "switch_1916_closed", "l_514"])
const TX_LEAD_A  = parse(Float64, get(ENV, "LV12_LEAD_A", string(round(TX_KVA * 1000 / (3 * 230.94), digits=1))))
branch_rating_a(br) = lowercase(get(br, "name", "")) in TX_LEADS ? TX_LEAD_A : THERMAL_A

# Bisection on uniform kVA per customer
const HC_LO  = 0.0
const HC_HI  = 6.0
const HC_CAP = 20.0
const HC_TOL = 0.01

const STAGES       = Set([1, 2, 3, 4, 5, 6, 11])   # 7-10 = uniform-size device studies (done 13 Oct); 11 = adoption, refuse-and-continue (~54 solves per order per config)
const STOP_ON_FAIL = true

USE_COMPILED_CORE = true
const VALIDATE_CORE = false   # true = compare compiled vs include() once (slow: include() path)

const HC1_MODE      = :vvw
const STATCOM_MODES = [(:qonly, "B", "Q-only"), (:qcap, "Bc", "Q-only capacitive"),
                       (:ponly, "C", "P-only"), (:pq,   "D", "P+Q")]
const STATCOM_HC_OBJECTIVE = get(ENV, "LV12_DISPATCH", "minvmax")   # "minvmax" (HC-maximising) or "cost" (loss-minimising, as Small_Network_PV.jl)

# OpenDSS targets: reference_rebuild/Master_pmd_opendss_check.dss, draw 45, 13:00
const T_LOAD_KW   = 56.829   # sum of Loads_h13_draw045.dss
const T_LOAD_KVAR = 18.679
const T_NOPV = (src_kw=57.304,   src_kvar=18.842, vmin=238.853, vmax=245.070, nev=2.211, vuf=0.318)
const T_PV3  = (src_kw=-300.411, src_kvar=21.398, vmin=243.299, vmax=252.966, nev=3.122, vuf=0.582)
const T_HC0  = 3.025          # OpenDSS bisection with the same 253.05 V pass test
const TOL_V   = 0.10   # V
const TOL_NEV = 0.05   # V
const TOL_P   = 0.10   # kW / kvar
const TOL_VUF = 0.01   # percentage points
const TOL_HC  = 0.02   # kVA/customer (bisection step is 0.01)

# -----------------------------------------------------------------------
# kW/kVAr <-> pu  (verbatim)
# -----------------------------------------------------------------------
kw_to_pu(p_kw, sbase_kva)     = p_kw / sbase_kva
kvar_to_pu(q_kvar, sbase_kva) = q_kvar / sbase_kva
pu_to_kw(p_pu, sbase_kva)     = p_pu * sbase_kva
pu_to_kvar(q_pu, sbase_kva)   = q_pu * sbase_kva

# =======================================================================
# VERBATIM FROM Small_Network_PV.jl
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

get_sbase_kva(dm) = dm["settings"]["sbase"] * dm["settings"]["power_scale_factor"] / 1000

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

const THERMAL_RATINGS = Dict{Int,Vector{Float64}}()

function record_and_maybe_relax_thermal!(dm; relax::Bool)
    empty!(THERMAL_RATINGS)
    found = Set(lowercase(get(br, "name", "")) for br in values(dm["branch"])) ∩ TX_LEADS
    found == TX_LEADS || error("transformer leads not found in the math model: $(setdiff(TX_LEADS, found))")
    for (k, br) in dm["branch"]
        n = length(br["c_rating_a"])
        virtual = startswith(get(br, "name", ""), "_virtual")
        # assumed uniform rating (see THERMAL_A), in amps, for the post-solve check
        THERMAL_RATINGS[parse(Int, k)] = virtual ? fill(Inf, n) : fill(branch_rating_a(br), n)
        haskey(br, "rate_a") && (br["rate_a"] = fill(Inf, length(br["rate_a"])))
        if relax || virtual
            br["c_rating_a"] = fill(Inf, n)
        else   # enforced route: the same rating, in per unit, as a constraint
            ibase = get_sbase_kva(dm) / get(dm["bus"][string(br["f_bus"])], "vbase", 0.4 / sqrt(3))
            br["c_rating_a"] = fill(branch_rating_a(br) / ibase, n)
        end
    end
end

# --- PV control (verbatim) ---------------------------------------------
_sp(t::Real)   = t > 0 ? t + log1p(exp(-t)) : log1p(exp(t))
_dsp(t::Real)  = t >= 0 ? 1 / (1 + exp(-t)) : exp(t) / (1 + exp(t))
_d2sp(t::Real) = (s = _dsp(t); s * (1 - s))

function pwl_exact(U, pts)
    U <= pts[1][1]   && return pts[1][2]
    U >= pts[end][1] && return pts[end][2]
    for k in 1:length(pts)-1
        (x0, y0), (x1, y1) = pts[k], pts[k+1]
        U <= x1 && return y0 + (y1 - y0) * (U - x0) / (x1 - x0)
    end
end

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

const PV_VM = Dict{Int,Any}()

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

        p_head = mode == "vvw" ? smin(pav, S * curve_expr(U, vw_pu, reluV)) : min(pav, S)

        if PRIORITY == :var
            Qexpr = S * fvv
            Pexpr = smin(p_head, S * sqrt(1 - fvv^2))
        elseif PRIORITY == :watt
            Pexpr = p_head
            qh    = sqrt(S^2 - Pexpr^2 + eP^2)
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

# --- criteria (verbatim, except thermal is reported but optional) ------
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

# As Small_Network_PV.jl's thermal_worst, plus the largest PHASE current in A.
function thermal_worst(sbase_kva)
    worst, at = 0.0, ""
    in_max_a, iph_max_a, iph_at = 0.0, 0.0, ""
    for (l, br) in ref[:branch]
        startswith(get(br, "name", ""), "_virtual") && continue
        rating = get(THERMAL_RATINGS, l, Float64[])
        ibase = sbase_kva / get(ref[:bus][br["f_bus"]], "vbase", 0.4 / sqrt(3))   # A
        n = length(br["f_connections"])
        for c in 1:n
            I = hypot(JuMP.value(csr[c, l]), JuMP.value(csi[c, l]))
            c == 4 && (in_max_a = max(in_max_a, I * ibase))
            c <= 3 && I * ibase > iph_max_a && ((iph_max_a, iph_at) = (I * ibase, get(br, "name", string(l))))
            c <= length(rating) && isfinite(rating[c]) && rating[c] > 0 || continue
            u = I * ibase / rating[c]          # THERMAL_RATINGS are in amps (record_and_maybe_relax_thermal!)
            u > worst && ((worst, at) = (u, "$(get(br, "name", string(l))) c$c"))
        end
    end
    return (util=worst, at=at, in_max_a=in_max_a, iph_max_a=iph_max_a, iph_at=iph_at)
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
            iph_max_a=th.iph_max_a, iph_at=th.iph_at,
            loss_kw=network_losses(sbase_kva),
            src_kva=hypot(source_power_kw(sbase_kva)...),
            pv=pvd, statcom=sd)
end

function violations(m)
    v = Tuple{String,Float64}[]
    m.vmax_V > V_MAX_V + V_TOL_V     && push!(v, ("vmax "  * m.vmax_at, (m.vmax_V - V_MAX_V) / V_MAX_V))
    m.vmin_V < V_MIN_V - V_TOL_V     && push!(v, ("vmin "  * m.vmin_at, (V_MIN_V - m.vmin_V) / V_MIN_V))
    m.vuf    > VUF_LIMIT + VUF_TOL   && push!(v, ("vuf "   * m.vuf_at,  (m.vuf - VUF_LIMIT) / VUF_LIMIT))
    THERMAL_ON[] && m.therm > 1 + THERM_TOL && push!(v, ("therm " * m.therm_at, m.therm - 1))
    THERMAL_ON[] && m.src_kva > TX_KVA * (1 + THERM_TOL) &&
        push!(v, ("tx $(round(Int, TX_KVA)) kVA", m.src_kva / TX_KVA - 1))
    sort!(v, by = x -> -x[2])
    return v
end

function active_set(m)
    isnothing(m) && return ""
    a = String[]
    m.vmax_V >= V_MAX_V - 0.5          && push!(a, "vmax")
    m.vmin_V <= V_MIN_V + 0.5          && push!(a, "vmin")
    m.vuf    >= VUF_LIMIT - 5e-4       && push!(a, "vuf")
    THERMAL_ON[] && m.therm >= 0.98    && push!(a, "therm")
    THERMAL_ON[] && m.src_kva >= 0.98 * TX_KVA && push!(a, "tx")
    !isnothing(m.statcom) && m.statcom.worst_util >= 0.98 && push!(a, "device")
    return isempty(a) ? "none" : join(a, "+")
end

# =======================================================================
# NEW FOR LV12_248 -- network, loads, PV per customer
# =======================================================================

# Small CSV reader (no CSV.jl dependency). Files written by make_lv12_loads.py
# and by the OpenDSS reference step contain no quoted commas.
function read_csv_rows(path)
    lines = filter(l -> !isempty(strip(l)), readlines(path))
    hdr = String.(split(lines[1], ','))
    return [Dict(zip(hdr, String.(split(l, ',')))) for l in lines[2:end]]
end

const _PROFILES = Ref{Union{Nothing,Dict{String,Vector{Float64}}}}(nothing)
function profiles()
    if isnothing(_PROFILES[])
        _PROFILES[] = Dict(r["profile_id"] => [parse(Float64, r["h$(lpad(h, 2, '0'))"]) for h in 0:23]
                           for r in read_csv_rows(PROF_CSV))
    end
    return _PROFILES[]
end

load_name(load) = lowercase(get(load, "name", split(get(load, "source_id", ""), '.')[end]))

# Phase-domain source impedance from SRC_Z1/SRC_Z0, written over whatever
# PMD parsed (see the note at SRC_Z1). Keeps any neutral row/column as parsed.
# Returns the parsed (R11, X11) in ohm so stage 2 can report it.
const SRC_PARSED = Ref((NaN, NaN))
function set_source_impedance!(data_eng)
    vs = data_eng["voltage_source"]["source"]
    R = Matrix{Float64}(vs["rs"]); X = Matrix{Float64}(vs["xs"])
    size(R, 1) >= 3 || error("voltage source impedance is $(size(R)), expected at least 3x3")
    SRC_PARSED[] = (R[1, 1], X[1, 1])
    zs, zm = (SRC_Z0 + 2 * SRC_Z1) / 3, (SRC_Z0 - SRC_Z1) / 3
    for i in 1:3, j in 1:3
        R[i, j] = real(i == j ? zs : zm)
        X[i, j] = imag(i == j ? zs : zm)
    end
    vs["rs"] = R; vs["xs"] = X
    return SRC_PARSED[]
end

# Same structure as Small_Network_PV.jl's load_base_network, without the
# 4-bus lot bookkeeping: one MEN bond per customer bus, 10 ohm. The source
# impedance is set before anything else (set_source_impedance!).
function load_base_network(data_path; enforce_bounds=true, sbase_kva=100, neutral_ref=true)
    data_eng = PMD.parse_file(data_path, transformations=[PMD.transform_loops!])
    set_source_impedance!(data_eng)
    add_men!(data_eng; r_lot_ohm=MEN_OHM)

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

# Loads for any (draw, hour), from the SGSC files. kW -> pu as in
# Small_Network_PV.jl's set_phase_allocation!.
function set_loads!(dm, sbase_kva; draw, hour, pf=LOAD_PF)
    prof = profiles()
    asg = Dict(lowercase(r["load"]) => r["profile_id"]
               for r in read_csv_rows(ASSIGN_CSV) if parse(Int, r["draw"]) == draw)
    isempty(asg) && error("draw $draw not in $ASSIGN_CSV")
    tanphi = tan(acos(pf))
    n = 0
    for (_, load) in dm["load"]
        nm = load_name(load)
        haskey(asg, nm) || error("load '$nm' has no profile in draw $draw")
        kw = prof[asg[nm]][hour + 1]
        load["pd"] = [kw_to_pu(kw, sbase_kva)]
        load["qd"] = [kw_to_pu(kw * tanphi, sbase_kva)]
        n += 1
    end
    n == 121 || error("set $n loads, expected 121")
    return nothing
end

total_load_kw(dm, sbase_kva) = sum(pu_to_kw(sum(l["pd"]), sbase_kva) for l in values(dm["load"]))

# One single-phase PV unit per customer (bus and phase taken from the
# customer's load), rated s_kva. customers = nothing -> every customer;
# otherwise a collection of load names (count-based HC, later).
# Fields identical to Small_Network_PV.jl's add_pv_fleet! (type "PV_VVW").
function add_pv_customers!(dm, sbase_kva; s_kva, mode::Symbol, customers=nothing)
    ids = Int[]
    sel = isnothing(customers) ? nothing : Set(lowercase.(String.(collect(customers))))
    for (_, load) in sort(collect(dm["load"]), by = kv -> load_name(kv[2]))
        nm = load_name(load)
        isnothing(sel) || nm in sel || continue
        ph = filter(c -> c != 4, load["connections"])
        length(ph) == 1 || error("load $nm: expected single-phase, got $(load["connections"])")
        ph = ph[1]
        S   = kw_to_pu(s_kva, sbase_kva)
        pav = S * PAVAIL_FRAC

        gid = maximum(parse.(Int, collect(keys(dm["gen"])))) + 1
        g = deepcopy(dm["gen"]["1"])
        g["index"]         = gid
        g["gen_bus"]       = load["load_bus"]
        g["type"]          = "PV_VVW"
        g["name"]          = "pv_$(nm)"
        g["cost"]          = [0.0, 0.0]
        g["connections"]   = [ph, 4]
        g["configuration"] = PMD.WYE
        g["pmin"] = [0.0];       g["pmax"] = [2.0 * pav]
        g["qmin"] = [-1.2 * S];  g["qmax"] = [1.2 * S]
        g["pv_s"]      = S
        g["pv_pavail"] = pav
        g["pv_mode"]   = String(mode)
        g["pv_phase"]  = ph
        g["pv_lots"]   = 1.0
        dm["gen"][string(gid)] = g
        push!(ids, gid)
    end
    return ids
end

# =======================================================================
# BUILD + SOLVE (verbatim pipeline)
# =======================================================================
# HC-maximising device dispatch: minimise the highest phase-to-neutral |V|^2 on any
# real bus (epigraph variable t). The inverters' volt-var/volt-watt are fixed control
# laws in the constraints, so the device dispatch is the only freedom. Always feasible;
# a point passes iff the minimised maximum is within V_MAX_V + V_TOL_V (checked post hoc
# with the same violations() test as every other route), and nested modes give
# HC(D) >= HC(B), HC(D) >= HC(C) up to solver tolerance.
function add_minvmax_objective!()
    t = JuMP.@variable(model, base_name = "vpn2_max", lower_bound = 0.0, start = 1.2)
    for (i, bus) in ref[:bus]
        is_real_bus(bus) || continue
        terms = bus["terminals"]
        4 in terms || continue
        for p in 1:3
            p in terms || continue
            JuMP.@constraint(model, (vr[p, i] - vr[4, i])^2 + (vi[p, i] - vi[4, i])^2 <= t)
        end
    end
    JuMP.@objective(model, Min, t)
    return t
end

# LOCAL STATCOM control (realistic tier): each device minimises the negative-sequence
# voltage at its OWN bus (its point of connection), plus a tiny penalty on its output so
# that it does the least needed when the bus is already balanced. Uses only the device's
# terminal voltages; the OPF finds the steady state a local feedback controller would
# settle at. V2 = (Va + a^2 Vb + a Vc)/3 is the same phase-to-ground or phase-to-neutral
# (a neutral shift is common to all phases, i.e. zero sequence). With two devices the
# objective is the sum of their own terms (they sit on different circuits).
#
# Further LOCAL objectives (same rules: only quantities measurable at the device's own
# connection, no network limits imposed -- hc_point never sends a "local_*" case to the
# bounds-enforced route, so pass/fail is judged post hoc on the operating point the local
# controller produces):
#   local_vuf      |V2|^2                                  negative sequence only (as above)
#   local_vbal     |V0pn|^2 + |V2|^2                      phase-to-neutral voltages driven towards
#                  a balanced set, i.e. the three V_pn magnitudes equalised. V0pn is the zero
#                  sequence of the phase-to-NEUTRAL voltages (= V0pg - Vn), the part VUF ignores.
#   local_vbal_vv  local_vbal + positive-sequence voltage regulation with a deadband:
#                  LOCAL_VV_W * s^2, s >= |V1|^2 - Vset^2, s >= 0  (absorbs Q only when the
#                  device's own V1 exceeds Vset; default Vset = 240 V, the AS/NZS 4777.2
#                  volt-var knee, so the device and the inverters start absorbing together).
# V1 and V2 are the same phase-to-ground or phase-to-neutral (a neutral shift is pure zero
# sequence). With two devices the objective is the SUM of their own terms; they sit on
# different circuits (coupling ~0.01 V, see siting screen), so this is not a hidden
# coordination between them.
const LOCAL_EFFORT_W = 1e-6     # tie-breaker only: full output (~0.03 pu^2) costs ~3e-8, vs |V2|^2 ~1e-5 pu^2 at an unbalanced PCC
const LOCAL_SCALE    = 1e4      # objective scaling for Ipopt
const LOCAL_VSET_V = parse(Float64, get(ENV, "LV12_LOCAL_VSET", "240"))
const LOCAL_VV_W   = parse(Float64, get(ENV, "LV12_LOCAL_VV_W", "1.0"))
const LOCAL_OBJECTIVES = ("local_vuf", "local_vbal", "local_vbal_vv")
is_local(obj) = obj in LOCAL_OBJECTIVES

function add_local_objective!(kind::String)
    sc = [(id, g) for (id, g) in ref[:gen] if get(g, "type", "") == "STATCOM"]
    isempty(sc) && return nothing
    c1, s1 = cos(2pi / 3), sin(2pi / 3)      # a   = 1 at 120 deg
    c2, s2 = cos(4pi / 3), sin(4pi / 3)      # a^2 = 1 at 240 deg
    # sequence components of three phasors given as (re, im) expressions
    seq2(xr, xi) = ((xr[1] + c2 * xr[2] - s2 * xi[2] + c1 * xr[3] - s1 * xi[3]) / 3,
                    (xi[1] + s2 * xr[2] + c2 * xi[2] + s1 * xr[3] + c1 * xi[3]) / 3)
    seq1(xr, xi) = ((xr[1] + c1 * xr[2] - s1 * xi[2] + c2 * xr[3] - s2 * xi[3]) / 3,
                    (xi[1] + s1 * xr[2] + c1 * xi[2] + s2 * xr[3] + c2 * xi[3]) / 3)
    seq0(xr, xi) = ((xr[1] + xr[2] + xr[3]) / 3, (xi[1] + xi[2] + xi[3]) / 3)
    terms  = JuMP.QuadExpr[]
    effort = JuMP.QuadExpr[]
    vset2  = (LOCAL_VSET_V / (1000 * ref[:bus][sc[1][2]["gen_bus"]]["vbase"]))^2
    for (id, g) in sc
        i = g["gen_bus"]
        vpr = [vr[p, i] for p in 1:3]; vpi = [vi[p, i] for p in 1:3]
        if kind == "local_vuf"
            a, b = seq2(vpr, vpi); push!(terms, a^2 + b^2)
        elseif kind in ("local_vbal", "local_vbal_vv")
            a, b = seq2(vpr, vpi)
            zr, zi = seq0([vr[p, i] - vr[4, i] for p in 1:3], [vi[p, i] - vi[4, i] for p in 1:3])
            push!(terms, a^2 + b^2 + zr^2 + zi^2)
            if kind == "local_vbal_vv"
                u1r, u1i = seq1(vpr, vpi)
                s = JuMP.@variable(model, lower_bound = 0.0, base_name = "vv_excess_$id")
                JuMP.@constraint(model, u1r^2 + u1i^2 - vset2 <= s)
                push!(terms, LOCAL_VV_W * s^2)
            end
        else
            error("unknown local objective $kind")
        end
        push!(effort, sum(pg[c, id]^2 + qg[c, id]^2 for c in 1:3))
    end
    JuMP.@objective(model, Min, LOCAL_SCALE * (sum(terms) + LOCAL_EFFORT_W * sum(effort)))
    return nothing
end
add_local_vuf_objective!() = add_local_objective!("local_vuf")

const TIMING = Dict{Symbol,Float64}(:parse => 0.0, :build => 0.0, :optimize => 0.0,
                                    :ipopt => 0.0, :n => 0.0)

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
        _core_variables!()
        _core_constraints!()
    else
        include(joinpath(CORE_DIR, "variables.jl"))
        include(joinpath(CORE_DIR, "constraints_unified.jl"))
    end

    add_pv_control!()
    bounds_active && add_neutral_referenced_bounds!()

    # "minvmax" is defined here, not in the core: build the core's "cost" objective,
    # then replace it (every core objective is a plain JuMP.@objective, so this overrides).
    global objective = (objective_mode in ("minvmax", "none") || is_local(objective_mode)) ? "cost" : objective_mode
    global objective_target_gens = [i for (i, g) in ref[:gen] if get(g, "type", "") == "STATCOM"]
    if USE_COMPILED_CORE
        _core_objectives!()
    else
        include(joinpath(CORE_DIR, "objectives_FIXED.jl"))
    end
    objective_mode == "minvmax"   && add_minvmax_objective!()
    is_local(objective_mode)      && add_local_objective!(objective_mode)

    TIMING[:build] += time() - t0
    t1 = time()
    JuMP.optimize!(model)
    TIMING[:optimize] += time() - t1
    TIMING[:ipopt]    += JuMP.solve_time(model)
    TIMING[:n]        += 1
    return JuMP.termination_status(model)
end

const _BASE_CACHE = Dict{Bool,Any}()
function base_network(bounds_active::Bool)
    if !haskey(_BASE_CACHE, bounds_active)
        t = time()
        _BASE_CACHE[bounds_active] = load_base_network(NET_FILE; enforce_bounds=bounds_active,
                                                       neutral_ref=true)
        TIMING[:parse] += time() - t
    end
    t = time()
    dm = deepcopy(_BASE_CACHE[bounds_active])
    TIMING[:parse] += time() - t
    return dm
end

# One case. pv_kva = uniform size per customer (0 = no PV).
function run_case(; pv_kva, pv_mode::Symbol, pv_customers=nothing, statcom_mode::Symbol=:none,
                  statcom_kva=STATCOM_KVA, statcom_bus=STATCOM_BUS, objective="cost",
                  bounds_active=false, draw=LOAD_DRAW, hour=HOUR)
    dm = base_network(bounds_active)
    record_and_maybe_relax_thermal!(dm; relax = !THERMAL_ON[] || !bounds_active)
    (draw, hour) == (DSS_DRAW, DSS_HOUR) || set_loads!(dm, SBASE_KVA; draw=draw, hour=hour)
    pv_kva > 0 && add_pv_customers!(dm, SBASE_KVA; s_kva=pv_kva, mode=pv_mode, customers=pv_customers)

    # statcom_bus: one bus name, or a vector of names (one device of statcom_kva at each)
    sc_buses = statcom_bus isa AbstractString ? [statcom_bus] : collect(statcom_bus)
    gids = statcom_mode == :none ? Int[] :
        reduce(vcat, [add_statcom!(dm, SBASE_KVA; bus_name=b, rating_kva_total=statcom_kva,
                                   mode=statcom_mode) for b in sc_buses])

    st = try
        build_and_solve_pv(dm; objective_mode=objective, bounds_active=bounds_active,
                           vuf_cap = bounds_active ? VUF_LIMIT : Inf)
    catch e
        println("      ERROR: ", sprint(showerror, e)); :ERROR
    end
    return st, gids
end

# =======================================================================
# HOSTING CAPACITY -- uniform size per customer (logic as Small_Network_PV.jl)
# =======================================================================
const POINT_HEADER = ["stage", "pv_mode", "statcom", "x_kva_per_customer", "route", "status", "pass",
    "binding", "vmax_v", "vmax_at", "vmin_v", "vuf_pct", "vuf_at", "nev_pu", "iph_max_a", "iph_at",
    "in_max_a", "loss_kw", "pv_p_kw", "pv_pavail_kw", "pv_curt_kw", "pv_q_kvar", "droop_err_frac",
    "sc_p1_kw", "sc_p2_kw", "sc_p3_kw", "sc_q1_kvar", "sc_q2_kvar", "sc_q3_kvar", "sc_util"]
point_rows = Vector{Any}[]

function point_row(stage, pv_mode, sc, x, route, st, pass, binding, m)
    head = Any[stage, String(pv_mode), sc, x, route, string(st), pass, binding]
    isnothing(m) && return Any[head; fill(nothing, length(POINT_HEADER) - length(head))]
    s = m.statcom
    sv(f, k) = isnothing(s) ? nothing : f(s)[k]
    return Any[head...,
        m.vmax_V, m.vmax_at, m.vmin_V, 100 * m.vuf, m.vuf_at, m.nev, m.iph_max_a, m.iph_at,
        m.in_max_a, m.loss_kw, m.pv.p_kw, m.pv.pavail_kw, m.pv.curt_kw, m.pv.q_kvar, m.pv.droop_err,
        sv(x -> x.net_p, 1), sv(x -> x.net_p, 2), sv(x -> x.net_p, 3),
        sv(x -> x.net_q, 1), sv(x -> x.net_q, 2), sv(x -> x.net_q, 3),
        isnothing(s) ? nothing : s.worst_util]
end

function hc_point(x; pv_mode, statcom_mode=:none, stage="", enforce_fallback=true,
                  statcom_kva=STATCOM_KVA, statcom_bus=STATCOM_BUS, pv_customers=nothing,
                  objective=STATCOM_HC_OBJECTIVE)
    sc = statcom_mode == :none ? "none" : String(statcom_mode)
    relaxed_only = statcom_mode == :none || !enforce_fallback || is_local(objective)

    st, gids = run_case(pv_kva=x, pv_mode=pv_mode, pv_customers=pv_customers,
                        statcom_mode=statcom_mode,
                        statcom_kva=statcom_kva, statcom_bus=statcom_bus,
                        objective=objective, bounds_active=false)
    if solved(st)
        m = evaluate_point(gids, SBASE_KVA, statcom_kva)
        v = violations(m)
        pass = isempty(v)
        push!(point_rows, point_row(stage, pv_mode, sc, x, "relaxed", st, pass, pass ? "" : v[1][1], m))
        # Under "minvmax" the device already minimised the worst voltage, so a vmax failure
        # cannot be rescued by the bounds-enforced route; only a VUF failure can.
        skip = objective == "minvmax" && !pass && startswith(v[1][1], "vmax")
        (pass || relaxed_only || skip) && return (pass=pass, binding=pass ? "" : v[1][1], status=string(st), m=m)
        binding = v[1][1]
    else
        push!(point_rows, point_row(stage, pv_mode, sc, x, "relaxed", st, false, "not solved", nothing))
        binding = "relaxed solve: $st"
        relaxed_only && return (pass=false, binding=binding, status=string(st), m=nothing)
    end

    st2, gids2 = run_case(pv_kva=x, pv_mode=pv_mode, pv_customers=pv_customers,
                          statcom_mode=statcom_mode,
                          statcom_kva=statcom_kva, statcom_bus=statcom_bus,
                          objective=objective, bounds_active=true)
    if solved(st2)
        m2 = evaluate_point(gids2, SBASE_KVA, statcom_kva)
        v2 = violations(m2)
        pass2 = isempty(v2)
        push!(point_rows, point_row(stage, pv_mode, sc, x, "enforced", st2, pass2, pass2 ? "" : v2[1][1], m2))
        return (pass=pass2, binding=pass2 ? "" : v2[1][1], status=string(st2), m=m2)
    end
    push!(point_rows, point_row(stage, pv_mode, sc, x, "enforced", st2, false, "infeasible: $binding", nothing))
    return (pass=false, binding=binding, status=string(st2), m=nothing)
end

function hc_bisect(; pv_mode, statcom_mode=:none, stage="", lo_start=HC_LO,
                   enforce_fallback=true, statcom_kva=STATCOM_KVA, statcom_bus=STATCOM_BUS)
    n = 0
    pt(x) = (n += 1; hc_point(x; pv_mode, statcom_mode, stage, enforce_fallback, statcom_kva, statcom_bus))

    r_lo = pt(lo_start)
    if !r_lo.pass && lo_start > HC_LO
        println("      note: lo_start = $(round(lo_start, digits=3)) failed ($(r_lo.binding)); restarting from $(HC_LO)")
        lo_start = HC_LO
        r_lo = pt(lo_start)
    end
    r_lo.pass || return (hc=NaN, binding="fails at x=$(lo_start): $(r_lo.binding)", n=n,
                         capped=false, active="", m=nothing)

    lo, m_lo = lo_start, r_lo.m
    hi = max(HC_HI, lo_start + 1.0)
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

fmt_hc(r) = isnan(r.hc) ? "   --  " : (r.capped ? @sprintf(">%6.3f", r.hc) : @sprintf("%7.3f", r.hc))

# =======================================================================
# COMPILE-ONCE LOADING OF THE CORE FILES (verbatim)
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
    end
    return names
end

function _collect_assigned!(names::Set{Symbol}, ex)
    ex isa Expr || return names
    h = ex.head
    if h in (:(=), Symbol("+="), Symbol("-="), Symbol("*="), Symbol("/="))
        _is_funcdef(ex) || _lhs_names!(names, ex.args[1])
        _collect_assigned!(names, ex.args[2])
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
            Core.eval(Main, st)
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

function validate_compiled_core(; x=3.0, tol=1e-7)
    global USE_COMPILED_CORE
    keep = USE_COMPILED_CORE
    try
        for obj in ("cost", "VUF")
            snap = Dict{Bool,Any}()
            for compiled in (false, true)
                USE_COMPILED_CORE = compiled
                st, _ = run_case(pv_kva=x, pv_mode=:vvw, statcom_mode=:pq, objective=obj)
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
# STAGE CHECKS
# =======================================================================
const STAGE_ROWS = Vector{Any}[]

function check(stage, what, value, target, tol)
    pass = isfinite(value) && abs(value - target) <= tol
    @printf("    %-40s %11.3f   target %11.3f ± %-6.3f %s\n", what, value, target, tol, pass ? "PASS" : "FAIL")
    push!(STAGE_ROWS, Any[stage, what, value, target, tol, pass])
    return pass
end

function check_bool(stage, what, ok, detail="")
    @printf("    %-40s %-42s %s\n", what, detail, ok ? "PASS" : "FAIL")
    push!(STAGE_ROWS, Any[stage, what, detail, nothing, nothing, ok])
    return ok
end

function stage_gate(stage, ok)
    ok && return
    msg = "stage $stage failed -- fix this before trusting anything after it"
    write_csv(joinpath(OUTDIR, "lv12_stage_checks.csv"),
              ["stage", "check", "value", "target", "tolerance", "pass"], STAGE_ROWS)
    STOP_ON_FAIL ? error(msg) : println("  !! ", msg)
end

function source_power_kw(sbase_kva)
    P = Q = 0.0
    for (i, g) in ref[:gen]
        ref[:bus][g["gen_bus"]]["bus_type"] == 3 || continue
        for c in filter(c -> c != 4, g["connections"])
            P += JuMP.value(pg[c, i]); Q += JuMP.value(qg[c, i])
        end
    end
    return pu_to_kw(P, sbase_kva), pu_to_kvar(Q, sbase_kva)
end

# Every bus's phase-to-neutral voltages and NEV against the OpenDSS reference.
function compare_to_reference(case::String)
    rows = [r for r in read_csv_rows(REF_CSV) if r["case"] == case]
    byname = Dict(lowercase(get(b, "name", "")) => i for (i, b) in ref[:bus])
    dv_max, dv_at, dnev_max, unmatched = 0.0, "", 0.0, String[]
    for r in rows
        nm = lowercase(r["bus"])
        haskey(byname, nm) || (push!(unmatched, nm); continue)
        i = byname[nm]; bus = ref[:bus][i]
        vb = get(bus, "vbase", 0.4 / sqrt(3)) * 1000
        t = bus["terminals"]
        vn = 4 in t ? JuMP.value(vr[4, i]) + im * JuMP.value(vi[4, i]) : 0.0im
        for (p, col) in zip(1:3, ("Van_V", "Vbn_V", "Vcn_V"))
            v = abs(JuMP.value(vr[p, i]) + im * JuMP.value(vi[p, i]) - vn) * vb
            d = abs(v - parse(Float64, r[col]))
            d > dv_max && ((dv_max, dv_at) = (d, "$nm ph$p"))
        end
        dnev_max = max(dnev_max, abs(abs(vn) * vb - parse(Float64, r["Vng_V"])))
    end
    return (dv_max=dv_max, dv_at=dv_at, dnev_max=dnev_max, n=length(rows), unmatched=unmatched)
end

function check_power_flow(stage, case, T)
    P, Q = source_power_kw(SBASE_KVA)
    m = evaluate_point(Int[], SBASE_KVA, STATCOM_KVA)
    c = compare_to_reference(case)
    vb_V = 400 / sqrt(3)
    ok = true
    ok &= check_bool(stage, "buses matched to OpenDSS", isempty(c.unmatched),
                     "$(c.n - length(c.unmatched))/$(c.n)" * (isempty(c.unmatched) ? "" : "  e.g. $(first(c.unmatched))"))
    ok &= check(stage, "source P (kW)",            P, T.src_kw,   TOL_P)
    ok &= check(stage, "source Q (kvar)",          Q, T.src_kvar, TOL_P)
    ok &= check(stage, "V_pn max (V)",      m.vmax_V, T.vmax,     TOL_V)
    ok &= check(stage, "V_pn min (V)",      m.vmin_V, T.vmin,     TOL_V)
    ok &= check(stage, "max NEV (V)",  m.nev * vb_V, T.nev,      TOL_NEV)
    ok &= check(stage, "max VUF (%)",    100 * m.vuf, T.vuf,      TOL_VUF)
    ok &= check(stage, "max |ΔV_pn| any bus (V)  [$(c.dv_at)]", c.dv_max, 0.0, TOL_V)
    ok &= check(stage, "max |ΔNEV| any bus (V)", c.dnev_max, 0.0, TOL_NEV)
    @printf("    (info) losses %.3f kW; largest phase current %.1f A on %s; neutral %.1f A\n",
            m.loss_kw, m.iph_max_a, m.iph_at, m.in_max_a)
    return ok
end

# =======================================================================
# RUN
# =======================================================================
USE_COMPILED_CORE && load_compiled_core!()

_probe = load_base_network(NET_FILE)
const SBASE_KVA = get_sbase_kva(_probe)
_src_vm = [b["vm"][1:3] for (_, b) in _probe["bus"] if b["bus_type"] == 3]

println("=" ^ 96)
println(" LV12_248 OPF -- staged port check, then HC0 / HC1 / HC2")
println(" sbase = $SBASE_KVA kVA | loads: draw $LOAD_DRAW, $(HOUR):00, pf $LOAD_PF | source vm (pu) = $(_src_vm) | source model: $SOURCE")
println(" criteria: $(V_MIN_V)-$(V_MAX_V) V phase-to-neutral (+$(V_TOL_V) V), VUF $(100*VUF_LIMIT)%",
        THERMAL_ON[] ? ", thermal $(round(Int, THERMAL_A)) A per conductor (assumed), transformer leads $(round(Int, TX_LEAD_A)) A, transformer $(round(Int, TX_KVA)) kVA" : ", thermal NOT checked")
println(" STATCOM: $(STATCOM_KVA) kVA at $STATCOM_BUS | core files: ",
        USE_COMPILED_CORE ? "compiled once" : "include() per solve")
println("=" ^ 96)

if USE_COMPILED_CORE && VALIDATE_CORE
    println("\n Validating compiled core against include() (slow: the include() side)...")
    validate_compiled_core()
end
foreach(k -> TIMING[k] = 0.0, collect(keys(TIMING)))

# --- STAGE 1: parse ----------------------------------------------------
if 1 in STAGES
    println("\n STAGE 1 -- parse $(basename(NET_FILE))")
    eng = PMD.parse_file(NET_FILE, transformations=[PMD.transform_loops!])
    ok = true
    nb  = length(eng["bus"])
    nl  = length(get(eng, "line", Dict()))
    nsw = length(get(eng, "switch", Dict()))
    ok &= check_bool(1, "buses", nb == 248, "$nb (expect 248)")
    ok &= check_bool(1, "lines incl. 7 ex-switches", nl == 247, "$nl (expect 247)")
    ok &= check_bool(1, "no switch objects", nsw == 0, "$nsw")
    nc = sort(unique(length(l["f_connections"]) for l in values(eng["line"])))
    ok &= check_bool(1, "every line four-wire", nc == [4], "conductors per line: $nc")
    maxterm = maximum(maximum(b["terminals"]) for b in values(eng["bus"]))
    ok &= check_bool(1, "no fifth conductor", maxterm <= 4, "highest terminal id $maxterm")
    nld = length(eng["load"])
    okc = all(length(l["connections"]) == 2 && l["connections"][2] == 4 for l in values(eng["load"]))
    ok &= check_bool(1, "121 loads, phase-to-neutral", nld == 121 && okc, "$nld loads")
    ok &= check(1, "total load P (kW)",   sum(sum(l["pd_nom"]) for l in values(eng["load"])), T_LOAD_KW, 0.01)
    ok &= check(1, "total load Q (kvar)", sum(sum(l["qd_nom"]) for l in values(eng["load"])), T_LOAD_KVAR, 0.01)
    src = get(eng["bus"], "b908", nothing)
    gr = isnothing(src) ? Int[] : get(src, "grounded", Int[])
    ok &= check_bool(1, "source neutral grounded at b908", 4 in gr, "grounded = $gr")
    stage_gate(1, ok)

    # --- STAGE 2: earthing ---------------------------------------------
    if 2 in STAGES
        println("\n STAGE 2 -- MEN bonds, source impedance")
        eng2 = deepcopy(eng)
        men = sort(lowercase.(add_men!(eng2; r_lot_ohm=MEN_OHM)))
        loadb = sort(unique(lowercase(l["bus"]) for l in values(eng["load"])))
        csvb = sort([lowercase(r["bus"]) for r in read_csv_rows(MEN_CSV)])
        ok2 = check_bool(2, "MEN bonds", length(men) == 121, "$(length(men)) bonds at $(MEN_OHM) Ω")
        ok2 &= check_bool(2, "MEN buses = load buses = men_buses.csv", men == loadb == csvb, "")
        ok2 &= check_bool(2, "math model loads", length(_probe["load"]) == 121, "$(length(_probe["load"]))")
        ok2 &= check(2, "math model total load (kW)", total_load_kw(_probe, SBASE_KVA), T_LOAD_KW, 0.01)
        # set_loads! self-test: re-applying draw 45 / 13:00 from the CSVs must not change anything
        _t = deepcopy(_probe); set_loads!(_t, SBASE_KVA; draw=DSS_DRAW, hour=DSS_HOUR)
        ok2 &= check(2, "set_loads!(draw 45, 13:00) total (kW)", total_load_kw(_t, SBASE_KVA), T_LOAD_KW, 0.01)
        # source impedance as it reached the math model (the virtual source branch), in ohm
        rb = parse(Int, first(k for (k, b) in _probe["bus"] if b["bus_type"] == 3))
        vbr = [br for br in values(_probe["branch"]) if rb in (br["f_bus"], br["t_bus"])]
        zs_set = (SRC_Z0 + 2 * SRC_Z1) / 3
        okz, detail = length(vbr) == 1, "$(length(vbr)) branches at the source bus"
        if okz
            zb = (_probe["bus"][string(rb)]["vbase"] * 1e3)^2 / (SBASE_KVA * 1e3)
            Rm, Xm = vbr[1]["br_r"][1, 1] * zb, vbr[1]["br_x"][1, 1] * zb
            okz = isapprox(Rm, real(zs_set); rtol=1e-3) && isapprox(Xm, imag(zs_set); rtol=1e-3)
            pr, px = SRC_PARSED[]
            detail = @sprintf("self %.4f + j%.4f mOhm in the model (set %.4f + j%.4f; PMD had parsed %.4g + j%.4g Ohm)",
                              1e3 * Rm, 1e3 * Xm, 1e3 * real(zs_set), 1e3 * imag(zs_set), pr, px)
        end
        ok2 &= check_bool(2, "source impedance", okz, detail)
        stage_gate(2, ok2)
    end
end

# --- STAGE 3: power flow, no PV ----------------------------------------
if 3 in STAGES && PORT_CHECK
    println("\n STAGE 3 -- power flow, no PV, no device (draw 45, 13:00)")
    st, _ = run_case(pv_kva=0.0, pv_mode=:unity, draw=DSS_DRAW, hour=DSS_HOUR)
    ok = check_bool(3, "solve", solved(st), string(st))
    ok && (ok &= check_power_flow(3, "no_pv", T_NOPV))
    stage_gate(3, ok)
end

# --- STAGE 4: 3.0 kW PV per customer, unity PF ---------------------------
if 4 in STAGES && PORT_CHECK
    println("\n STAGE 4 -- 121 x 3.0 kW PV, unity PF")
    st, _ = run_case(pv_kva=3.0, pv_mode=:unity, draw=DSS_DRAW, hour=DSS_HOUR)
    ok = check_bool(4, "solve", solved(st), string(st))
    ok && (ok &= check_power_flow(4, "pv_3kW_unity", T_PV3))
    stage_gate(4, ok)
end

HC = Dict{Symbol,Any}()
hc_rows = Vector{Any}[]

# --- STAGE 5: HC0 --------------------------------------------------------
if 5 in STAGES
    println("\n STAGE 5 -- HC0: uniform PV per customer, unity PF, no device")
    th_saved = THERMAL_ON[]
    THERMAL_ON[] = false            # the OpenDSS comparison is voltage-only
    r = hc_bisect(pv_mode=:unity, stage="5v")
    THERMAL_ON[] = th_saved
    @printf("    HC0 (voltage + VUF only) = %s kVA/customer   first fail: %s  (%d solves)\n",
            fmt_hc(r), r.binding, r.n)
    PORT_CHECK && stage_gate(5, check(5, "HC0 voltage-only (kVA/customer)", r.hc, T_HC0, TOL_HC))
    if THERMAL_ON[]
        r = hc_bisect(pv_mode=:unity, stage="5")
        @printf("    HC0 (with thermal + transformer) = %s kVA/customer   first fail: %s  (%d solves)\n",
                fmt_hc(r), r.binding, r.n)
    end
    HC[:hc0] = r
    push!(hc_rows, Any["5", "unity", "none", "A", r.hc, r.capped, r.active, r.binding, r.n, nothing])
end

# --- STAGE 6: HC1 --------------------------------------------------------
if 6 in STAGES
    println("\n STAGE 6 -- HC1: AS/NZS 4777.2 volt-var + volt-watt, no device")
    r = hc_bisect(pv_mode=HC1_MODE, stage="6")
    HC[:hc1] = r
    push!(hc_rows, Any["6", String(HC1_MODE), "none", "A", r.hc, r.capped, r.active, r.binding, r.n, nothing])
    @printf("    HC1 = %s kVA/customer   active: %-12s first fail: %s  (%d solves)\n",
            fmt_hc(r), r.active, r.binding, r.n)
    if !isnothing(r.m)
        @printf("    at HC1: Vmax %.2f V (%s), VUF %.3f%%, PV %.1f kW of %.1f available (curtailed %.1f), Q %.1f kvar\n",
                r.m.vmax_V, r.m.vmax_at, 100 * r.m.vuf, r.m.pv.p_kw, r.m.pv.pavail_kw, r.m.pv.curt_kw, r.m.pv.q_kvar)
        check(6, "droop error at HC1 (fraction of S)", r.m.pv.droop_err, 0.0, 0.005)
    end
end

# --- STAGE 7: HC2 per STATCOM mode ----------------------------------------
if 7 in STAGES && haskey(HC, :hc1) && !isnan(HC[:hc1].hc)
    hc1 = HC[:hc1].hc
    println("\n STAGE 7 -- HC2 on top of HC1 = $(round(hc1, digits=3)), $(STATCOM_KVA) kVA at $STATCOM_BUS")
    for (smode, code, label) in STATCOM_MODES
        r7 = hc_bisect(pv_mode=HC1_MODE, statcom_mode=smode, stage="7", lo_start=hc1)
        HC[smode] = r7
        push!(hc_rows, Any["7", String(HC1_MODE), String(smode), code, r7.hc, r7.capped, r7.active,
                           r7.binding, r7.n, r7.hc - hc1])
        util7 = (isnothing(r7.m) || isnothing(r7.m.statcom)) ? NaN : r7.m.statcom.worst_util
        @printf("    %-2s %-18s HC2 = %s   dHC = %+6.3f   device util %4.2f   active: %-18s (%2d)\n",
                code, label, fmt_hc(r7), r7.hc - hc1, util7, r7.active, r7.n)
    end
end

# --- STAGE 8: siting sweep ------------------------------------------------
# b733 sits on a different branch from the HC1 binding bus (b3278); the two
# share only the first 2 m from the transformer, so a device there cannot move
# b3278 (OpenDSS: -0.009 V per 30 kVA). Sweep the weak branch instead.
const SITING_SITES = ["b733", "b2351", "b1142", "b2593", "b3278", "b740", "b2848"]
const SITING_MODES = [(:qonly, "B"), (:ponly, "C"), (:pq, "D")]   # Bc cannot lower voltage
const SITING_KVA   = [30.0]                                       # add 60.0, 100.0 for rating curves
siting_rows = Vector{Any}[]
if 8 in STAGES && haskey(HC, :hc1) && !isnan(HC[:hc1].hc)
    hc1 = HC[:hc1].hc
    println("\n STAGE 8 -- siting sweep on top of HC1 = $(round(hc1, digits=3))")
    for kva in SITING_KVA, site in SITING_SITES, (smode, code) in SITING_MODES
        r8 = hc_bisect(pv_mode=HC1_MODE, statcom_mode=smode, stage="8", lo_start=hc1,
                       statcom_kva=kva, statcom_bus=site)
        util8 = (isnothing(r8.m) || isnothing(r8.m.statcom)) ? NaN : r8.m.statcom.worst_util
        push!(siting_rows, Any[site, kva, code, r8.hc, r8.hc - hc1, r8.capped, r8.active,
                               r8.binding, util8, r8.n])
        @printf("    %-6s %5.0f kVA  %-2s  HC2 = %s   dHC = %+6.3f   util %4.2f   next limit: %-22s (%2d)\n",
                site, kva, code, fmt_hc(r8), r8.hc - hc1, util8, r8.binding, r8.n)
    end
end

# --- STAGE 9: one device vs two (one per branch) --------------------------
# Stage 8's ceiling (4.788) is set by b3162 ph B on the OTHER branch (b733's,
# shared from b2620). Two 30 kVA units, one per branch, against one 60 kVA unit
# (same total rating) at the best single site.
const PAIR_CASES = [(["b3278"], 60.0), (["b3278", "b3162"], 30.0),
                    (["b3278", "b2620"], 30.0), (["b1142", "b3162"], 30.0)]
const PAIR_MODES = [(:qonly, "B"), (:pq, "D")]
pair_rows = Vector{Any}[]
if 9 in STAGES && haskey(HC, :hc1) && !isnan(HC[:hc1].hc)
    hc1 = HC[:hc1].hc
    println("\n STAGE 9 -- one device vs two, on top of HC1 = $(round(hc1, digits=3))")
    for (sites, kva) in PAIR_CASES, (smode, code) in PAIR_MODES
        r9 = hc_bisect(pv_mode=HC1_MODE, statcom_mode=smode, stage="9", lo_start=hc1,
                       statcom_kva=kva, statcom_bus=sites)
        util9 = (isnothing(r9.m) || isnothing(r9.m.statcom)) ? NaN : r9.m.statcom.worst_util
        label = join(sites, "+")
        push!(pair_rows, Any[label, length(sites), kva, length(sites) * kva, code, r9.hc,
                             r9.hc - hc1, r9.capped, r9.active, r9.binding, util9, r9.n])
        @printf("    %-12s %d x %4.0f kVA  %-2s  HC2 = %s   dHC = %+6.3f   util %4.2f   next limit: %-22s (%2d)\n",
                label, length(sites), kva, code, fmt_hc(r9), r9.hc - hc1, util9, r9.binding, r9.n)
    end
end

# --- STAGE 10: STATCOM on the unity-PF baseline (no volt-var) -------------
# The device on top of HC0 instead of HC1: legacy inverters, or no volt-var mandate.
const BASE0_CASES = [(["b2848"], 30.0), (["b3278"], 30.0), (["b3278", "b3162"], 30.0)]
const BASE0_MODES = [(:qonly, "B"), (:ponly, "C"), (:pq, "D")]
base0_rows = Vector{Any}[]
if 10 in STAGES && haskey(HC, :hc0) && !isnan(HC[:hc0].hc)
    hc0 = HC[:hc0].hc
    println("\n STAGE 10 -- STATCOM on the unity-PF baseline, on top of HC0 = $(round(hc0, digits=3))")
    for (sites, kva) in BASE0_CASES, (smode, code) in BASE0_MODES
        r10 = hc_bisect(pv_mode=:unity, statcom_mode=smode, stage="10", lo_start=hc0,
                        statcom_kva=kva, statcom_bus=sites)
        util10 = (isnothing(r10.m) || isnothing(r10.m.statcom)) ? NaN : r10.m.statcom.worst_util
        label = join(sites, "+")
        push!(base0_rows, Any[label, length(sites), kva, code, r10.hc, r10.hc - hc0, r10.capped,
                              r10.active, r10.binding, util10, r10.n])
        @printf("    %-12s %d x %4.0f kVA  %-2s  HC = %s   dHC vs HC0 = %+6.3f   util %4.2f   next limit: %-22s (%2d)\n",
                label, length(sites), kva, code, fmt_hc(r10), r10.hc - hc0, util10, r10.binding, r10.n)
    end
end

# --- STAGE 11: sequential random adoption from an existing fleet ----------
# HC as a COUNT (claude/lv12_hc_procedure.md). ADOPT_PEN of the customers already
# have a 5 kVA inverter (chosen at random, one seed per "init" draw); the remaining
# customers apply one at a time in a random order (one seed per "new" draw). HC_new =
# how many applicants connect before the first limit breach. The same orders are used
# for every configuration, so differences are paired order by order. The nested design
# (fleets x applicant orders) separates the effect of the initial fleet from the
# effect of who applies next.
const ADOPT_PEN    = parse(Float64, get(ENV, "LV12_PEN", "0.55"))   # existing penetration (0.55 base; e.g. 0.80 new estates)
const ADOPT_KVA    = 5.0                                             # every system (base case)
const ADOPT_N_INIT = parse(Int, get(ENV, "LV12_NINIT", "10"))        # initial-fleet draws
const ADOPT_N_NEW  = parse(Int, get(ENV, "LV12_NNEW", "10"))         # applicant orders per fleet
# :sequential = test applicants one at a time in order and stop at the first breach (the
# definition; a DNSP refuses that applicant). :bisect is faster but assumes pass/fail is
# monotone in k, which is FALSE here (13 Oct quick test: 7 crossings under HC1 -- a new
# system on another phase shifts the neutral and can lower a phase-to-neutral voltage).
# :refuse_continue (default) = every applicant is tested in turn against everyone accepted
# so far; a breach means that applicant is refused and the next one is still tested (a DNSP
# assesses each application on its own). Gives BOTH the count before the first refusal and
# the total accepted, in exactly (n - k0) solves per order.
const ADOPT_METHOD = Symbol(get(ENV, "LV12_ADOPT_METHOD", "refuse_continue"))
const ADOPT_SEED   = 2026
const ADOPT_SCAN   = ADOPT_METHOD == :bisect   # crossing check only matters for bisection
# Each configuration applies to the WHOLE fleet (existing + new): HC0 = every inverter at
# unity PF; HC1 = every inverter with AS/NZS 4777.2 volt-var + volt-watt; HC2 = HC1 + device.
# Dispatch: "minvmax" = COORDINATED (upper bound: minimise the feeder's highest voltage,
# needs every bus); "local_vuf" = LOCAL (realistic: each device minimises negative sequence
# at its own bus; see add_local_objective! for local_vbal / local_vbal_vv); "none" = nothing to dispatch (no device: loads, source and every PV
# inverter are fixed by their control laws, so each solve is a plain power flow and the
# objective cannot change the answer; the cheapest objective is built).
const ADOPT_CONFIGS = [   # (label, pv_mode, statcom_mode, kVA per device, sites, dispatch)
    ("HC0",              :unity,   :none, 0.0,  String[],           "none"),
    ("HC1",              HC1_MODE, :none, 0.0,  String[],           "none"),
    ("HC2_D_one_coord",  HC1_MODE, :pq,   30.0, ["b3278"],          "minvmax"),
    ("HC2_D_one_local",  HC1_MODE, :pq,   30.0, ["b3278"],          "local_vuf"),
    ("HC2_D_pair_coord", HC1_MODE, :pq,   30.0, ["b3278", "b3162"], "minvmax"),
    ("HC2_D_pair_local", HC1_MODE, :pq,   30.0, ["b3278", "b3162"], "local_vuf"),
    # further local controllers (15 Oct): voltage balancing incl. zero sequence, and
    # voltage balancing + positive-sequence volt-var (current balancing dropped 15 Oct)
    ("HC2_D_one_vbal",     HC1_MODE, :pq, 30.0, ["b3278"],          "local_vbal"),
    ("HC2_D_one_vbalvv",   HC1_MODE, :pq, 30.0, ["b3278"],          "local_vbal_vv"),
    ("HC2_D_pair_vbal",    HC1_MODE, :pq, 30.0, ["b3278", "b3162"], "local_vbal"),
    ("HC2_D_pair_vbalvv",  HC1_MODE, :pq, 30.0, ["b3278", "b3162"], "local_vbal_vv"),
    # device at b800, the LV busbar just downstream of the transformer leads: coordinated
    # upper bound for a substation-sited unit (relevant where thermal limits bind, 70-80%)
    ("HC2_D_head_coord",   HC1_MODE, :pq, 30.0, ["b800"],           "minvmax"),
]

function adoption_orders(n; pen=ADOPT_PEN, n_init=ADOPT_N_INIT, n_new=ADOPT_N_NEW, seed=ADOPT_SEED)
    k0 = round(Int, pen * n)
    out = Tuple{Int,Int,Int,Vector{Int}}[]
    # NESTED design (Fred, 14 Oct: "only add PV systems to that baseline, not start with a
    # complete new sample"). Fleet i has one master order M_i of all customers; the fleet at
    # penetration p is M_i[1:round(p*n)], so a higher-penetration fleet CONTAINS the lower one.
    # Applicant order j is set by a fixed random key per customer, u_ij, so applicants keep
    # the same relative order whatever p is. Changing p therefore only moves the line between
    # "existing" and "applicant"; nothing is re-sampled.
    for i in 1:n_init
        master = randperm(MersenneTwister(seed + 1000 * i), n)
        existing = master[1:k0]
        rest = master[k0+1:end]
        for j in 1:n_new
            key = rand(MersenneTwister(seed + 1000 * i + j), n)
            applicants = rest[sortperm(key[rest])]
            push!(out, (i, j, k0, vcat(existing, applicants)))
        end
    end
    return out
end

# Does the network pass with exactly the customers `idx` carrying PV?
function adopt_pass_set(idx, names, cfg, tag)
    _, pvm, scm, kva, sites, disp = cfg
    return hc_point(ADOPT_KVA; pv_mode=pvm, statcom_mode=scm, stage=tag,
                    statcom_kva = scm == :none ? STATCOM_KVA : kva,
                    statcom_bus = scm == :none ? STATCOM_BUS : sites,
                    pv_customers = names[idx], objective = disp)
end

# Refuse and continue: returns (already_over, first = applicants accepted before the first
# refusal, accepted = indices of accepted applicants, refused = indices, first refusal binding).
function adopt_refuse_continue(perm, names, k0, cfg, tag)
    ns = 0
    have = perm[1:k0]
    r0 = (ns += 1; adopt_pass_set(have, names, cfg, "$tag base"))
    r0.pass || return (already_over=true, first=0, accepted=Int[], refused=perm[k0+1:end],
                       binding=r0.binding, n=ns, m=r0.m)
    accepted, refused, first, bind, mlast = Int[], Int[], -1, "all connect", r0.m
    for (q, a) in enumerate(perm[k0+1:end])
        r = (ns += 1; adopt_pass_set(vcat(have, accepted, [a]), names, cfg, "$tag app=$q"))
        if r.pass
            push!(accepted, a); mlast = r.m
        else
            push!(refused, a)
            if first < 0
                first, bind = length(accepted), r.binding
            end
        end
    end
    first < 0 && (first = length(accepted))
    return (already_over=false, first=first, accepted=accepted, refused=refused, binding=bind, n=ns, m=mlast)
end

# Does the network pass with the first k customers of `perm` carrying PV?
function adopt_pass(k, perm, names, cfg, tag)
    _, pvm, scm, kva, sites, disp = cfg
    return hc_point(ADOPT_KVA; pv_mode=pvm, statcom_mode=scm, stage=tag,
                    statcom_kva = scm == :none ? STATCOM_KVA : kva,
                    statcom_bus = scm == :none ? STATCOM_BUS : sites,
                    pv_customers = names[perm[1:k]], objective = disp)
end

# Number of customers connected: existing fleet k0, then applicants in order until the
# first one whose connection breaches a limit (ADOPT_METHOD = :sequential), or the
# largest passing k by bisection (:bisect; only valid if pass/fail is monotone in k).
function adopt_count(perm, names, k0, cfg, tag)
    n, ns = length(perm), 0
    t(k) = (ns += 1; adopt_pass(k, perm, names, cfg, "$tag k=$k"))
    r0 = t(k0)
    r0.pass || return (k=k0, already_over=true, binding=r0.binding, n=ns)
    if ADOPT_METHOD == :sequential
        for k in k0+1:n
            r = t(k)
            r.pass || return (k=k-1, already_over=false, binding=r.binding, n=ns)
        end
        return (k=n, already_over=false, binding="all connect", n=ns)
    end
    rN = t(n)
    rN.pass && return (k=n, already_over=false, binding="all connect", n=ns)
    lo, hi, bind = k0, n, rN.binding
    while hi - lo > 1
        mid = (lo + hi) ÷ 2
        r = t(mid)
        if r.pass
            lo = mid
        else
            hi, bind = mid, r.binding
        end
    end
    return (k=lo, already_over=false, binding=bind, n=ns)
end

adopt_rows = Vector{Any}[]
scan_rows  = Vector{Any}[]
if 11 in STAGES
    dmn   = base_network(false)
    lds   = sort(collect(values(dmn["load"])), by = load_name)
    cust_names = [load_name(l) for l in lds]
    ph_of  = Dict(load_name(l) => only(filter(c -> c != 4, l["connections"])) for l in lds)
    bus_of = Dict(load_name(l) => get(dmn["bus"][string(l["load_bus"])], "name", "") for l in lds)
    orders = adoption_orders(length(cust_names))
    k0 = orders[1][3]
    println("\n STAGE 11 -- sequential adoption ($(ADOPT_METHOD)): $k0 of $(length(cust_names)) customers start with ",
            "$(ADOPT_KVA) kVA ($(round(Int, 100 * ADOPT_PEN))%); $(length(orders)) orders = ",
            "$ADOPT_N_INIT fleets x $ADOPT_N_NEW applicant orders")
    write_csv(joinpath(OUTDIR, "lv12_adoption_orders.csv"), ["init", "order", "k0", "perm_load_names"],
              Vector{Any}[Any[i, j, kk, join(cust_names[p], " ")] for (i, j, kk, p) in orders])
    phase_counts(set) = [count(nm -> ph_of[nm] == p, set) for p in 1:3]

    for cfg in ADOPT_CONFIGS
        label = cfg[1]
        isempty(ADOPT_ONLY) || label in ADOPT_ONLY || continue
        println("\n   configuration $label")
        if ADOPT_SCAN
            (i, j, kk, perm) = orders[1]
            flags = [adopt_pass(k, perm, cust_names, cfg, "11scan $label k=$k").pass for k in kk:length(perm)]
            crossings = count(q -> flags[q] != flags[q + 1], 1:length(flags) - 1)
            first_fail = findfirst(!, flags)
            push!(scan_rows, Any[label, i, j, kk, crossings,
                                 isnothing(first_fail) ? length(perm) : kk + first_fail - 2])
            println("     single-crossing scan, order (1,1): $crossings crossing(s) over k = $kk..$(length(perm))")
        end
        hc_new = Int[]; over = 0
        if ADOPT_METHOD == :refuse_continue
            n_acc = Int[]; n_first = Int[]; over = 0
            for (i, j, kk, perm) in orders
                r = adopt_refuse_continue(perm, cust_names, kk, cfg, "11 $label i=$i j=$j")
                acc_nm, ref_nm = cust_names[r.accepted], cust_names[r.refused]
                fr = isempty(r.refused) ? "" : cust_names[r.refused[1]]
                push!(adopt_rows, Any[label, i, j, kk, r.already_over, r.first, length(r.accepted),
                                      length(r.refused), phase_counts(cust_names[perm[1:kk]])...,
                                      phase_counts(acc_nm)..., phase_counts(ref_nm)...,
                                      fr, fr == "" ? "" : bus_of[fr], fr == "" ? "" : ph_of[fr],
                                      r.binding, join(ref_nm, " "),
                                      isnothing(r.m) ? nothing : round(r.m.iph_max_a, digits=1),
                                      isnothing(r.m) ? "" : r.m.iph_at,
                                      isnothing(r.m) ? nothing : round(r.m.src_kva, digits=1),
                                      isnothing(r.m) ? nothing : round(r.m.vmax_V, digits=2), r.n])
                push!(n_acc, length(r.accepted)); push!(n_first, r.first); over += r.already_over
                @printf("     fleet %2d order %2d: %s  accepted %2d of %d  (before first refusal: %2d)   first refused: %-14s %s\n",
                        i, j, r.already_over ? "ALREADY OVER" : "            ", length(r.accepted),
                        length(perm) - kk, r.first, fr == "" ? "-" : "$(bus_of[fr]) ph$(ph_of[fr])", r.binding)
            end
            qa(p) = quantile(n_acc, p); qf(p) = quantile(n_first, p)
            @printf("   %s: accepted median %.1f (p10 %.1f, p90 %.1f); before first refusal median %.1f (p10 %.1f, p90 %.1f); already over: %d of %d\n",
                    label, qa(0.5), qa(0.1), qa(0.9), qf(0.5), qf(0.1), qf(0.9), over, length(orders))
            continue
        end
        for (i, j, kk, perm) in orders
            r = adopt_count(perm, cust_names, kk, cfg, "11 $label i=$i j=$j")
            new_k = r.already_over ? 0 : r.k - kk
            refused = (!r.already_over && r.k < length(perm)) ? cust_names[perm[r.k + 1]] : ""
            push!(adopt_rows, Any[label, i, j, kk, r.already_over, r.already_over ? nothing : r.k, new_k,
                                  phase_counts(cust_names[perm[1:kk]])...,
                                  phase_counts(cust_names[perm[kk+1:kk+new_k]])...,
                                  refused, refused == "" ? "" : bus_of[refused],
                                  refused == "" ? "" : ph_of[refused], r.binding, r.n])
            push!(hc_new, new_k); over += r.already_over
            @printf("     fleet %2d order %2d: %s  +%2d connect   refused: %-14s %s\n", i, j,
                    r.already_over ? "ALREADY OVER" : "            ", new_k,
                    refused == "" ? "-" : "$(bus_of[refused]) ph$(ph_of[refused])", r.binding)
        end
        q(p) = quantile(hc_new, p)
        @printf("   %s: new connections median %.1f (p10 %.1f, p90 %.1f) of %d applicants; already over at %d%%: %d of %d orders\n",
                label, q(0.5), q(0.1), q(0.9), length(cust_names) - k0, round(Int, 100 * ADOPT_PEN), over, length(orders))
    end
end

# --- OUTPUT --------------------------------------------------------------
println()
isempty(adopt_rows) || write_csv(joinpath(OUTDIR, "lv12_adoption.csv"),
    ADOPT_METHOD == :refuse_continue ?
    ["config", "init", "order", "k0", "already_over", "n_before_first_refusal", "n_accepted",
     "n_refused", "exist_A", "exist_B", "exist_C", "acc_A", "acc_B", "acc_C", "ref_A", "ref_B",
     "ref_C", "first_refused", "first_refused_bus", "first_refused_phase", "first_fail",
     "refused_list", "final_iph_max_a", "final_iph_at", "final_src_kva", "final_vmax_v", "solves"] :
    ["config", "init", "order", "k0", "already_over", "k_total", "k_new",
     "exist_A", "exist_B", "exist_C", "new_A", "new_B", "new_C",
     "first_refused", "refused_bus", "refused_phase", "first_fail", "solves"], adopt_rows)
isempty(scan_rows) || write_csv(joinpath(OUTDIR, "lv12_adoption_scan.csv"),
    ["config", "init", "order", "k0", "crossings", "k_scan"], scan_rows)
isempty(base0_rows) || write_csv(joinpath(OUTDIR, "lv12_base0.csv"),
    ["sites", "n_devices", "kva_each", "mode", "hc_kva_per_customer", "dhc_vs_hc0", "capped",
     "active", "first_fail", "util", "solves"], base0_rows)
isempty(pair_rows) || write_csv(joinpath(OUTDIR, "lv12_pairs.csv"),
    ["sites", "n_devices", "kva_each", "kva_total", "mode", "hc_kva_per_customer", "dhc_vs_hc1",
     "capped", "active", "first_fail", "util", "solves"], pair_rows)
isempty(siting_rows) || write_csv(joinpath(OUTDIR, "lv12_siting.csv"),
    ["site", "kva", "mode", "hc_kva_per_customer", "dhc_vs_hc1", "capped", "active",
     "first_fail", "util", "solves"], siting_rows)
write_csv(joinpath(OUTDIR, "lv12_stage_checks.csv"),
          ["stage", "check", "value", "target", "tolerance", "pass"], STAGE_ROWS)
write_csv(joinpath(OUTDIR, "lv12_hc.csv"),
          ["stage", "pv_mode", "statcom", "code", "hc_kva_per_customer", "capped", "active",
           "first_fail", "solves", "dhc_vs_hc1"], hc_rows)
write_csv(joinpath(OUTDIR, "lv12_points.csv"), POINT_HEADER, point_rows)

@printf("\n solves %d | parse %.1f s | build %.1f s | optimize %.1f s (Ipopt %.1f s)\n",
        Int(TIMING[:n]), TIMING[:parse], TIMING[:build], TIMING[:optimize], TIMING[:ipopt])
println("Run finished: ", Libc.strftime("%Y-%m-%d %H:%M:%S", time()),
        @sprintf("   (%.1f min)", (time() - RUN_T0) / 60))