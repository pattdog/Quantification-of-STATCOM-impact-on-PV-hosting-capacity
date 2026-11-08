#=
==============================================================================
TIER 2 TEST: STATCOM inter-phase P exchange on a minimal unbalanced network
==============================================================================
Purpose:
    Tier 1 proved the constraint FORM is sound in isolation (free pg/qg vars).
    This test embeds the same two constraints into a minimal real network --
    slack source, series line impedance, unbalanced constant-impedance loads,
    and a STATCOM at the PCC using the SAME pg/qg-from-voltage-and-current
    relations as your real constraints.jl:

        pg[p] ==  vr[p]*crg[p] + vi[p]*cig[p]
        qg[p] == -vr[p]*cig[p] + vi[p]*crg[p]

    We solve the SAME network twice:
        Scenario A (baseline): pg[p] == 0 for all phases (today's model --
                                STATCOM is Q-only, independent per phase)
        Scenario B (extension): S^2 capability + sum(pg)==0 coupling (the
                                 new physics), pg free to move between phases

    and compare the negative-sequence voltage at the PCC -- using the EXACT
    same symmetric-component transform (T, Tre, Tim) from your real
    constraints.jl -- to confirm Scenario B achieves equal-or-better voltage
    balancing than Scenario A. If it doesn't improve, something is wrong
    before this ever touches the real codebase.

Network (single PCC bus, 3-wire, no neutral, phases decoupled for simplicity):

    slack (fixed V) --[R+jX]-- PCC (unknown V) --- unbalanced load
                                    |
                                STATCOM (crg, cig)

    KCL per phase at PCC:  Y*(V1 - V2) + Igen = Iload
    (branch current in, plus generator current in, equals load current out)

Run with:
    julia test_statcom_tier2_network.jl
==============================================================================
=#

using JuMP
using Ipopt
using LinearAlgebra

# -----------------------------
# Symmetric component transform -- IDENTICAL to your real constraints.jl
# -----------------------------
alpha = exp(im*2/3*pi)
T = 1/3 * [1 1 1 ; 1 alpha alpha^2 ; 1 alpha^2 alpha]
Tre = real.(T)
Tim = imag.(T)

# -----------------------------
# Network parameters
# -----------------------------
# Slack (source) voltage -- balanced positive sequence, 1.0 pu
Vslack_mag = 1.0
theta = [0.0, -2pi/3, 2pi/3]
vr1 = Vslack_mag .* cos.(theta)
vi1 = Vslack_mag .* sin.(theta)

# Series line impedance (same per phase -- decoupled phases, crude simplification)
R = 0.02
X = 0.08
Zline = R + im*X
Yline = 1 / Zline
Yr, Yi = real(Yline), imag(Yline)

# Unbalanced constant-impedance load, specified as apparent power at ~1.0 pu
# voltage (heavy load on phase 1, light on 2 and 3 -- the imbalance we want
# the STATCOM to help correct)
Sload = [0.9 + 0.3im, 0.30 + 0.10im, 0.30 + 0.10im]   # [pu], phase 1,2,3
Yload = conj.(Sload)                                   # since S=V*conj(I), |V|~1 => Y=conj(S)
Yload_r = real.(Yload)
Yload_i = imag.(Yload)

# STATCOM rating
s_rated = 0.15   # pu per phase -- deliberately modest vs. the ~0.6pu imbalance

# -----------------------------
# Model builder
# -----------------------------
function build_and_solve(; allow_p_exchange::Bool)
    model = Model(Ipopt.Optimizer)
    set_silent(model)

    @variable(model, vr2[1:3])
    @variable(model, vi2[1:3])
    @variable(model, -s_rated <= crg[1:3] <= s_rated)   # loose safety box
    @variable(model, -s_rated <= cig[1:3] <= s_rated)   # loose safety box

    # pg/qg as expressions of V and I -- SAME relations as your real constraints.jl
    @expression(model, pg[p=1:3],  vr2[p]*crg[p] + vi2[p]*cig[p])
    @expression(model, qg[p=1:3], -vr2[p]*cig[p] + vi2[p]*crg[p])

    # --- KCL at PCC, per phase ---
    # branch current in (from slack through line impedance)
    @expression(model, cr_branch[p=1:3], Yr*(vr1[p]-vr2[p]) - Yi*(vi1[p]-vi2[p]))
    @expression(model, ci_branch[p=1:3], Yr*(vi1[p]-vi2[p]) + Yi*(vr1[p]-vr2[p]))
    # load current out (constant impedance)
    @expression(model, cr_load[p=1:3], Yload_r[p]*vr2[p] - Yload_i[p]*vi2[p])
    @expression(model, ci_load[p=1:3], Yload_r[p]*vi2[p] + Yload_i[p]*vr2[p])

    @constraint(model, kcl_re[p=1:3], cr_branch[p] + crg[p] == cr_load[p])
    @constraint(model, kcl_im[p=1:3], ci_branch[p] + cig[p] == ci_load[p])

    if allow_p_exchange
        # --- NEW physics: joint S^2 capability + DC-bus coupling ---
        @constraint(model, cap[p=1:3], pg[p]^2 + qg[p]^2 <= s_rated^2)
        @constraint(model, dc_coupling, sum(pg[p] for p in 1:3) == 0)
    else
        # --- baseline: STATCOM is Q-only, independent per phase (today's model) ---
        @constraint(model, p_zero[p=1:3], pg[p] == 0)
        @constraint(model, qcap[p=1:3], qg[p]^2 <= s_rated^2)
    end

    # --- Objective: drive PCC negative-sequence voltage to zero ---
    @expression(model, v_neg_re, Tre[3,:]'*vr2 - Tim[3,:]'*vi2)
    @expression(model, v_neg_im, Tre[3,:]'*vi2 + Tim[3,:]'*vr2)
    @objective(model, Min, v_neg_re^2 + v_neg_im^2)

    optimize!(model)

    return (
        status = termination_status(model),
        vr2 = value.(vr2), vi2 = value.(vi2),
        pg = value.(pg), qg = value.(qg),
        vneg2 = value(v_neg_re)^2 + value(v_neg_im)^2,
    )
end

# -----------------------------
# Run both scenarios
# -----------------------------
resA = build_and_solve(allow_p_exchange=false)   # baseline: Q-only
resB = build_and_solve(allow_p_exchange=true)    # extension: P-exchange enabled

function vmag(vr, vi)
    return sqrt.(vr.^2 .+ vi.^2)
end

println("="^70)
println("SCENARIO A -- baseline (STATCOM Q-only, pg forced to 0)")
println("Solver status: ", resA.status)
vmA = vmag(resA.vr2, resA.vi2)
for p in 1:3
    println("  Phase $p:  |V| = ", round(vmA[p], digits=4),
            "   pg = ", round(resA.pg[p], digits=4),
            "   qg = ", round(resA.qg[p], digits=4))
end
println("  |V_neg|^2 (objective) = ", round(resA.vneg2, digits=6))
println()

println("="^70)
println("SCENARIO B -- extension (P-exchange enabled: S^2 cap + DC coupling)")
println("Solver status: ", resB.status)
vmB = vmag(resB.vr2, resB.vi2)
for p in 1:3
    println("  Phase $p:  |V| = ", round(vmB[p], digits=4),
            "   pg = ", round(resB.pg[p], digits=4),
            "   qg = ", round(resB.qg[p], digits=4))
end
println("  sum(pg) = ", round(sum(resB.pg), digits=6), "   (should be ~0)")
println("  |V_neg|^2 (objective) = ", round(resB.vneg2, digits=6))
println()

println("="^70)
println("COMPARISON")
improvement = resA.vneg2 - resB.vneg2
pct = resA.vneg2 > 0 ? 100*improvement/resA.vneg2 : 0.0
println("  Negative-sequence voltage^2:  A = ", round(resA.vneg2, digits=6),
        "   B = ", round(resB.vneg2, digits=6))
println("  Improvement from P-exchange:  ", round(improvement, digits=6),
        "  (", round(pct, digits=2), "% reduction)")
if resB.vneg2 <= resA.vneg2 + 1e-8
    println("  RESULT: PASS -- P-exchange achieved equal-or-better voltage balancing.")
else
    println("  RESULT: FAIL -- P-exchange did WORSE than Q-only baseline.")
    println("          This should be impossible (B is a relaxation of A's feasible")
    println("          region), so check: is dc_coupling actually active? Is the")
    println("          objective using the same v_neg expression in both builds?")
end
println("="^70)