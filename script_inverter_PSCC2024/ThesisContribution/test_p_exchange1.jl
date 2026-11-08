#=
==============================================================================
CRUDE STANDALONE TEST: STATCOM inter-phase active power exchange
==============================================================================
Purpose:
    Isolate and verify the TWO new constraints needed for the STATCOM P-exchange
    extension, with NO network physics, NO PMD, NO ref/variables/constraints.jl
    plumbing. Just pg/qg decision variables and the constraints themselves.

    (1) Per-phase joint capability constraint:
            pg[p]^2 + qg[p]^2 <= s_rated^2      for p in 1:3

    (2) DC-bus coupling constraint (Eq. 18 in your derivation):
            pg[1] + pg[2] + pg[3] == 0

    If this script behaves as expected, we know the constraint FORM is correct
    before splicing it into the real constraints.jl generator loop.

Test design:
    We give the optimizer an objective that WANTS an unbalanced P profile
    (e.g. draw power on phase 1, inject it on phase 3) and a Q profile that
    competes for the same S budget. We check:
        - the solver finds a feasible, unbalanced P allocation
        - sum(pg) == 0 to numerical tolerance
        - no phase exceeds S_rated
        - Q gets traded off against P when S is tight (confirms S^2 is really
          coupling P and Q, not acting as two independent boxes)

Run with:
    julia test_statcom_p_exchange.jl
==============================================================================
=#

using JuMP
using Ipopt

# -----------------------------
# Parameters
# -----------------------------
s_rated = 0.7                    # per-phase converter rating [pu]
p_target = [-0.6, 0.0, 0.6]      # what we're ASKING the STATCOM to do:
                                  # draw 0.6 pu off phase 1, push 0.6 pu into phase 3
                                  # (sums to zero by construction -- realistic ask)
q_target = [0.5, 0.5, 0.5]       # simultaneously wanting balanced Q support on all 3
                                  # phases -- large enough to compete with P for the
                                  # S budget on phase 1 and phase 3 (0.6^2+0.5^2=0.61>0.5)
                                  # i.e. phase 1 and 3 CANNOT hit both targets exactly --
                                  # this is deliberate, to see how the tradeoff resolves.

w_p = 1.0     # objective weight on tracking p_target
w_q = 1.0     # objective weight on tracking q_target

# -----------------------------
# Model
# -----------------------------
model = Model(Ipopt.Optimizer)
set_silent(model)   # comment out if you want Ipopt's solver log

@variable(model, -s_rated <= pg[1:3] <= s_rated)   # loose outer box (safety only)
@variable(model, -s_rated <= qg[1:3] <= s_rated)   # loose outer box (safety only)

# --- Constraint (1): per-phase joint S^2 capability ---
@constraint(model, cap[p=1:3], pg[p]^2 + qg[p]^2 <= s_rated^2)

# --- Constraint (2): DC-bus coupling across the three phases ---
@constraint(model, dc_coupling, sum(pg[p] for p in 1:3) == 0)

# --- Objective: track the (infeasible-as-stated) P and Q targets as closely
#     as possible, forcing the solver to reveal how it resolves the S^2 tradeoff ---
@objective(model, Min,
    w_p * sum((pg[p] - p_target[p])^2 for p in 1:3) +
    w_q * sum((qg[p] - q_target[p])^2 for p in 1:3)
)

optimize!(model)

# -----------------------------
# Report
# -----------------------------
println("Solver status: ", termination_status(model))
println()

pg_val = value.(pg)
qg_val = value.(qg)
S_val  = sqrt.(pg_val.^2 .+ qg_val.^2)

println("Phase |   P target |    P sol   |   Q target |   Q sol    |  S used  | S limit")
for p in 1:3
    println(rpad("  $p", 6),
        "| ", rpad(round(p_target[p], digits=4), 11),
        "| ", rpad(round(pg_val[p], digits=4), 11),
        "| ", rpad(round(q_target[p], digits=4), 11),
        "| ", rpad(round(qg_val[p], digits=4), 11),
        "| ", rpad(round(S_val[p], digits=4), 9),
        "| ", s_rated)
end

println()
p_sum = sum(pg_val)
println("sum(pg) across phases = ", p_sum, "   (should be ~0)")
println("max |sum(pg)| tolerance check: ", isapprox(p_sum, 0.0; atol=1e-6) ? "PASS" : "FAIL")

for p in 1:3
    ok = S_val[p] <= s_rated + 1e-6
    println("Phase $p capability check (S <= $s_rated): ", ok ? "PASS" : "FAIL", "  (S=$(round(S_val[p],digits=4)))")
end

println()
println("Interpretation:")
println("- If sum(pg) != 0, the coupling constraint isn't binding -- check constraint")
println("  is actually attached to `model` and uses the SAME pg container as the objective.")
println("- If any phase shows S > s_rated, the capability constraint isn't binding --")
println("  check it's a <=, not <, and that pg/qg in the constraint are the same")
println("  variables referenced in the objective (not stale copies).")
println("- Phase 1 and 3 should show Q pulled below 0.5 (their target) because hitting")
println("  both |P|=0.6 and Q=0.5 simultaneously violates S<=1.0 (0.6^2+0.5^2=0.61>1).")
println("  This confirms P and Q are genuinely coupled through a shared S budget,")
println("  not two independent boxes.")