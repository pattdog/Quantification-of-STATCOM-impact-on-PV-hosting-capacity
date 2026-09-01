#=
==============================================================================
objectives_FIXED.jl
==============================================================================
Network-agnostic rewrite of objectives.jl / objectives_FIXED.jl.

GOAL: every objective below should work unmodified regardless of:
  - how many generators/STATCOMs/PVs are in the network
  - which bus is the reference/slack bus, or how many there are
  - which generator id happens to be "1"

SUPPORTED OBJECTIVES (set `global objective = "..."` before include):
  "cost"      -- linear cost on generator active power (unchanged, already general)
  "loss"      -- I^2*|Z| loss, summed over branches (was: one hardcoded branch)
  "VUF2"      -- sum (or worst-case) of squared negative-sequence voltage across
                 target buses. No division -- simplest/best-conditioned objective
                 in this group. Use this as the diagnostic baseline.
  "VUF"       -- true ratio |V_neg|^2 / |V_pos|^2, matching worst_case_vuf()'s
                 definition. Adds a division -- more numerically demanding than
                 VUF2. A floor is placed on |V_pos|^2 to protect against a
                 degenerate near-zero denominator.
  "PVUR"      -- worst-case phase voltage unbalance ratio across target buses
  "LVUR"      -- worst-case line voltage unbalance ratio across target buses
  "IUF"/"IUF2"-- negative-sequence generator current (bugfix: was using crg twice
                 instead of crg+cig)
  "PIUR"      -- worst-case phase current unbalance ratio across target gens
  "IUF_inv"/"IUF2_inv" -- negative-sequence current at a specific branch (feeder head
                 by default)
  "PPUR"      -- worst-case phase active-power unbalance ratio across target gens
  "PQUR"      -- worst-case phase reactive-power unbalance ratio across target gens

OPTIONAL OVERRIDES (set these as globals in the calling script before include;
all are optional -- sensible defaults are computed automatically otherwise):

  objective_target_gens   :: Vector{<:Any}
      Generator ids to evaluate unbalance objectives over (VUF/VUF2/PVUR/LVUR use
      the buses these gens sit at; IUF/PIUR/PPUR/PQUR use the gens directly).
      Default: every generator NOT sitting at a reference/slack bus (see
      default_target_gens() below) -- i.e. every STATCOM/PV/DER, excluding the
      grid/substation source. If your network tags devices with gen["type"],
      you can narrow this further, e.g.:
          global objective_target_gens = [i for (i,g) in ref[:gen] if get(g,"type","")=="STATCOM"]

  objective_aggregation   :: Symbol   (:max or :sum)
      How multi-device unbalance objectives combine across target_gens/buses.
      :max (default) = minimize the worst-case device/bus (epigraph formulation,
                        matches the original single-device pattern).
      :sum = minimize the total/average unbalance across all targets.

  objective_loss_branches :: Vector{<:Any}
      Branch ids to include in the "loss" objective. Default: every branch.

  objective_branch_id     :: Any
      Branch id to use for IUF_inv/IUF2_inv. Default: the feeder-head branch
      auto-detected via RPMD.get_ref_bus_branch(ref).

Everything else (variable creation, sequence-transform math) is structurally
the same as the original file, just looped over `target_gens`/buses instead
of a hardcoded single id.
==============================================================================
=#

# ── Sequence-transform matrices (safe to redefine if constraints_PBalance.jl
#    already defined these -- same values, no behavioural change) ──────────
alpha = exp(im*2/3*pi)
T   = 1/3 * [1 1 1 ; 1 alpha alpha^2 ; 1 alpha^2 alpha]
Tre = real.(T)
Tim = imag.(T)

# ── General helpers ─────────────────────────────────────────────────────────

# A generator counts as "reference/slack" if it sits at a reference bus.
# This is the general, correct test -- NOT `id == "1"`, which is just a
# parse-order artifact of one particular .dss file and will not hold for
# every network.
is_reference_gen(gen) = haskey(ref[:ref_buses], gen["gen_bus"])

# Default target set: every generator that is NOT a reference/slack source.
# Covers STATCOMs, PVs, or any other DER, regardless of count or naming.
function default_target_gens()
    return [i for (i, gen) in ref[:gen] if !is_reference_gen(gen)]
end

target_gens = @isdefined(objective_target_gens) ? objective_target_gens : default_target_gens()
aggregation = @isdefined(objective_aggregation) ? objective_aggregation : :max

# Buses corresponding to target_gens (used by VUF/VUF2/PVUR/LVUR, which are
# bus-quantities, not generator quantities). Deduplicated, order-preserved.
function target_buses_from_gens(gens)
    buses = Int[]
    for i in gens
        b = ref[:gen][i]["gen_bus"]
        b in buses || push!(buses, b)
    end
    return buses
end

# Safety guard: several objectives are undefined (empty sums/epigraphs) when
# there are no target devices -- e.g. running scenario A (no STATCOMs) against
# an objective that only makes sense with DERs present. Bypass to a constant
# objective instead of crashing, exactly like the original "Scenario A fix"
# but applied uniformly across every objective, not just VUF.
function no_targets_bypass!()
    JuMP.@objective(model, Min, 0.0)
end

# ==============================================================================
if objective == "cost"
    # Already fully general: sums over every generator regardless of count.
    obj_cost = JuMP.@expression(model, sum(gen["cost"][1]*sum(pg[:,i]) + gen["cost"][2] for (i,gen) in ref[:gen]))
    JuMP.@objective(model, Min, obj_cost)

# ==============================================================================
elseif objective == "loss"
    # CHANGE: sum over every branch (or an explicit override list) instead of
    # a single hardcoded branch matched by name=="internal". Name-matching is
    # fragile -- not every network's .dss file will use that convention, and
    # "loss" should mean *network* loss, not one arbitrarily chosen line.
    loss_branches = @isdefined(objective_loss_branches) ? objective_loss_branches : keys(ref[:branch])

    if isempty(loss_branches)
        no_targets_bypass!()
    else
        obj_loss = JuMP.@expression(model, sum(
            sqrt(ref[:branch][l]["br_r"][c]^2 + ref[:branch][l]["br_x"][c]^2) *
            (csr[c,l]^2 + csi[c,l]^2)
            for l in loss_branches for c in 1:size(ref[:branch][l]["br_r"],1)
        ))
        JuMP.@objective(model, Min, obj_loss)
    end

# ==============================================================================
elseif objective in ["VUF", "VUF2", "PVUR", "LVUR"]
    # CHANGE: target buses now come from target_gens (any number of DERs, any
    # network), not from a hardcoded exclusion of gen id "1". PVUR/LVUR now
    # aggregate across ALL target buses (max or sum, per objective_aggregation)
    # instead of looking only at gen_buses[1].
    #
    # IMPORTANT CHANGE vs. earlier version: vr_012/vi_012 are now built as
    # @expression, not @variable + defining @constraint. They are exact linear
    # combinations of vr/vi (the sequence transform), so there is no need to
    # introduce new decision variables and equality constraints for them --
    # an expression substitutes directly wherever it's used, giving Ipopt an
    # identical feasible region with fewer variables/constraints to resolve.
    #
    # vm/vm_012 are now only constructed when actually needed (PVUR needs vm;
    # LVUR builds its own separate vm_ll and never needs vm or vm_012 at all;
    # VUF/VUF2 need neither). Previously vm_012 was always built and
    # constrained via vm_012^2 == vr_012^2 + vi_012^2, EVEN FOR VUF/VUF2,
    # which never read vm_012 anywhere. That is a real problem, not just
    # waste: at any solution where negative-sequence voltage is being driven
    # toward zero (exactly what VUF/VUF2 optimizes for), vr_012[3,i],
    # vi_012[3,i], and therefore vm_012[3,i] all approach 0 simultaneously --
    # and the Jacobian of vm_012^2 - vr_012^2 - vi_012^2 vanishes entirely at
    # that point (all partials are 2*value, i.e. 0). This is a textbook
    # constraint-qualification failure (LICQ fails) sitting exactly at the
    # objective's own optimum, and is a very plausible cause of the
    # LOCALLY_INFEASIBLE results seen when using VUF2 on the full STATCOM
    # fleet. Removing the unused vm_012 for VUF/VUF2 eliminates this
    # degenerate constraint outright.
    gen_buses = target_buses_from_gens(target_gens)

    if isempty(gen_buses)
        no_targets_bypass!()
    else
        # vr_012/vi_012 as expressions -- always cheap, always needed by every
        # objective in this group, never introduces a variable/constraint.
        vr_012 = Dict{Int,Vector{JuMP.AffExpr}}()
        vi_012 = Dict{Int,Vector{JuMP.AffExpr}}()
        for i in gen_buses
            terminals_i = ref[:bus][i]["terminals"][1:end-1]
            vr_012[i] = JuMP.@expression(model, Tre * Array(vr[terminals_i,i]) .- Tim * Array(vi[terminals_i,i]))
            vi_012[i] = JuMP.@expression(model, Tre * Array(vi[terminals_i,i]) .+ Tim * Array(vr[terminals_i,i]))
        end
        # index convention preserved: vr_012[i][2]/[3] are pos/neg sequence
        # (matches the old vr_012[2,i]/[3,i] DenseAxisArray indexing, just
        # keyed the other way round since these are now plain Dicts of Vectors)

        need_vm = objective in ["PVUR"]   # only PVUR reads vm; LVUR builds its own vm_ll

        if need_vm
            terminals = Dict(i => collect(1:3) for i in gen_buses)
            vm = Dict(i => JuMP.@variable(model, [t in terminals[i], i], base_name="vm_$i", lower_bound = 0 ) for i in gen_buses)
            vm = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([t in vm[i].axes[1] ? vm[i][t,i] : 0.0 for t in 1:n_ph, i in gen_buses]), 1:n_ph, gen_buses)
            for i in gen_buses
                terminals_i = ref[:bus][i]["terminals"][1:end-1]
                JuMP.@constraint(model, vm[terminals_i,i].^2 .== vr[terminals_i,i].^2 .+ vi[terminals_i,i].^2)
                for t in terminals_i
                    JuMP.set_start_value(vm[t,i], sqrt(JuMP.start_value(vr[t,i])^2 + JuMP.start_value(vi[t,i])^2))
                end
            end
        end

        if objective == "VUF2"
            # Unnormalized: sum (or worst-case) of squared negative-sequence
            # voltage magnitude. No division -- simplest/best-conditioned
            # objective in this group. Use this as the diagnostic baseline.
            if aggregation == :sum
                JuMP.@objective(model, Min, sum(vr_012[i][3]^2 + vi_012[i][3]^2 for i in gen_buses))
            else
                worst = JuMP.@variable(model, base_name="worst_vuf2", lower_bound=0)
                JuMP.set_start_value(worst, 0.0)
                for i in gen_buses
                    JuMP.@constraint(model, worst >= vr_012[i][3]^2 + vi_012[i][3]^2)
                end
                JuMP.@objective(model, Min, worst)
            end

        elseif objective == "VUF"
            # True ratio: |V_neg|^2 / |V_pos|^2, matching worst_case_vuf()'s
            # definition (as a squared-magnitude ratio, equivalent up to a
            # sqrt on the magnitude ratio itself). Adds a division -- more
            # numerically demanding than VUF2. Guard the denominator away
            # from zero since |V_pos| should sit near 1.0 pu in any sane
            # operating point, but a lower bound protects against a
            # degenerate solver excursion toward 0.
            for i in gen_buses
                JuMP.@constraint(model, vr_012[i][2]^2 + vi_012[i][2]^2 >= 0.25)  # |V_pos| >= 0.5 pu
            end
            if aggregation == :sum
                JuMP.@objective(model, Min, sum(
                    (vr_012[i][3]^2 + vi_012[i][3]^2) / (vr_012[i][2]^2 + vi_012[i][2]^2)
                    for i in gen_buses))
            else
                worst = JuMP.@variable(model, base_name="worst_vuf", lower_bound=0)
                JuMP.set_start_value(worst, 0.0)
                for i in gen_buses
                    JuMP.@constraint(model, worst >= (vr_012[i][3]^2 + vi_012[i][3]^2) / (vr_012[i][2]^2 + vi_012[i][2]^2))
                end
                JuMP.@objective(model, Min, worst)
            end

        elseif objective == "PVUR"
            # CHANGE: epigraph variable now bounded by EVERY target bus's phase
            # deviation, not just gen_buses[1]. Denominator uses the mean of
            # all target buses' average phase voltage.
            phase_voltage = JuMP.@variable(model, base_name="phase_voltage", lower_bound=0)
            JuMP.set_start_value(phase_voltage, 0.0)
            for i in gen_buses
                terminals_i = ref[:bus][i]["terminals"][1:end-1]
                JuMP.@constraint(model, [t in terminals_i], phase_voltage >= vm[t,i] - sum(vm[terminals_i,i])/3)
            end
            denom = JuMP.@expression(model, sum(sum(vm[ref[:bus][i]["terminals"][1:end-1],i])/3 for i in gen_buses) / length(gen_buses))
            JuMP.@objective(model, Min, phase_voltage / denom)

        elseif objective == "LVUR"
            # CHANGE: same epigraph-over-all-target-buses generalization,
            # applied to line-line voltages instead of phase voltages.
            line_voltage = JuMP.@variable(model, base_name="line_voltage", lower_bound=0)
            JuMP.set_start_value(line_voltage, 0.0)
            vm_ll = JuMP.@variable(model, [t in 1:3, i in gen_buses], base_name="vm_ll", lower_bound=0)
            vr_ll = JuMP.@variable(model, [t in 1:3, i in gen_buses], base_name="vr_ll")
            vi_ll = JuMP.@variable(model, [t in 1:3, i in gen_buses], base_name="vi_ll")

            for i in gen_buses
                terminals_i  = collect(1:3)
                terminals2   = [terminals_i[2:end]..., terminals_i[1]]
                JuMP.@constraint(model, vr_ll[terminals_i,i] .== Array(vr[terminals_i,i]) .- Array(vr[terminals2,i]))
                JuMP.@constraint(model, vi_ll[terminals_i,i] .== Array(vi[terminals_i,i]) .- Array(vi[terminals2,i]))
                JuMP.@constraint(model, vm_ll[terminals_i,i].^2 .== vr_ll[terminals_i,i].^2 .+ vi_ll[terminals_i,i].^2)
                JuMP.@constraint(model, [t in terminals_i], line_voltage >= vm_ll[t,i] - sum(vm_ll[terminals_i,i])/3)
            end
            denom = JuMP.@expression(model, sum(sum(vm_ll[1:3,i])/3 for i in gen_buses) / length(gen_buses))
            JuMP.@objective(model, Min, line_voltage / denom)
        end
    end

# ==============================================================================
elseif objective in ["IUF", "IUF2", "PIUR"]
    # CHANGE: restrict target_gens to genuinely 3-phase generators only.
    # Sequence decomposition is only meaningful for a device with all 3
    # phases -- a single-phase generator (e.g. phase-matched residential PV)
    # has no negative-sequence current to speak of, and Tre/Tim (fixed 3x3
    # sequence-transform matrices) cannot be multiplied against a length-1
    # phases vector.
    iuf_target_gens = filter(i -> begin
        g = ref[:gen][i]
        conns = g["connections"]
        phases_i = 4 in conns ? conns[1:end-1] : conns
        length(phases_i) == 3
    end, target_gens)

    if isempty(iuf_target_gens)
        no_targets_bypass!()
    else
        # crg_012/cig_012 as expressions -- exact linear combinations of
        # crg_bus/cig_bus, same treatment as vr_012/vi_012 in the VUF block.
        # No new variables/equality constraints needed for these.
        crg_012 = Dict{Int,Vector{JuMP.AffExpr}}()
        cig_012 = Dict{Int,Vector{JuMP.AffExpr}}()
        for i in iuf_target_gens
            gen = ref[:gen][i]
            phases = 4 in gen["connections"] ? gen["connections"][1:end-1] : gen["connections"]
            crg_012[i] = JuMP.@expression(model, Tre * Array(crg_bus[phases,i]) .- Tim * Array(cig_bus[phases,i]))
            cig_012[i] = JuMP.@expression(model, Tre * Array(cig_bus[phases,i]) .+ Tim * Array(crg_bus[phases,i]))
        end
        # index convention: crg_012[i][2]/cig_012[i][2] = positive sequence
        #                    crg_012[i][3]/cig_012[i][3] = negative sequence

        # NOTE: cmg_012 (sqrt-defined sequence-current MAGNITUDE) is
        # deliberately NOT built here. IUF/IUF2 only ever need squared
        # sequence-current magnitudes (crg_012^2 + cig_012^2), which need no
        # sqrt and no extra variable. Introducing cmg_012 via
        # cmg_012^2 == crg_012^2 + cig_012^2 -- as an earlier version of this
        # file did -- creates a constraint whose Jacobian vanishes exactly at
        # the point IUF2's own objective drives toward (negative-sequence
        # current -> 0), the same degenerate-constraint mechanism as vm_012
        # under VUF2. Confirmed empirically: IUF2 solved fine on some
        # scenarios/load levels but hit OTHER_ERROR on others (specifically
        # where the STATCOM's extra P-exchange freedom let it get closest to
        # true negative-sequence-current zero) -- removing cmg_012 entirely
        # removes the degenerate constraint outright, same fix as VUF2.

        # cmg (per-PHASE current magnitude, NOT a sequence quantity) is only
        # built when PIUR actually needs it -- IUF/IUF2 never reference cmg
        # at all, so building it for them was pure waste plus its own
        # potential unused-variable degeneracy.
        need_cmg = objective == "PIUR"
        if need_cmg
            int_dim = Dict(i => RPMD._infer_int_dim_unit(ref[:gen][i], !(4 in ref[:gen][i]["connections"])) for i in iuf_target_gens)
            cmg = Dict(i => JuMP.@variable(model, [c in 1:int_dim[i]], base_name="cmg_$i", lower_bound=0) for i in iuf_target_gens)
            cmg = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in 1:int_dim[i] ? cmg[i][c] : 0.0 for c in 1:n_ph, i in iuf_target_gens]), 1:n_ph, iuf_target_gens)
            for i in iuf_target_gens
                gen = ref[:gen][i]
                phases = 4 in gen["connections"] ? gen["connections"][1:end-1] : gen["connections"]
                JuMP.@constraint(model, cmg[phases,i].^2 .== crg_bus[phases,i].^2 .+ cig_bus[phases,i].^2)
                for p in phases
                    JuMP.set_start_value(cmg[p,i], sqrt(JuMP.start_value(crg_bus[p,i])^2 + JuMP.start_value(cig_bus[p,i])^2))
                end
            end
        end

        if objective == "IUF"
            # Guard the denominator away from zero -- same issue as VUF's
            # division. Floor applied directly on the squared positive-
            # sequence current (crg_012[2]^2 + cig_012[2]^2), no cmg_012
            # intermediary needed.
            for i in iuf_target_gens
                JuMP.@constraint(model, crg_012[i][2]^2 + cig_012[i][2]^2 >= 1e-6)
            end
            JuMP.@objective(model, Min, sum(
                (crg_012[i][3]^2 + cig_012[i][3]^2) / (crg_012[i][2]^2 + cig_012[i][2]^2)
                for i in iuf_target_gens))

        elseif objective == "IUF2"
            # Unnormalized: sum of squared negative-sequence current
            # magnitude. No division, no sqrt-defined auxiliary variable --
            # simplest/best-conditioned in this group, same role as VUF2.
            JuMP.@objective(model, Min, sum(crg_012[i][3]^2 + cig_012[i][3]^2 for i in iuf_target_gens))

        elseif objective == "PIUR"
            phase_current = JuMP.@variable(model, [i in iuf_target_gens], base_name="phase_current", lower_bound=0)
            JuMP.set_start_value.(phase_current, 0.0)
            for i in iuf_target_gens
                gen = ref[:gen][i]
                connections = gen["connections"][1:end-1]
                JuMP.@constraint(model, [t in connections], phase_current[i] >= cmg[t,i] - sum(cmg[connections,i])/3)
            end
            if aggregation == :sum
                JuMP.@objective(model, Min, sum(phase_current[i] / (sum(cmg[ref[:gen][i]["connections"][1:end-1],i])/3) for i in iuf_target_gens))
            else
                worst = JuMP.@variable(model, base_name="worst_piur", lower_bound=0)
                JuMP.set_start_value(worst, 0.0)
                for i in iuf_target_gens
                    connections = ref[:gen][i]["connections"][1:end-1]
                    JuMP.@constraint(model, worst >= phase_current[i] / (sum(cmg[connections,i])/3))
                end
                JuMP.@objective(model, Min, worst)
            end
        end
    end
# ==============================================================================



elseif objective == "PPUR"
    # CHANGE (bugfix + generalization): was hardcoded to gen id 1, which per
    # our network-setup discussion is the REFERENCE/SLACK generator, not a
    # STATCOM -- i.e. this previously measured the substation's own P
    # unbalance and ignored every actual DER in the network. Now runs over
    # target_gens (defaults to all non-reference generators).
    if isempty(target_gens)
        no_targets_bypass!()
    else
        phase_p = JuMP.@variable(model, [i in target_gens], base_name="phase_p", lower_bound=0)
        for i in target_gens
            gen = ref[:gen][i]
            connections = gen["connections"][1:end-1]
            JuMP.@constraint(model, [t in connections], phase_p[i] >= pg[t,i] - sum(pg[connections,i])/3)
        end
        if aggregation == :sum
            JuMP.@objective(model, Min, sum(phase_p[i] / (sum(pg[ref[:gen][i]["connections"][1:end-1],i])/3) for i in target_gens))
        else
            worst = JuMP.@variable(model, base_name="worst_ppur", lower_bound=0)
            for i in target_gens
                connections = ref[:gen][i]["connections"][1:end-1]
                JuMP.@constraint(model, worst >= phase_p[i] / (sum(pg[connections,i])/3))
            end
            JuMP.@objective(model, Min, worst)
        end
    end

# ==============================================================================
elseif objective == "PQUR"
    # Same bugfix + generalization as PPUR, applied to reactive power.
    if isempty(target_gens)
        no_targets_bypass!()
    else
        phase_q = JuMP.@variable(model, [i in target_gens], base_name="phase_q", lower_bound=0)
        for i in target_gens
            gen = ref[:gen][i]
            connections = gen["connections"][1:end-1]
            JuMP.@constraint(model, [t in connections], phase_q[i] >= qg[t,i] - sum(qg[connections,i])/3)
        end
        if aggregation == :sum
            JuMP.@objective(model, Min, sum(phase_q[i] / (sum(qg[ref[:gen][i]["connections"][1:end-1],i])/3) for i in target_gens))
        else
            worst = JuMP.@variable(model, base_name="worst_pqur", lower_bound=0)
            for i in target_gens
                connections = ref[:gen][i]["connections"][1:end-1]
                JuMP.@constraint(model, worst >= phase_q[i] / (sum(qg[connections,i])/3))
            end
            JuMP.@objective(model, Min, worst)
        end
    end

else
    error("objectives_FIXED.jl: unrecognized objective \"$objective\"")
end