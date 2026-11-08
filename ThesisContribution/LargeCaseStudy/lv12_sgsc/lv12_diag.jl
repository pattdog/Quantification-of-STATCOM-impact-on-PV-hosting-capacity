#=
lv12_diag.jl -- why does stage 3 (a plain power flow) come back LOCALLY_INFEASIBLE?

Run in the SAME Julia session, right after lv12_opf.jl stopped (its functions and globals are
still loaded):

    include("ThesisContribution/LargeCaseStudy/lv12_sgsc/lv12_diag.jl")

D0  Source setup in the math model (no solve): reference buses, grounding, vm/va, terminals.
D1  Line impedances as PMD parsed them, in ohms, against OpenDSS's own values. The LV12
    linecodes are per kft and the lines are in km (the 4-bus used km for both), so a unit
    slip would show here as a ratio of 3.28 or 1000.
D2  PMD's OWN IVRENPowerModel (no core files) on the identical model.
      solves      -> network, MEN and loads are fine; the problem is in the core-file build
      infeasible  -> the problem is in the network model itself (see D0/D1)
D3  The core files with every load set to zero.
      solves while D4 fails -> the load constraints, or the loading level, are the problem
D4  The core files as in stage 3, Ipopt log on, then the worst-violated constraints at
    Ipopt's final point: WHICH constraint cannot be met.
D5  The core files again, warm-started from D2's voltages (only if D2 solved).
      solves      -> it was the flat start; the fix is a warm start
      still fails -> a constraint in the core files is genuinely inconsistent here (see D4)

Send back everything printed from the "D0" banner on. The Ipopt log in D4 can be trimmed to
its last ~30 lines.
=#

banner(s) = println("\n", "="^96, "\n ", s, "\n", "="^96)

# terminals that a PMD solution vector covers (all, or only the ungrounded ones)
function _sol_terms(bus, n)
    t = bus["terminals"]
    n == length(t) && return t
    ug = t[(!).(bus["grounded"])]
    n == length(ug) && return ug
    return nothing
end

# --- D0 ------------------------------------------------------------------------------
banner("D0 -- source setup in the math model")
let dm0 = base_network(false)
    for (k, b) in sort(collect(dm0["bus"]), by = x -> parse(Int, x[1]))
        (b["bus_type"] == 3 || lowercase(get(b, "name", "")) == "b908") || continue
        println("  bus $k  name=$(get(b, "name", ""))  type=$(b["bus_type"])  terminals=$(b["terminals"])  grounded=$(b["grounded"])")
        vm = get(b, "vm", Float64[]); va = get(b, "va", Float64[])
        println("        vm=$(round.(vm, digits=4))  va(deg)=$(round.(rad2deg.(va), digits=1))  vbase=$(get(b, "vbase", "?"))  vmax=$(b["vmax"])")
    end
    for (k, g) in dm0["gen"]
        println("  gen $k  bus=$(g["gen_bus"])  connections=$(g["connections"])  configuration=$(g["configuration"])")
    end
    ntype3 = count(b -> b["bus_type"] == 3, values(dm0["bus"]))
    ngnd = [get(b, "name", k) for (k, b) in dm0["bus"] if any(b["grounded"])]
    println("  buses: $(length(dm0["bus"]))  type-3 buses: $ntype3  buses with a grounded terminal: $ngnd")
    println("  shunts: $(length(get(dm0, "shunt", Dict())))   branches: $(length(dm0["branch"]))   loads: $(length(dm0["load"]))")
    odd = [(get(b, "name", k), b["terminals"]) for (k, b) in dm0["bus"] if b["terminals"] != [1, 2, 3, 4]]
    println("  buses whose terminals are not [1,2,3,4]: ", isempty(odd) ? "none" : odd)
    vbr = [(get(br, "name", k), size(br["br_r"])) for (k, br) in dm0["branch"] if size(br["br_r"]) != (4, 4)]
    println("  branches whose impedance is not 4x4: ", isempty(vbr) ? "none" : vbr)
end

# --- D1 ------------------------------------------------------------------------------
banner("D1 -- line impedance as parsed by PMD vs OpenDSS (ohms, conductor 1 self terms)")
let dm0 = base_network(false)
    sb = get_sbase_kva(dm0)
    # OpenDSS values from Master_pmd_opendss_check.dss (R, X in ohm; C in nF)
    expect = Dict("l_4360" => (0.175518, 0.097652, 3.2993),
                  "l_10" => (0.146095, 0.081282, 2.7463),
                  "switch_114_closed" => (0.000150, 0.000731, 0.0330))
    for (k, br) in dm0["branch"]
        nm = lowercase(get(br, "name", ""))
        haskey(expect, nm) || continue
        vb_kV = dm0["bus"][string(br["f_bus"])]["vbase"]
        zb = (vb_kV * 1e3)^2 / (sb * 1e3)                          # ohm
        R, X = br["br_r"][1, 1] * zb, br["br_x"][1, 1] * zb
        Cn = 2 * br["b_fr"][1, 1] / zb / (2pi * 50) * 1e9            # nF (b_fr is half the total)
        eR, eX, eC = expect[nm]
        @printf("  %-19s R11 %.6f (OpenDSS %.6f, ratio %.3f)   X11 %.6f (ratio %.3f)   C11 %.4f nF (ratio %.3f)\n",
                nm, R, eR, R / eR, X, X / eX, Cn, Cn / eC)
    end
    println("  (zbase used: $(round((0.4e3 / sqrt(3))^2 / (sb * 1e3), digits=5)) ohm at 230.9 V, sbase $(sb) kVA)")
end

# --- D2 ------------------------------------------------------------------------------
banner("D2 -- PMD's own IVRENPowerModel on the same model (no core files)")

function d2_compare(res, dm)
    rows = [r for r in read_csv_rows(REF_CSV) if r["case"] == "no_pv"]
    byname = Dict(lowercase(get(b, "name", "")) => k for (k, b) in dm["bus"])
    dvmax, at, dnev, n = 0.0, "", 0.0, 0
    for r in rows
        nm = lowercase(r["bus"]); haskey(byname, nm) || continue
        k = byname[nm]; bus = dm["bus"][k]
        s = get(res["solution"]["bus"], k, nothing)
        isnothing(s) && continue
        terms = _sol_terms(bus, length(s["vr"]))
        isnothing(terms) && return (NaN, "solution vr has $(length(s["vr"])) entries, bus terminals $(bus["terminals"])", NaN, 0)
        vb = bus["vbase"] * 1000
        vc = Dict(t => complex(s["vr"][j], s["vi"][j]) for (j, t) in enumerate(terms))
        vn = get(vc, 4, 0.0im)
        for (p, col) in zip(1:3, ("Van_V", "Vbn_V", "Vcn_V"))
            d = abs(abs(vc[p] - vn) * vb - parse(Float64, r[col]))
            d > dvmax && ((dvmax, at) = (d, "$nm ph$p"))
        end
        dnev = max(dnev, abs(abs(vn) * vb - parse(Float64, r["Vng_V"])))
        n += 1
    end
    return (dvmax, at, dnev, n)
end

res1, st1 = nothing, :not_run
for (label, f) in (("solve_mc_pf", PMD.solve_mc_pf), ("solve_mc_opf", PMD.solve_mc_opf))
    global res1, st1
    dm_x = base_network(false)
    record_and_maybe_relax_thermal!(dm_x; relax=true)
    res_x = try
        f(dm_x, PMD.IVRENPowerModel, ipopt_solver)
    catch e
        println("  $label: ERROR ", first(sprint(showerror, e), 200)); nothing
    end
    isnothing(res_x) && continue
    st_x = res_x["termination_status"]
    println("  $label: ", st_x, "   (", round(res_x["solve_time"], digits=2), " s)")
    if solved(st_x)
        dv_x, at_x, dnev_x, n_x = d2_compare(res_x, dm_x)
        println("    vs OpenDSS ($n_x buses): max |ΔV_pn| = $(round(dv_x, digits=3)) V at $at_x;  max |ΔNEV| = $(round(dnev_x, digits=3)) V")
        if !solved(st1)
            res1, st1 = res_x, st_x
        end
    end
end

# --- core-file builder (same steps as build_and_solve_pv, bounds off) -------------------
function build_core!(dm; solver)
    pm = PMD.instantiate_mc_model(dm, PMD.IVRENPowerModel, (pm) -> nothing;
        jump_model = JuMP.Model(solver), ref_extensions = Function[])
    global ref = pm.ref[:it][:pmd][:nw][0]
    global model = pm.model
    global VUF_CAP = Inf
    _core_variables!()
    return nothing
end

function finish_core!()
    _core_constraints!()
    add_pv_control!()
    global objective = "cost"
    global objective_target_gens = Int[]
    _core_objectives!()
    return nothing
end

function worst_violations(; n=12)
    try
        rep = JuMP.primal_feasibility_report(model; atol=1e-6)
        viol = sort(collect(rep), by = x -> -x[2])
        println("  constraints violated by more than 1e-6: ", length(viol))
        for (c, v) in first(viol, n)
            println("   ", rpad(string(round(v, sigdigits=3)), 11), first(sprint(show, c), 170))
        end
    catch e
        println("  primal_feasibility_report failed: ", first(sprint(showerror, e), 300))
    end
end

function final_point_summary()
    JuMP.has_values(model) || return println("  no primal values available")
    v = vpn_extremes_V()
    P, Q = source_power_kw(SBASE_KVA)
    @printf("  at Ipopt's final point: V_pn %.2f V (%s) to %.2f V (%s); source %.3f kW / %.3f kvar\n",
            v.vmin, v.vmin_at, v.vmax, v.vmax_at, P, Q)
end

# --- D3 ------------------------------------------------------------------------------
banner("D3 -- core files, every load set to zero")
dm3 = base_network(false)
record_and_maybe_relax_thermal!(dm3; relax=true)
for (_, l) in dm3["load"]
    l["pd"] = zero.(l["pd"]); l["qd"] = zero.(l["qd"])
end
build_core!(dm3; solver=ipopt_solver)
finish_core!()
JuMP.optimize!(model)
st3 = JuMP.termination_status(model)
println("  status: ", st3)
final_point_summary()

# --- D4 ------------------------------------------------------------------------------
banner("D4 -- core files as in stage 3, Ipopt log on, worst violations at the final point")
loud = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 5, "max_iter" => 3000)
dm4 = base_network(false)
record_and_maybe_relax_thermal!(dm4; relax=true)
build_core!(dm4; solver=loud)
finish_core!()
JuMP.optimize!(model)
st4 = JuMP.termination_status(model)
println("\n  status: ", st4)
worst_violations()
final_point_summary()

# --- D5 ------------------------------------------------------------------------------
banner("D5 -- core files warm-started from D2's PMD solution")

function warm_start_from!(res)
    nset = 0
    for (i, bus) in ref[:bus]
        s = get(res["solution"]["bus"], string(i), nothing)
        isnothing(s) && continue
        terms = _sol_terms(bus, length(s["vr"]))
        isnothing(terms) && continue
        for (j, t) in enumerate(terms), (arr, key) in ((vr, "vr"), (vi, "vi"))
            e = arr[t, i]
            (e isa JuMP.AffExpr && !isempty(e.terms)) || continue
            JuMP.set_start_value(first(keys(e.terms)), s[key][j])
            nset += 1
        end
    end
    return nset
end

if solved(st1)
    dm5 = base_network(false)
    record_and_maybe_relax_thermal!(dm5; relax=true)
    build_core!(dm5; solver=ipopt_solver)
    println("  warm-started $(warm_start_from!(res1)) voltage variables from D2")
    finish_core!()
    JuMP.optimize!(model)
    st5 = JuMP.termination_status(model)
    println("  status: ", st5)
    if solved(st5)
        c = compare_to_reference("no_pv")
        println("  vs OpenDSS: max |ΔV_pn| = $(round(c.dv_max, digits=3)) V at $(c.dv_at);  max |ΔNEV| = $(round(c.dnev_max, digits=3)) V")
    else
        worst_violations()
        final_point_summary()
    end
else
    println("  skipped: PMD's own model did not solve either, so the network model is the place to look (D0, D1)")
end
println("\n diagnostics done")
