#=
==============================================================================
TOY UNBALANCE SWEEP -- FIGURES (CairoMakie)
==============================================================================
Reads the CSVs written by toy_unbalance_sweep.jl and emits figures. Never
solves anything, so re-running after a cosmetic change costs seconds.

Deliberately dependency-light: CairoMakie only. The CSV reader below is ~20
lines of stdlib rather than a CSV.jl/DataFrames.jl dependency, matching the
writer in the sweep script.

Run with:
    julia --project=. toy_figures.jl

Makie version note: `Figure(size=...)` requires Makie >= 0.20. On older
versions replace every `size=` with `resolution=`.
==============================================================================
=#

using CairoMakie
using Printf
using Statistics

CairoMakie.activate!(type = "png")

# Depth-independent anchoring: the results directory is found by searching
# upward for the OUTERMOST Project.toml (thesis_STATCOM), matching how
# toy_unbalance_sweep.jl resolves THESIS_ROOT. Moving this script deeper
# (e.g. into No_PV/) requires no edits.
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

const THESIS_ROOT = _find_up(@__DIR__, d -> isfile(joinpath(d, "Project.toml")); outermost=true)
isnothing(THESIS_ROOT) && error("could not find a Project.toml above $(@__DIR__)")

const INDIR  = joinpath(THESIS_ROOT, "results", "toy_unbalance")
const OUTDIR = joinpath(INDIR, "figures")
mkpath(OUTDIR)

# Distance from the substation, metres. Keys must match bus_name in the CSV.
# The sweep script prints the bus names it found -- check them against this.
const BUS_DISTANCE_M = Dict(
    "sourcebus" => 0.0,
    "midbus"    => 250.0,
    "loadbus"   => 400.0,
)

const PHASE_COLOR = Dict(1 => :firebrick, 2 => :goldenrod3, 3 => :royalblue)
const PHASE_NAME  = Dict(1 => "A", 2 => "B", 3 => "C")
const SCN_NAME = Dict("A" => "No STATCOM", "B" => "STATCOM Q-only", "C" => "STATCOM P-exchange")
const SCN_COLOR = Dict("A" => :grey30, "B" => :darkorange2, "C" => :seagreen)

# -----------------------------------------------------------------------
# Minimal CSV reader. Returns a Vector of Dict{String,String}; empty fields
# stay as "". Parsing is done at the point of use so a malformed numeric
# field fails where it is read, not silently at load time.
# -----------------------------------------------------------------------
function read_csv(path)
    isfile(path) || error("missing $path -- run toy_unbalance_sweep.jl first")
    lines = readlines(path)
    isempty(lines) && error("$path is empty")
    header = split(lines[1], ',')
    rows = Dict{String,String}[]
    for ln in lines[2:end]
        isempty(strip(ln)) && continue
        f = split(ln, ',', keepempty=true)
        length(f) == length(header) || error("ragged row in $path: $ln")
        push!(rows, Dict(String(h) => String(v) for (h, v) in zip(header, f)))
    end
    return rows
end

num(r, k) = (s = get(r, k, ""); isempty(s) ? NaN : parse(Float64, s))
str(r, k) = get(r, k, "")

where(rows, pairs...) = filter(r -> all(str(r, k) == v for (k, v) in pairs), rows)

function distance_of(name)
    haskey(BUS_DISTANCE_M, name) && return BUS_DISTANCE_M[name]
    @warn "bus name '$name' not in BUS_DISTANCE_M -- plotted at 0 m. Update the dict."
    return 0.0
end

summary  = read_csv(joinpath(INDIR, "toy_summary.csv"))
profiles = read_csv(joinpath(INDIR, "toy_profiles.csv"))

solved = filter(r -> str(r, "status") == "SOLVED", summary)
isempty(solved) && error("no solved rows in toy_summary.csv")

# The most unbalanced split that solved -- used for the single-operating-point
# figures. Picked from the data rather than hard-coded, so trimming SPLITS in
# the sweep script does not silently produce a figure of the wrong case.
worst_ubi   = maximum(num(r, "unbalance_index") for r in solved)
worst_split = str(first(filter(r -> num(r, "unbalance_index") == worst_ubi, solved)), "split_kw")
@info "single-point figures use split $worst_split kW (unbalance index $(round(worst_ubi, digits=3)))"

# =======================================================================
# FIGURE 1 -- VUF vs unbalance, all three scenarios. The headline trend.
# Solid = VUF objective (what the device can do), dashed = cost objective
# (what it does when only a constraint forces it).
# =======================================================================
let
    fig = Figure(size = (820, 480))
    ax = Axis(fig[1, 1],
        xlabel = "Load unbalance index  (max − min) / mean",
        ylabel = "Worst-case VUF  [%]",
        title  = "Voltage unbalance vs load unbalance",
        subtitle = "18 kW total demand held constant; STATCOM rated 6 kVAr (33% of load)")

    for scn in ["A", "B", "C"]
        for (obj, ls, alpha) in [("VUF", :solid, 1.0), ("cost", :dash, 0.45)]
            rows = where(solved, "scenario" => scn, "objective" => obj)
            isempty(rows) && continue
            ord = sortperm([num(r, "unbalance_index") for r in rows])
            x = [num(rows[i], "unbalance_index") for i in ord]
            y = [100 * num(rows[i], "vuf_neutral") for i in ord]
            lines!(ax, x, y; color = (SCN_COLOR[scn], alpha), linestyle = ls, linewidth = 2.5)
            scatter!(ax, x, y; color = (SCN_COLOR[scn], alpha), markersize = 9)
        end
    end

    elems = [LineElement(color = SCN_COLOR[s], linewidth = 2.5) for s in ["A", "B", "C"]]
    labs  = [SCN_NAME[s] for s in ["A", "B", "C"]]
    push!(elems, LineElement(color = :grey50, linestyle = :dash, linewidth = 2.5))
    push!(labs, "cost objective")
    Legend(fig[1, 2], elems, labs, framevisible = false)

    save(joinpath(OUTDIR, "fig1_vuf_vs_unbalance.png"), fig, px_per_unit = 3)
    save(joinpath(OUTDIR, "fig1_vuf_vs_unbalance.pdf"), fig)
    fig
end

# =======================================================================
# FIGURE 2 -- voltage profile along the feeder at the worst split.
# Three phases, no-STATCOM vs P-exchange. Phase-to-NEUTRAL, which is what
# a customer's appliance actually sees on a MEN network.
# =======================================================================
let
    fig = Figure(size = (960, 460))
    panels = [("A", "No STATCOM"), ("C", "STATCOM P-exchange")]

    for (k, (scn, ttl)) in enumerate(panels)
        ax = Axis(fig[1, k],
            xlabel = "Distance from substation  [m]",
            ylabel = k == 1 ? "Phase-to-neutral voltage  [pu]" : "",
            title = ttl)

        hlines!(ax, [0.90, 1.10]; color = :grey60, linestyle = :dot, linewidth = 1.5)

        obj = scn == "A" ? "cost" : "VUF"
        for p in 1:3
            rows = where(profiles, "split_kw" => worst_split, "scenario" => scn,
                         "objective" => obj, "phase" => string(p))
            isempty(rows) && continue
            d = [distance_of(str(r, "bus_name")) for r in rows]
            v = [num(r, "vm_neutral_pu") for r in rows]
            ord = sortperm(d)
            lines!(ax, d[ord], v[ord]; color = PHASE_COLOR[p], linewidth = 2.5,
                   label = "Phase $(PHASE_NAME[p])")
            scatter!(ax, d[ord], v[ord]; color = PHASE_COLOR[p], markersize = 10)
        end
        k == 2 && axislegend(ax; position = :lb, framevisible = false)
    end

    linkyaxes!(contents(fig[1, 1])[1], contents(fig[1, 2])[1])
    Label(fig[0, :], "Voltage profile at $worst_split kW split", fontsize = 16, font = :bold)

    save(joinpath(OUTDIR, "fig2_voltage_profile.png"), fig, px_per_unit = 3)
    save(joinpath(OUTDIR, "fig2_voltage_profile.pdf"), fig)
    fig
end

# =======================================================================
# FIGURE 3 -- neutral-to-earth voltage along the feeder.
# Should be ~0 everywhere when balanced and rise with distance when not.
# =======================================================================
let
    fig = Figure(size = (820, 460))
    ax = Axis(fig[1, 1],
        xlabel = "Distance from substation  [m]",
        ylabel = "Neutral-to-earth voltage  [V]",
        title  = "Neutral displacement at $worst_split kW split",
        subtitle = "10 Ω MEN bond at the load bus; midbus unbonded")

    for scn in ["A", "B", "C"]
        obj = scn == "A" ? "cost" : "VUF"
        rows = where(profiles, "split_kw" => worst_split, "scenario" => scn,
                     "objective" => obj, "phase" => "1")   # NEV is per-bus
        isempty(rows) && continue
        d = [distance_of(str(r, "bus_name")) for r in rows]
        v = [num(r, "nev_v") for r in rows]
        ord = sortperm(d)
        lines!(ax, d[ord], v[ord]; color = SCN_COLOR[scn], linewidth = 2.5, label = SCN_NAME[scn])
        scatter!(ax, d[ord], v[ord]; color = SCN_COLOR[scn], markersize = 10)
    end
    axislegend(ax; position = :lt, framevisible = false)

    save(joinpath(OUTDIR, "fig3_nev_profile.png"), fig, px_per_unit = 3)
    save(joinpath(OUTDIR, "fig3_nev_profile.pdf"), fig)
    fig
end

# =======================================================================
# FIGURE 4 -- STATCOM per-phase dispatch. THE MECHANISM FIGURE.
# The point a reader must take away: in scenario C the per-phase real powers
# are large and sum to ~zero. No net energy is imported; power is moved
# between phases. Scenario B has no P bars at all.
# =======================================================================
let
    fig = Figure(size = (900, 460))

    for (k, scn) in enumerate(["B", "C"])
        rows = where(solved, "split_kw" => worst_split, "scenario" => scn, "objective" => "VUF")
        isempty(rows) && continue
        r = first(rows)
        P = [num(r, "statcom_p_ph$(p)_kw")   for p in 1:3]
        Q = [num(r, "statcom_q_ph$(p)_kvar") for p in 1:3]

        ax = Axis(fig[1, k],
            xlabel = "Phase",
            ylabel = k == 1 ? "Injection  [kW / kVAr]" : "",
            title  = SCN_NAME[scn],
            xticks = (1:3, ["A", "B", "C"]))

        barplot!(ax, (1:3) .- 0.18, P; width = 0.34, color = :steelblue, label = "P  [kW]")
        barplot!(ax, (1:3) .+ 0.18, Q; width = 0.34, color = :darkorange2, label = "Q  [kVAr]")
        hlines!(ax, [0.0]; color = :black, linewidth = 1)

        text!(ax, 0.5, 0.02; text = @sprintf("ΣP = %+.3f kW", sum(P)),
              space = :relative, align = (:center, :bottom), fontsize = 14)
        k == 2 && axislegend(ax; position = :rt, framevisible = false)
    end

    Label(fig[0, :], "STATCOM per-phase dispatch at $worst_split kW split",
          fontsize = 16, font = :bold)

    save(joinpath(OUTDIR, "fig4_statcom_dispatch.png"), fig, px_per_unit = 3)
    save(joinpath(OUTDIR, "fig4_statcom_dispatch.pdf"), fig)
    fig
end

# =======================================================================
# FIGURE 5 -- summary table, rendered rather than printed. Fred asked for
# figures, not terminal output; that includes the numbers.
# =======================================================================
let
    rows = filter(r -> str(r, "objective") == "VUF" || str(r, "scenario") == "A", solved)
    sort!(rows, by = r -> (num(r, "unbalance_index"), str(r, "scenario")))

    headers = ["Split [kW]", "Scn", "VUF [%]", "Vmin [pu]", "Vmax [pu]", "NEV [pu]", "ΣP [kW]"]
    cells = [[
        str(r, "split_kw"),
        str(r, "scenario"),
        @sprintf("%.3f", 100 * num(r, "vuf_neutral")),
        @sprintf("%.4f", num(r, "vmin_pu")),
        @sprintf("%.4f", num(r, "vmax_pu")),
        @sprintf("%.4f", num(r, "nev_pu")),
        isempty(str(r, "statcom_p_sum_kw")) ? "—" : @sprintf("%+.3f", num(r, "statcom_p_sum_kw")),
    ] for r in rows]

    nrow = length(cells)
    fig = Figure(size = (860, 60 + 26 * (nrow + 1)))
    ax = Axis(fig[1, 1])
    hidedecorations!(ax); hidespines!(ax)
    xlims!(ax, 0, length(headers)); ylims!(ax, -nrow - 0.5, 1.0)

    for (j, h) in enumerate(headers)
        text!(ax, j - 0.5, 0.3; text = h, align = (:center, :center), font = :bold, fontsize = 13)
    end
    hlines!(ax, [0.0]; color = :black, linewidth = 1.2)

    for (i, row) in enumerate(cells), (j, c) in enumerate(row)
        text!(ax, j - 0.5, -i + 0.5; text = c, align = (:center, :center), fontsize = 12)
    end

    Label(fig[0, :], "Summary — VUF objective (scenario A: no controllable device)",
          fontsize = 15, font = :bold)

    save(joinpath(OUTDIR, "fig5_summary_table.png"), fig, px_per_unit = 3)
    save(joinpath(OUTDIR, "fig5_summary_table.pdf"), fig)
    fig
end

# =======================================================================
# FIGURES 6 & 7 -- RATING SWEEP at the most unbalanced split.
#
# The point of these two panels together: Q-only is not merely weaker than
# P-exchange, it is categorically unable to finish the job. Its VUF floors
# out while pinned at util = 1.0, and it degrades neutral displacement
# monotonically as rating grows. P-exchange drives VUF to ~0 and comes off
# the capability limit once it has enough headroom.
#
# Marker fill encodes saturation: FILLED = at the converter's limit
# (util >= 0.999), HOLLOW = headroom remaining. That distinction is what
# separates "device too small" from "device cannot".
#
# The two scenarios hit different limits, verified in constraints_unified.jl:
#   B (Q-only)     -- no s_rated field, so it takes the ORIGINAL box path:
#                     qmin <= qg[p] <= qmax per phase, plus a sum bound, with
#                     pg pinned to zero by pmin=pmax=0. A linear box on Q.
#   C (P-exchange) -- (pg/s_rated)^2 + (qg/s_rated)^2 <= 1 per leg, plus
#                     sum(pg) == p_loss. A quadratic capability circle.
# util = sqrt(P^2+Q^2)/S_rated is the right metric for C and coincides with
# |Q|/qmax for B only because P is pinned there.
#
# The 0.999 threshold is arbitrary: Ipopt solves to ~1e-3, so anything above
# it is on the boundary within solver tolerance rather than provably at it.
# =======================================================================
let
    ratepath = joinpath(INDIR, "toy_rating.csv")
    if !isfile(ratepath)
        @warn "toy_rating.csv not found -- skipping rating figures. Re-run the sweep."
    else
        rate = filter(r -> str(r, "status") == "SOLVED", read_csv(ratepath))

        baseline = filter(r -> str(r, "scenario") == "A", rate)
        base_vuf = isempty(baseline) ? NaN : 100 * num(first(baseline), "vuf_neutral")
        base_nev = isempty(baseline) ? NaN : num(first(baseline), "nev_pu")

        function series(scn, field, scale)
            rows = filter(r -> str(r, "scenario") == scn, rate)
            ord  = sortperm([num(r, "rating_kvar") for r in rows])
            rows = rows[ord]
            x    = [num(r, "rating_kvar") for r in rows]
            y    = [scale * num(r, field) for r in rows]
            u    = [num(r, "worst_leg_utilisation") for r in rows]
            return x, y, u
        end

        # Rating at which P-exchange first comes off the capability limit.
        cx, _, cu = series("C", "vuf_neutral", 100)
        unpin_idx = findfirst(u -> u < 0.999, cu)
        unpin     = isnothing(unpin_idx) ? nothing : cx[unpin_idx]

        fig = Figure(size = (1000, 470))

        for (k, (field, scale, ylab, ttl, baseval)) in enumerate([
                ("vuf_neutral", 100.0, "Worst-case VUF  [%]",
                 "Negative sequence", base_vuf),
                ("nev_pu", 1.0, "Neutral-to-earth voltage  [pu]",
                 "Zero sequence", base_nev)])

            ax = Axis(fig[1, k], xlabel = "STATCOM rating  [kVAr]",
                      ylabel = ylab, title = ttl)

            isnan(baseval) || hlines!(ax, [baseval]; color = :grey40,
                                      linestyle = :dash, linewidth = 2)

            for scn in ["B", "C"]
                x, y, u = series(scn, field, scale)
                isempty(x) && continue
                lines!(ax, x, y; color = SCN_COLOR[scn], linewidth = 2.5)
                scatter!(ax, x, y;
                    color       = [ui >= 0.999 ? SCN_COLOR[scn] : (:white) for ui in u],
                    strokecolor = SCN_COLOR[scn], strokewidth = 2, markersize = 11)
            end

            if !isnothing(unpin)
                vlines!(ax, [unpin]; color = SCN_COLOR["C"], linestyle = :dot, linewidth = 1.5)
                if k == 1
                    # Data coordinates, not relative. Makie has no per-axis
                    # relative space (:relative applies to BOTH x and y), so
                    # anchoring an annotation to a data-space x and a relative
                    # y needs the y computed from the data.
                    allvals = vcat(series("B", field, scale)[2],
                                   series("C", field, scale)[2])
                    isnan(baseval) || push!(allvals, baseval)
                    ylo, yhi = extrema(allvals)
                    text!(ax, unpin, yhi - 0.04 * (yhi - ylo);
                        text = "P-exchange leaves\nthe capability limit",
                        align = (:left, :top), fontsize = 11,
                        offset = (6, 0), color = SCN_COLOR["C"])
                end
            end
        end

        elems = [LineElement(color = SCN_COLOR["B"], linewidth = 2.5),
                 LineElement(color = SCN_COLOR["C"], linewidth = 2.5),
                 LineElement(color = :grey40, linestyle = :dash, linewidth = 2),
                 MarkerElement(color = :grey30, marker = :circle, markersize = 11),
                 MarkerElement(color = :white, strokecolor = :grey30, strokewidth = 2,
                               marker = :circle, markersize = 11)]
        # Q-only and P-exchange hit DIFFERENT limits, so the shared marker
        # fill is labelled generically and each series names its own.
        # B: box bound  qmin <= qg <= qmax  (no s_rated set, original path).
        # C: capability circle  (P/S)^2 + (Q/S)^2 <= 1  per leg.
        labs = ["STATCOM Q-only — Q box limit",
                "STATCOM P-exchange — S² capability circle",
                "No STATCOM",
                "at converter limit", "headroom remaining"]
        Legend(fig[2, :], elems, labs; orientation = :horizontal,
               framevisible = false, nbanks = 1)

        Label(fig[0, :], "Effect of converter rating at the $worst_split kW split",
              fontsize = 16, font = :bold)

        save(joinpath(OUTDIR, "fig6_rating_sweep.png"), fig, px_per_unit = 3)
        save(joinpath(OUTDIR, "fig6_rating_sweep.pdf"), fig)

        # ---- FIGURE 7 -- utilisation, the evidence behind the marker fill ----
        f7 = Figure(size = (760, 420))
        ax7 = Axis(f7[1, 1], xlabel = "STATCOM rating  [kVAr]",
            ylabel = "Worst leg utilisation  S / S_rated",
            title = "Converter loading vs rating",
            subtitle = "Q-only pinned at its Q box; P-exchange leaves its S² circle once sufficient")
        hlines!(ax7, [1.0]; color = :grey40, linestyle = :dash, linewidth = 2)
        for scn in ["B", "C"]
            x, y, _ = series(scn, "worst_leg_utilisation", 1.0)
            isempty(x) && continue
            lines!(ax7, x, y; color = SCN_COLOR[scn], linewidth = 2.5, label = SCN_NAME[scn])
            scatter!(ax7, x, y; color = SCN_COLOR[scn], markersize = 10)
        end
        ylims!(ax7, 0.6, 1.06)
        axislegend(ax7; position = :lb, framevisible = false)
        save(joinpath(OUTDIR, "fig7_utilisation.png"), f7, px_per_unit = 3)
        save(joinpath(OUTDIR, "fig7_utilisation.pdf"), f7)
    end
end

# =======================================================================
# FIGURE 8 -- LOAD-SCALE STRESS. The headline compliance figure.
#
# The 2% line is a REFERENCE, never a constraint: VUF_CAP_PU = Inf in the
# sweep, so every point here is a freely-returned operating point. A curve
# crossing back below the line is a result, not something the solver was
# ordered to achieve.
# =======================================================================
let
    spath = joinpath(INDIR, "toy_loadscale.csv")
    if !isfile(spath)
        @warn "toy_loadscale.csv not found -- skipping. Re-run the sweep."
    else
        sc = filter(r -> str(r, "status") == "SOLVED", read_csv(spath))

        fig = Figure(size = (900, 500))
        ax = Axis(fig[1, 1], xlabel = "Total demand  [kW]",
            ylabel = "Worst-case VUF  [%]",
            title  = "Compliance under increasing demand",
            subtitle = "Limits relaxed in the model; 2% shown for reference only")

        hlines!(ax, [2.0]; color = :firebrick, linestyle = :dash, linewidth = 2.5)

        for scn in ["A", "B", "C"]
            rows = filter(r -> str(r, "scenario") == scn, sc)
            isempty(rows) && continue
            ord = sortperm([num(r, "total_kw") for r in rows])
            x = [num(rows[i], "total_kw") for i in ord]
            y = [100 * num(rows[i], "vuf_neutral") for i in ord]
            lines!(ax, x, y; color = SCN_COLOR[scn], linewidth = 2.5, label = SCN_NAME[scn])
            scatter!(ax, x, y; color = SCN_COLOR[scn], markersize = 10)
        end

        text!(ax, 0.02, 2.05; text = "2% statutory limit", space = :relative,
              align = (:left, :bottom), color = :firebrick, fontsize = 12)
        axislegend(ax; position = :lt, framevisible = false)

        save(joinpath(OUTDIR, "fig8_loadscale_compliance.png"), fig, px_per_unit = 3)
        save(joinpath(OUTDIR, "fig8_loadscale_compliance.pdf"), fig)
    end
end

# =======================================================================
# FIGURE 9 -- LOSSES. The result that converts to dollars.
#
# Left: total network loss vs converter rating. Q-only rises monotonically,
# P-exchange falls monotonically. Same device size, opposite sign.
#
# Right: the trade-off directly -- loss against the VUF actually achieved.
# Q-only tracks up and to the left (buys unbalance correction by spending
# loss); P-exchange tracks down and to the left (improves both). Reading
# right-to-left along each curve is reading increasing converter rating.
#
# Mechanism: unbalanced REAL currents create neutral current, and neutral
# I^2R is real loss that a 3-wire Kron-reduced model cannot see. P-exchange
# equalises the real currents so neutral current shrinks. Q-only adds
# circulating reactive current on top of an unchanged neutral current.
# =======================================================================
let
    ratepath = joinpath(INDIR, "toy_rating.csv")
    if !isfile(ratepath)
        @warn "toy_rating.csv not found -- skipping loss figure."
    else
        rate = filter(r -> str(r, "status") == "SOLVED", read_csv(ratepath))
        base = filter(r -> str(r, "scenario") == "A", rate)
        base_loss = isempty(base) ? NaN : num(first(base), "loss_total_kw")
        base_vuf  = isempty(base) ? NaN : 100 * num(first(base), "vuf_neutral")

        function ser(scn)
            rows = filter(r -> str(r, "scenario") == scn, rate)
            ord  = sortperm([num(r, "rating_kvar") for r in rows])
            rows = rows[ord]
            return ([num(r, "rating_kvar")        for r in rows],
                    [num(r, "loss_total_kw")      for r in rows],
                    [100 * num(r, "vuf_neutral")  for r in rows],
                    [num(r, "loss_neutral_kw")    for r in rows])
        end

        fig = Figure(size = (1000, 470))

        ax1 = Axis(fig[1, 1], xlabel = "STATCOM rating  [kVAr]",
                   ylabel = "Total network loss  [kW]",
                   title = "Loss vs converter rating")
        isnan(base_loss) || hlines!(ax1, [base_loss]; color = :grey40,
                                    linestyle = :dash, linewidth = 2)
        for scn in ["B", "C"]
            x, L, _, Ln = ser(scn)
            isempty(x) && continue
            lines!(ax1, x, L; color = SCN_COLOR[scn], linewidth = 2.5)
            scatter!(ax1, x, L; color = SCN_COLOR[scn], markersize = 10)
            lines!(ax1, x, Ln; color = (SCN_COLOR[scn], 0.45),
                   linestyle = :dot, linewidth = 2)
        end

        ax2 = Axis(fig[1, 2], xlabel = "Worst-case VUF achieved  [%]",
                   ylabel = "Total network loss  [kW]",
                   title = "The trade-off",
                   subtitle = "right to left = increasing rating")
        if !isnan(base_loss)
            scatter!(ax2, [base_vuf], [base_loss]; color = :grey40,
                     marker = :diamond, markersize = 15)
            text!(ax2, base_vuf, base_loss; text = " no STATCOM",
                  align = (:left, :center), fontsize = 11, color = :grey40)
        end
        for scn in ["B", "C"]
            _, L, V, _ = ser(scn)
            isempty(L) && continue
            lines!(ax2, V, L; color = SCN_COLOR[scn], linewidth = 2.5)
            scatter!(ax2, V, L; color = SCN_COLOR[scn], markersize = 10)
        end

        elems = [LineElement(color = SCN_COLOR["B"], linewidth = 2.5),
                 LineElement(color = SCN_COLOR["C"], linewidth = 2.5),
                 LineElement(color = :grey40, linestyle = :dash, linewidth = 2),
                 LineElement(color = :grey30, linestyle = :dot, linewidth = 2)]
        labs = ["STATCOM Q-only", "STATCOM P-exchange", "No STATCOM",
                "neutral-conductor share (indicative)"]
        Legend(fig[2, :], elems, labs; orientation = :horizontal, framevisible = false)

        Label(fig[0, :], "Network losses at the $worst_split kW split",
              fontsize = 16, font = :bold)

        save(joinpath(OUTDIR, "fig9_losses.png"), fig, px_per_unit = 3)
        save(joinpath(OUTDIR, "fig9_losses.pdf"), fig)
    end
end

println("\nFigures written to $OUTDIR")
println("Read fig4 first: if ΣP is not ≈ 0 in scenario C, the power-exchange")
println("constraint is not doing what the thesis claims it does.")