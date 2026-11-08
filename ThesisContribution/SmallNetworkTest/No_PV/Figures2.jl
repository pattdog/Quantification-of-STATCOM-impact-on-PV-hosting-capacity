#=
==============================================================================
FIGURES FOR SMALL_NETWORK_SWEEP2  (CairoMakie)
==============================================================================
Reads the CSVs written by Small_Network_Sweep2.jl and emits figures. Never
solves anything, so a cosmetic change costs seconds rather than a re-solve.

LANGUAGE. Titles and captions are written to be read aloud. Two conventions
are followed deliberately:

  - "Neutral-to-earth voltage", never "zero sequence". Zero-sequence CURRENT
    flows in the neutral and CAUSES neutral displacement, but the displacement
    itself is not a zero-sequence voltage. Say which quantity is meant, and
    say whether it is a voltage or a current.

  - No shorthand for constraint types. What matters to a reader is that a
    device is at its limit, not which algebraic form that limit takes.

Dependencies: CairoMakie only. The CSV reader below is stdlib, matching the
writer in the sweep script.

Makie note: Figure(size=...) needs Makie >= 0.20; on older versions replace
every size= with resolution=.

Run with:
    julia --project=<thesis root> figures2.jl
==============================================================================
=#

using CairoMakie
using Printf
using Statistics

CairoMakie.activate!(type = "png")

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

const THESIS_ROOT = _find_up(@__DIR__, d -> isfile(joinpath(d, "Project.toml")); outermost=true)
isnothing(THESIS_ROOT) && error("could not find a Project.toml above $(@__DIR__)")

const INDIR  = joinpath(THESIS_ROOT, "results", "small_network_2")
const OUTDIR = joinpath(INDIR, "figures")
mkpath(OUTDIR)

# Distance from the substation, metres. Must match the line lengths declared
# in 4_Bus_Small.txt (150 + 150 + 250 m). The sweep script prints the bus
# names it found; check them against these keys.
const DIST_M = Dict("sourcebus" => 0.0, "mid1" => 150.0, "mid2" => 300.0, "endbus" => 550.0)

const VUF_LIMIT  = 2.0     # %, drawn only, never enforced in blocks 1/2/4
const VMIN_LIMIT = 0.90
const VMAX_LIMIT = 1.10
const RATING_KVA = 30.0    # device rating used outside the rating sweep

const SCN_ORDER = ["A", "B", "Bc", "C", "D"]
const SCN_NAME  = Dict("A"  => "No STATCOM",
                       "B"  => "Reactive only (both signs)",
                       "Bc" => "Reactive only (capacitive)",
                       "C"  => "Real power exchange only",
                       "D"  => "Real + reactive")
const SCN_SHORT = Dict("A" => "None", "B" => "Q only", "Bc" => "Q (cap)",
                       "C" => "P only", "D" => "P + Q")
const SCN_COLOR = Dict("A"  => :grey35,
                       "B"  => :darkorange2,
                       "Bc" => :goldenrod3,
                       "C"  => :steelblue4,
                       "D"  => :seagreen4)
const PHASE_COLOR = Dict(1 => :firebrick, 2 => :goldenrod3, 3 => :royalblue)
const PHASE_NAME  = Dict(1 => "A", 2 => "B", 3 => "C")

# -----------------------------------------------------------------------
# Minimal CSV reader. Values stay as strings and are parsed at the point of
# use, so a malformed field fails where it is read rather than silently at
# load time.
# -----------------------------------------------------------------------
# Splits one CSV line, honouring RFC 4180 quoting. A plain split on "," breaks
# the moment any field legitimately contains a comma.
function split_csv_line(ln::AbstractString)
    out = String[]; buf = IOBuffer(); inq = false; i = 1
    chars = collect(ln)
    while i <= length(chars)
        c = chars[i]
        if inq
            if c == '"'
                if i < length(chars) && chars[i+1] == '"'
                    write(buf, '"'); i += 2; continue
                end
                inq = false
            else
                write(buf, c)
            end
        elseif c == '"'
            inq = true
        elseif c == ','
            push!(out, String(take!(buf)))
        else
            write(buf, c)
        end
        i += 1
    end
    push!(out, String(take!(buf)))
    return out
end

function read_csv(path)
    isfile(path) || error("missing $path -- run Small_Network_Sweep2.jl first")
    lines = readlines(path)
    isempty(lines) && error("$path is empty")
    header = split_csv_line(lines[1])
    rows = Dict{String,String}[]
    for (k, ln) in enumerate(lines[2:end])
        isempty(strip(ln)) && continue
        f = split_csv_line(ln)
        length(f) == length(header) ||
            error("ragged row $(k+1) in $path: got $(length(f)) fields, expected $(length(header))\n$ln")
        push!(rows, Dict(h => v for (h, v) in zip(header, f)))
    end
    return rows
end

num(r, k) = (s = get(r, k, ""); isempty(s) ? NaN : parse(Float64, s))
str(r, k) = get(r, k, "")
tf(r, k)  = lowercase(str(r, k)) == "true"
where(rows, pairs...) = filter(r -> all(str(r, k) == v for (k, v) in pairs), rows)

function distance_of(name)
    haskey(DIST_M, name) && return DIST_M[name]
    @warn "bus '$name' is not in DIST_M -- plotted at 0 m. Update the dict to match the .dss lengths."
    return 0.0
end

save2(fig, name) = (save(joinpath(OUTDIR, name * ".png"), fig, px_per_unit = 3);
                    save(joinpath(OUTDIR, name * ".pdf"), fig))

# -----------------------------------------------------------------------
summary  = read_csv(joinpath(INDIR, "sweep2_summary.csv"))
profiles = read_csv(joinpath(INDIR, "sweep2_profiles.csv"))
relaxed  = filter(r -> str(r, "status") == "SOLVED" && !endswith(str(r, "objective"), "-enforced"), summary)
enforced = filter(r -> endswith(str(r, "objective"), "-enforced"), summary)

isempty(relaxed) && error("no solved rows under relaxed limits in sweep2_summary.csv")

# Allocations in sweep order, taken from the data so trimming the sweep in
# the script cannot silently produce a figure of the wrong cases.
allocs = String[]
for r in relaxed
    a = str(r, "allocation"); a in allocs || push!(allocs, a)
end
sort!(allocs, by = a -> num(first(where(relaxed, "allocation" => a)), "spread_kw"))
const WORST = allocs[end]
spread_of(a) = num(first(where(relaxed, "allocation" => a)), "spread_kw")
@info "figures use allocation $WORST as the worst case (phase spread $(round(spread_of(WORST), digits=1)) kW)"

# The headline comparison uses the unbalance-minimising objective; scenario A
# has no controllable device so its objective is immaterial.
headline(scn) = filter(r -> str(r, "scenario") == scn &&
                        (scn == "A" || str(r, "objective") == "VUF"), relaxed)

function by_alloc(rows, field)
    d = Dict(str(r, "allocation") => num(r, field) for r in rows)
    return [get(d, a, NaN) for a in allocs]
end

xs = 1:length(allocs)
xticks_alloc = (collect(xs), allocs)

# =======================================================================
# FIGURE 1 -- voltage unbalance against how the lots are split across phases.
# The headline chart: it shows the problem appearing and each device's
# response to it, against the statutory limit.
# =======================================================================
let
    fig = Figure(size = (900, 520))
    ax = Axis(fig[1, 1],
        xlabel = "Lots on each phase at the end of the feeder",
        ylabel = "Worst voltage unbalance factor  [%]",
        title  = "Voltage unbalance as lots are allocated unevenly across phases",
        subtitle = "Total demand held at 89.8 kW throughout; STATCOM rated $(RATING_KVA) kVA",
        xticks = xticks_alloc)

    hlines!(ax, [VUF_LIMIT]; color = :firebrick, linestyle = :dash, linewidth = 2.5)
    text!(ax, 1.0, VUF_LIMIT + 0.04; text = "2% limit", color = :firebrick,
          align = (:left, :bottom), fontsize = 12)

    for scn in SCN_ORDER
        rows = headline(scn); isempty(rows) && continue
        y = by_alloc(rows, "vuf_pct")
        lines!(ax, xs, y; color = SCN_COLOR[scn], linewidth = 2.5, label = SCN_NAME[scn])
        scatter!(ax, xs, y; color = SCN_COLOR[scn], markersize = 10)
    end
    axislegend(ax; position = :lt, framevisible = false)
    save2(fig, "fig1_unbalance_vs_allocation")
end

# =======================================================================
# FIGURE 2 -- voltage along the feeder at the worst allocation.
#
# Phase-to-NEUTRAL, which is what a customer's appliance is connected
# across. The third panel is the change the device makes, on an axis fine
# enough to see it -- without it, a device working at its limit looks
# identical to no device at all.
# =======================================================================
let
    fig = Figure(size = (1180, 470))
    panels = [("A", "No STATCOM"), ("D", "Real + reactive")]
    axes_ = Axis[]

    for (k, (scn, ttl)) in enumerate(panels)
        ax = Axis(fig[1, k], xlabel = "Distance from substation  [m]",
                  ylabel = k == 1 ? "Phase-to-neutral voltage  [pu]" : "", title = ttl)
        push!(axes_, ax)
        hlines!(ax, [VMIN_LIMIT, VMAX_LIMIT]; color = :grey60, linestyle = :dot, linewidth = 1.5)
        obj = scn == "A" ? "cost" : "VUF"
        for p in 1:3
            rows = where(profiles, "allocation" => WORST, "scenario" => scn,
                         "objective" => obj, "phase" => string(p))
            isempty(rows) && continue
            d = [distance_of(str(r, "bus_name")) for r in rows]
            v = [num(r, "vm_neutral_pu") for r in rows]
            o = sortperm(d)
            lines!(ax, d[o], v[o]; color = PHASE_COLOR[p], linewidth = 2.5,
                   label = "Phase $(PHASE_NAME[p])")
            scatter!(ax, d[o], v[o]; color = PHASE_COLOR[p], markersize = 10)
        end
        k == 2 && axislegend(ax; position = :lb, framevisible = false)
    end
    length(axes_) == 2 && linkyaxes!(axes_[1], axes_[2])

    ax3 = Axis(fig[1, 3], xlabel = "Distance from substation  [m]",
               ylabel = "Change in voltage  [pu]",
               title = "What the device changes",
               subtitle = "Real + reactive, minus no STATCOM")
    hlines!(ax3, [0.0]; color = :black, linewidth = 1)
    for p in 1:3
        ra = where(profiles, "allocation" => WORST, "scenario" => "A",
                   "objective" => "cost", "phase" => string(p))
        rd = where(profiles, "allocation" => WORST, "scenario" => "D",
                   "objective" => "VUF", "phase" => string(p))
        (isempty(ra) || isempty(rd)) && continue
        da = Dict(str(r, "bus_name") => num(r, "vm_neutral_pu") for r in ra)
        dd = Dict(str(r, "bus_name") => num(r, "vm_neutral_pu") for r in rd)
        names = sort(collect(keys(da)), by = distance_of)
        d = [distance_of(n) for n in names]
        v = [get(dd, n, NaN) - da[n] for n in names]
        lines!(ax3, d, v; color = PHASE_COLOR[p], linewidth = 2.5)
        scatter!(ax3, d, v; color = PHASE_COLOR[p], markersize = 10)
    end

    Label(fig[0, :], "Voltage along the feeder — $WORST lots per phase at the end bus",
          fontsize = 16, font = :bold)
    save2(fig, "fig2_voltage_profile")
end

# =======================================================================
# FIGURE 3 -- neutral-to-earth voltage along the feeder.
#
# Zero at the substation, where the neutral is solidly earthed, and rising
# with distance as unbalanced return current flows along the neutral
# conductor. It kinks at each earthed bus. This is NOT zero-sequence
# voltage; it is the displacement that zero-sequence CURRENT produces.
# =======================================================================
let
    fig = Figure(size = (880, 500))
    ax = Axis(fig[1, 1], xlabel = "Distance from substation  [m]",
        ylabel = "Neutral-to-earth voltage  [V]",
        title  = "Neutral displacement along the feeder",
        subtitle = "$WORST lots per phase at the end bus; neutral earthed at the substation and at every customer")

    for scn in SCN_ORDER
        obj = scn == "A" ? "cost" : "VUF"
        rows = where(profiles, "allocation" => WORST, "scenario" => scn,
                     "objective" => obj, "phase" => "1")   # NEV is per bus
        isempty(rows) && continue
        d = [distance_of(str(r, "bus_name")) for r in rows]
        v = [num(r, "nev_v") for r in rows]
        o = sortperm(d)
        lines!(ax, d[o], v[o]; color = SCN_COLOR[scn], linewidth = 2.5, label = SCN_NAME[scn])
        scatter!(ax, d[o], v[o]; color = SCN_COLOR[scn], markersize = 10)
    end
    axislegend(ax; position = :lt, framevisible = false)
    save2(fig, "fig3_neutral_to_earth_voltage")
end

# =======================================================================
# FIGURE 4 -- what each converter actually dispatches, phase by phase.
#
# The mechanism figure. In the real-power scenarios the per-phase values are
# individually large and sum to zero: power is MOVED from the lightly loaded
# phases to the heavily loaded one through the shared DC link, never imported
# from the network. The reactive-only panels have no real-power bars at all.
# =======================================================================
let
    devices = [s for s in SCN_ORDER if s != "A"]
    fig = Figure(size = (1180, 470))
    for (k, scn) in enumerate(devices)
        rows = where(relaxed, "allocation" => WORST, "scenario" => scn, "objective" => "VUF")
        isempty(rows) && continue
        r = first(rows)
        P = [num(r, "p_ph$(p)_kw")   for p in 1:3]
        Q = [num(r, "q_ph$(p)_kvar") for p in 1:3]
        allv = vcat(P, Q); lo, hi = minimum(allv), maximum(allv)
        pad = 0.18 * max(hi - lo, 1e-6)

        ax = Axis(fig[1, k], xlabel = "Phase",
                  ylabel = k == 1 ? "Injection into the network  [kW / kvar]" : "",
                  title = SCN_SHORT[scn], xticks = (1:3, ["A", "B", "C"]))
        barplot!(ax, (1:3) .- 0.18, P; width = 0.34, color = :steelblue3, label = "Real  [kW]")
        barplot!(ax, (1:3) .+ 0.18, Q; width = 0.34, color = :darkorange2, label = "Reactive  [kvar]")
        hlines!(ax, [0.0]; color = :black, linewidth = 1)
        text!(ax, 2.0, lo - 0.7 * pad; text = @sprintf("total real power = %+.2f kW", sum(P)),
              align = (:center, :bottom), fontsize = 12)
        ylims!(ax, lo - pad, hi + pad)
        k == length(devices) && axislegend(ax; position = :rt, framevisible = false)
    end
    Label(fig[0, :], "Converter output at $WORST lots per phase — $(RATING_KVA) kVA device",
          fontsize = 16, font = :bold)
    save2(fig, "fig4_converter_dispatch")
end

# =======================================================================
# FIGURE 5 -- network losses.
#
# The only metric here denominated in kilowatts, and therefore in dollars.
# The dashed line is the balanced case: the floor, what it costs to deliver
# 89.8 kW through this feeder however the lots are allocated. Everything
# above it is the cost of unbalance. The dotted underlay is the share
# dissipated in the neutral conductor, which a three-wire model cannot see.
# =======================================================================
let
    base_rows = where(relaxed, "allocation" => allocs[1], "scenario" => "A")
    floor_kw  = isempty(base_rows) ? NaN : num(first(base_rows), "loss_total_kw")

    fig = Figure(size = (940, 520))
    ax = Axis(fig[1, 1], xlabel = "Lots on each phase at the end of the feeder",
        ylabel = "Network losses  [kW]",
        title  = "Losses as lots are allocated unevenly across phases",
        subtitle = "Dotted lines: the share dissipated in the neutral conductor",
        xticks = xticks_alloc)

    isnan(floor_kw) || hlines!(ax, [floor_kw]; color = :grey50, linestyle = :dash, linewidth = 2)

    for scn in SCN_ORDER
        rows = headline(scn); isempty(rows) && continue
        lines!(ax, xs, by_alloc(rows, "loss_total_kw"); color = SCN_COLOR[scn],
               linewidth = 2.5, label = SCN_NAME[scn])
        scatter!(ax, xs, by_alloc(rows, "loss_total_kw"); color = SCN_COLOR[scn], markersize = 10)
        lines!(ax, xs, by_alloc(rows, "loss_neutral_kw"); color = (SCN_COLOR[scn], 0.45),
               linestyle = :dot, linewidth = 2)
    end
    axislegend(ax; position = :lt, framevisible = false)
    save2(fig, "fig5_losses")
end

# =======================================================================
# FIGURES 6 & 7 -- converter rating.
#
# Marker fill: FILLED means the converter is running at its limit, HOLLOW
# means it has capacity to spare. That is the distinction between "the
# device is too small" and "the device cannot do this" -- a device that
# stays at its limit at EVERY rating and still falls short is the second.
# =======================================================================
let
    path = joinpath(INDIR, "sweep2_rating.csv")
    if !isfile(path)
        @warn "sweep2_rating.csv not found -- skipping rating figures."
    else
        rate = filter(r -> str(r, "status") == "SOLVED", read_csv(path))
        base = filter(r -> str(r, "scenario") == "A", rate)
        base_vuf = isempty(base) ? NaN : num(first(base), "vuf_pct")
        base_in  = isempty(base) ? NaN : num(first(base), "i_neutral_max_a")
        devices  = [s for s in SCN_ORDER if s != "A"]

        function ser(scn, field)
            rows = filter(r -> str(r, "scenario") == scn, rate)
            o = sortperm([num(r, "rating_kva") for r in rows]); rows = rows[o]
            return ([num(r, "rating_kva") for r in rows],
                    [num(r, field) for r in rows],
                    [num(r, "util") for r in rows])
        end

        fig = Figure(size = (1060, 500))
        for (k, (field, ylab, ttl, baseval)) in enumerate([
                ("vuf_pct", "Worst voltage unbalance factor  [%]",
                 "Voltage unbalance", base_vuf),
                ("i_neutral_max_a", "Largest neutral current  [A]",
                 "Current returning in the neutral", base_in)])
            ax = Axis(fig[1, k], xlabel = "Converter rating  [kVA]", ylabel = ylab, title = ttl)
            isnan(baseval) || hlines!(ax, [baseval]; color = :grey50, linestyle = :dash, linewidth = 2)
            k == 1 && hlines!(ax, [VUF_LIMIT]; color = :firebrick, linestyle = :dash, linewidth = 2)
            for scn in devices
                x, y, u = ser(scn, field); isempty(x) && continue
                lines!(ax, x, y; color = SCN_COLOR[scn], linewidth = 2.5)
                scatter!(ax, x, y;
                    color = [ui >= 0.999 ? SCN_COLOR[scn] : :white for ui in u],
                    strokecolor = SCN_COLOR[scn], strokewidth = 2, markersize = 11)
            end
        end
        elems = vcat([LineElement(color = SCN_COLOR[s], linewidth = 2.5) for s in devices],
                     [LineElement(color = :grey50, linestyle = :dash, linewidth = 2),
                      LineElement(color = :firebrick, linestyle = :dash, linewidth = 2),
                      MarkerElement(color = :grey30, marker = :circle, markersize = 11),
                      MarkerElement(color = :white, strokecolor = :grey30, strokewidth = 2,
                                    marker = :circle, markersize = 11)])
        labs = vcat([SCN_NAME[s] for s in devices],
                    ["No STATCOM", "2% limit", "running at its limit", "capacity to spare"])
        Legend(fig[2, :], elems, labs; orientation = :horizontal, framevisible = false, nbanks = 2)
        Label(fig[0, :], "Effect of converter rating — $WORST lots per phase",
              fontsize = 16, font = :bold)
        save2(fig, "fig6_rating")

        f7 = Figure(size = (820, 460))
        ax7 = Axis(f7[1, 1], xlabel = "Converter rating  [kVA]",
            ylabel = "Fraction of converter capacity in use",
            title = "How hard the converter is working",
            subtitle = "A device still at 1.0 with a larger rating is limited by what it can do, not by its size")
        hlines!(ax7, [1.0]; color = :grey50, linestyle = :dash, linewidth = 2)
        for (j, scn) in enumerate(devices)
            x, y, _ = ser(scn, "util"); isempty(x) && continue
            # Small vertical offset so series pinned at 1.0 do not hide each other.
            off = (j - (length(devices) + 1) / 2) * 0.004
            lines!(ax7, x, y .+ off; color = SCN_COLOR[scn], linewidth = 2.5, label = SCN_NAME[scn])
            scatter!(ax7, x, y .+ off; color = SCN_COLOR[scn], markersize = 9)
        end
        axislegend(ax7; position = :lb, framevisible = false)
        save2(f7, "fig7_converter_loading")
    end
end

# =======================================================================
# FIGURE 8 -- limits as references versus limits as constraints.
#
# Left: the limits are relaxed in the model and drawn only for comparison,
# so every point is an operating state the network would genuinely reach.
# Right: the limits are enforced, which is what an operator bound by them
# can dispatch. A missing point on the right means no compliant dispatch
# exists for that case -- which is a result, not a gap in the data.
# =======================================================================
let
    if isempty(enforced)
        @warn "no enforced-limit rows found -- skipping fig8."
    else
        enf_ok = filter(r -> str(r, "status") == "SOLVED", enforced)
        fig = Figure(size = (1060, 500))
        for (k, (rows, ttl, sub)) in enumerate([
                (relaxed,  "Limits drawn as references",
                 "every point is a state the network would reach"),
                (enf_ok,   "Limits enforced as constraints",
                 "missing points: no compliant dispatch exists")])
            ax = Axis(fig[1, k], xlabel = "Lots on each phase at the end of the feeder",
                      ylabel = k == 1 ? "Worst voltage unbalance factor  [%]" : "",
                      title = ttl, subtitle = sub, xticks = xticks_alloc)
            hlines!(ax, [VUF_LIMIT]; color = :firebrick, linestyle = :dash, linewidth = 2.5)
            for scn in SCN_ORDER
                sel = filter(r -> str(r, "scenario") == scn, rows)
                k == 1 && (sel = filter(r -> scn == "A" || str(r, "objective") == "VUF", sel))
                isempty(sel) && continue
                d = Dict(str(r, "allocation") => num(r, "vuf_pct") for r in sel)
                keep = [i for i in xs if haskey(d, allocs[i])]
                isempty(keep) && continue
                y = [d[allocs[i]] for i in keep]
                lines!(ax, keep, y; color = SCN_COLOR[scn], linewidth = 2.5, label = SCN_NAME[scn])
                scatter!(ax, keep, y; color = SCN_COLOR[scn], markersize = 10)
            end
            k == 1 && axislegend(ax; position = :lt, framevisible = false)
        end
        Label(fig[0, :], "The same sweep read two ways", fontsize = 16, font = :bold)
        save2(fig, "fig8_limits_reference_vs_enforced")
    end
end

# =======================================================================
# FIGURE 9 -- where the device is installed.
#
# Every other figure puts the device at the end bus, colocated with the
# unbalance it is correcting: best-case siting. This prices that assumption
# before the same question is asked of the full feeder, where the right
# location is not obvious.
# =======================================================================
let
    path = joinpath(INDIR, "sweep2_siting.csv")
    if !isfile(path)
        @warn "sweep2_siting.csv not found -- skipping fig9."
    else
        sit = filter(r -> str(r, "status") == "SOLVED", read_csv(path))
        buses = String[]
        for r in sit; b = str(r, "statcom_bus"); b in buses || push!(buses, b); end
        sort!(buses, by = distance_of)
        devices = [s for s in SCN_ORDER if s != "A"]

        fig = Figure(size = (1000, 470))
        for (k, (field, ylab, ttl)) in enumerate([
                ("vuf_pct", "Worst voltage unbalance factor  [%]", "Voltage unbalance"),
                ("loss_total_kw", "Network losses  [kW]", "Losses")])
            ax = Axis(fig[1, k], xlabel = "Where the device is installed",
                      ylabel = ylab, title = ttl,
                      xticks = (1:length(buses),
                                ["$b\n($(Int(distance_of(b))) m)" for b in buses]))
            for (j, scn) in enumerate(devices)
                y = Float64[]
                for b in buses
                    rows = filter(r -> str(r, "statcom_bus") == b && str(r, "scenario") == scn, sit)
                    push!(y, isempty(rows) ? NaN : num(first(rows), field))
                end
                barplot!(ax, (1:length(buses)) .+ (j - (length(devices)+1)/2) * 0.2, y;
                         width = 0.18, color = SCN_COLOR[scn], label = SCN_NAME[scn])
            end
            k == 1 && hlines!(ax, [VUF_LIMIT]; color = :firebrick, linestyle = :dash, linewidth = 2)
            k == 2 && axislegend(ax; position = :rt, framevisible = false)
        end
        Label(fig[0, :], "Siting — $WORST lots per phase, $(RATING_KVA) kVA device",
              fontsize = 16, font = :bold)
        save2(fig, "fig9_siting")
    end
end

# =======================================================================
# FIGURE 10 -- summary table, rendered rather than printed. Terminal output
# is not a deliverable; that includes the numbers.
# =======================================================================
let
    rows = vcat([headline(s) for s in SCN_ORDER]...)
    sort!(rows, by = r -> (findfirst(==(str(r, "allocation")), allocs),
                           findfirst(==(str(r, "scenario")), SCN_ORDER)))

    headers = ["Lots/phase", "Device", "VUF [%]", "Vmin [pu]", "Neutral [A]",
               "NEV [pu]", "Loss [kW]", "Capacity used"]
    cells = [[str(r, "allocation"),
              SCN_SHORT[str(r, "scenario")],
              @sprintf("%.3f", num(r, "vuf_pct")) * (tf(r, "breaches_vuf") ? " *" : ""),
              @sprintf("%.4f", num(r, "vmin_pu")) * (tf(r, "breaches_vmin") ? " *" : ""),
              @sprintf("%.1f", num(r, "i_neutral_max_a")),
              @sprintf("%.4f", num(r, "nev_pu")),
              @sprintf("%.3f", num(r, "loss_total_kw")),
              isempty(str(r, "util")) ? "—" : @sprintf("%.3f", num(r, "util"))] for r in rows]

    n = length(cells)
    fig = Figure(size = (940, 86 + 23 * (n + 1)))
    ax = Axis(fig[1, 1]); hidedecorations!(ax); hidespines!(ax)
    xlims!(ax, 0, length(headers)); ylims!(ax, -n - 0.5, 1.2)

    for (j, h) in enumerate(headers)
        text!(ax, j - 0.5, 0.35; text = h, align = (:center, :center), font = :bold, fontsize = 12)
    end
    hlines!(ax, [0.0]; color = :black, linewidth = 1.2)
    for (i, row) in enumerate(cells)
        j0 = findfirst(==(str(rows[i], "scenario")), SCN_ORDER)
        j0 == 1 && i > 1 && hlines!(ax, [-i + 1.0]; color = :grey80, linewidth = 0.8)
        for (j, c) in enumerate(row)
            text!(ax, j - 0.5, -i + 0.5; text = c, align = (:center, :center), fontsize = 11)
        end
    end
    Label(fig[0, :], "Summary — limits drawn as references, not enforced   (*  outside the limit)",
          fontsize = 14, font = :bold)
    save2(fig, "fig10_summary_table")
end

println("\nFigures written to $OUTDIR")
println()
println("Read fig4 first. If the per-phase real powers in the P-only and P+Q panels are")
println("individually large and sum to zero, that is the mechanism the thesis claims:")
println("real power moved between phases, not imported. If they are all near zero, the")
println("device is acting through reactive power alone and the claim does not hold.")
println()
println("Then fig6. Filled markers mean the converter is at its limit. A curve that stays")
println("filled at every rating and still falls short is limited by capability, not size.")