#=
==============================================================================
figures_pv.jl -- figures for the small-network PV hosting-capacity results
==============================================================================
Reads the CSVs written by Small_Network_PV.jl and draws eight figures. It
never solves: a change to a figure never needs a re-run of the study.

    pv_hc.csv            Fig 1  HC ladder (HC0 vs HC1)
                         Fig 2  HC2 by device mode
                         Fig 3  limit-aware vs minimum-grid-import control
    pv_points.csv        Fig 4  delivered vs curtailed at the HC point
    pv_substitution.csv  Fig 5  substitution (E7)
    pv_sweep.csv         Fig 6  worst voltage and VUF against PV size
                         Fig 7  heavy-phase magnitude, angle shift, |V2|
                         Fig 8  R/X test
Figs 6-8 need Block 4 of Small_Network_PV.jl; they are skipped with a note if
pv_sweep.csv is absent.

Each figure is saved as PDF (for the thesis) and PNG (for slides).

Run with:
    julia --project=<thesis root> figures_pv.jl            # default results folder
    julia --project=<thesis root> figures_pv.jl <folder>   # any folder holding the CSVs

Titles are off by default (the thesis caption carries the message). Set
WITH_TITLES = true for slides.

Dependencies: CairoMakie only. The CSVs are read by the small parser below,
so CSV.jl / DataFrames.jl are not needed.

Every figure describes ONE configuration: source 1.0 pu, load 25% of ADMD,
MEN 10 ohm per customer, 30 kVA STATCOM at endbus. Say so in each caption.

Not executed by its author (no Julia available when it was written): it
mirrors a Python preview that was run on the same CSVs. If a Makie call
fails on your version, the error will name the line.
==============================================================================
=#

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

const THESIS_ROOT = _find_up(SCRIPT_DIR, d -> isfile(joinpath(d, "Project.toml")); outermost=true)
isnothing(THESIS_ROOT) || Pkg.activate(THESIS_ROOT)

using CairoMakie
using Printf

const WITH_TITLES = false

const INDIR  = !isempty(ARGS) ? ARGS[1] :
               joinpath(something(THESIS_ROOT, SCRIPT_DIR), "results", "small_network_pv")
const OUTDIR = joinpath(INDIR, "figures")
mkpath(OUTDIR)

# =======================================================================
# CSV reading -- RFC 4180 quoting, no dependencies
# =======================================================================
function split_csv(line::AbstractString)
    out = String[]
    buf = IOBuffer()
    chars = collect(line)
    inq = false
    k = 1
    while k <= length(chars)
        c = chars[k]
        if inq
            if c == '"'
                if k < length(chars) && chars[k+1] == '"'
                    print(buf, '"'); k += 1
                else
                    inq = false
                end
            else
                print(buf, c)
            end
        elseif c == '"'
            inq = true
        elseif c == ','
            push!(out, String(take!(buf)))
        else
            print(buf, c)
        end
        k += 1
    end
    push!(out, String(take!(buf)))
    return out
end

function read_csv(path)
    lines = readlines(path)
    hdr = split_csv(lines[1])
    rows = Dict{String,String}[]
    for l in lines[2:end]
        isempty(strip(l)) && continue
        f = split_csv(l)
        push!(rows, Dict(hdr[i] => (i <= length(f) ? f[i] : "") for i in eachindex(hdr)))
    end
    return rows
end

num(r, k) = (v = get(r, k, ""); isempty(v) ? NaN : parse(Float64, v))
sel(rows; kw...) = filter(r -> all(get(r, String(k), "") == v for (k, v) in kw), rows)

function only_row(rows; kw...)
    m = sel(rows; kw...)
    length(m) == 1 || error("expected exactly one row for $(collect(kw)), found $(length(m))")
    return m[1]
end

# =======================================================================
# Style
# =======================================================================
hexcolor(s) = RGBf(parse(Int, s[2:3], base=16) / 255, parse(Int, s[4:5], base=16) / 255,
              parse(Int, s[6:7], base=16) / 255)
fade(c, a) = RGBAf(c.r, c.g, c.b, a)

const INK, INK2, GRID = hexcolor("#0b0b0b"), hexcolor("#52514e"), hexcolor("#e4e3df")
const C_HC0, C_HC1    = hexcolor("#a3a29b"), hexcolor("#52514e")     # baselines: neutral greys
const C_LIMIT         = hexcolor("#d03b3b")

# Device modes, most capable first. One colour per device in every figure.
# The four hues were checked for colour-blind separation in this order; two
# (green, yellow) are faint on white, so every figure carries a legend.
const MODES = [("pq", "D", "P+Q", hexcolor("#2a78d6")), ("ponly", "C", "P-only", hexcolor("#eb6834")),
               ("qonly", "B", "Q-only", hexcolor("#1baf7a")), ("qcap", "Bc", "Q-capacitive", hexcolor("#eda100"))]
const C_IDEAL, C_COST = hexcolor("#2a78d6"), hexcolor("#eb6834")

# Subscripts are drawn with rich text: the default Makie font is not
# guaranteed to carry the Unicode subscript digits.
HCn(n) = rich("HC", subscript(string(n)))

const ALLOCS = ["3/3/3", "4/3/2", "5/3/1", "6/2/1", "7/1/1"]
const XLAB   = "endbus phase allocation (customers on A/B/C)"

set_theme!(Theme(
    fontsize = 11,
    Axis = (xgridvisible = false, ygridcolor = GRID, topspinevisible = false,
            rightspinevisible = false, titlealign = :left, titlesize = 12,
            xlabelcolor = INK2, ylabelcolor = INK2),
    Legend = (framevisible = false,)))

mkfig(w, h) = try
    Figure(size = (w, h))
catch
    Figure(resolution = (w, h))        # Makie < 0.20
end

ttl(s) = WITH_TITLES ? s : ""

function save_fig(fig, name)
    save(joinpath(OUTDIR, name * ".pdf"), fig)
    save(joinpath(OUTDIR, name * ".png"), fig; px_per_unit = 2)
    println("  wrote $(joinpath(OUTDIR, name)).pdf / .png")
end

bottom_legend(pos, elems, labels; nbanks = 1) =
    Legend(pos, elems, labels; orientation = :horizontal, nbanks = nbanks,
           tellheight = true, tellwidth = false)

pelem(c)            = PolyElement(color = c)
lelem(c; ls=:solid) = LineElement(color = c, linewidth = 2, linestyle = ls)

# =======================================================================
# Data
# =======================================================================
hc  = read_csv(joinpath(INDIR, "pv_hc.csv"))
pts = read_csv(joinpath(INDIR, "pv_points.csv"))
sub = read_csv(joinpath(INDIR, "pv_substitution.csv"))
const SWEEP_PATH = joinpath(INDIR, "pv_sweep.csv")
sweep = isfile(SWEEP_PATH) ? read_csv(SWEEP_PATH) : nothing

hcval(block, alloc; kw...) = num(only_row(hc; block = block, allocation = alloc, kw...), "hc_kva_per_lot")

HC0 = [hcval("1", a; pv_mode = "unity") for a in ALLOCS]
HC1 = [hcval("1", a; pv_mode = "vvw")   for a in ALLOCS]
HCBAL = HC1[1]
X = collect(1:length(ALLOCS))

# The evaluated point that DEFINES a hosting capacity: largest passing size,
# and where both routes were run at that size, the last one (the enforced
# solve, whose dispatch is the compliant one).
function at_hc(rows)
    ok = filter(r -> r["pass"] == "true", rows)
    isempty(ok) && error("no passing point")
    xh = maximum(num(r, "x_kva_per_lot") for r in ok)
    return last(filter(r -> num(r, "x_kva_per_lot") == xh, ok))
end

# =======================================================================
# Fig 1 -- HC ladder
# =======================================================================
let
    fig = mkfig(640, 380)
    ax = Axis(fig[1, 1]; xticks = (X, ALLOCS), xlabel = XLAB,
              ylabel = "hosting capacity (kVA per customer)",
              title = ttl("Hosting capacity falls with allocation imbalance; volt-var recovers part of it"))
    w = 0.36
    barplot!(ax, X .- w/2 .- 0.01, HC0; width = w, color = C_HC0)
    barplot!(ax, X .+ w/2 .+ 0.01, HC1; width = w, color = C_HC1)
    for (off, mode, vals) in ((-w/2 - 0.01, "unity", HC0), (w/2 + 0.01, "vvw", HC1))
        labs = [@sprintf("%.2f\n%s", vals[i], only_row(hc; block = "1", allocation = ALLOCS[i], pv_mode = mode)["active_at_hc"])
                for i in X]
        text!(ax, X .+ off, vals .+ 0.2; text = labs, align = (:center, :bottom), fontsize = 8, color = INK2)
    end
    cal = [(i, num(only_row(hc; block = "1", allocation = ALLOCS[i], pv_mode = "unity"), "calib_hc0")) for i in X]
    cal = filter(c -> !isnan(c[2]), cal)
    scatter!(ax, [c[1] - w/2 - 0.01 for c in cal], [c[2] for c in cal];
             marker = :diamond, color = INK, markersize = 9)
    ylims!(ax, 0, 17)
    bottom_legend(fig[2, 1],
        [pelem(C_HC0), pelem(C_HC1), MarkerElement(marker = :diamond, color = INK, markersize = 9)],
        [rich(HCn(0), "  unity PF"), rich(HCn(1), "  volt-var + volt-watt"), "independent solver"])
    save_fig(fig, "fig1_hc_ladder")
end

# =======================================================================
# Fig 2 -- HC2 by device mode (limit-aware control)
# =======================================================================
let
    fig = mkfig(700, 400)
    ax = Axis(fig[1, 1]; xticks = (X, ALLOCS), xlabel = XLAB, ylabel = rich(HCn(2), " (kVA per customer)"),
              title = ttl("With limit-aware control, real-power exchange recovers most of the lost capacity"))
    w = 0.19
    for (k, (sc, code, lab, col)) in enumerate(MODES)
        vals = [hcval("2", a; pv_mode = "vvw", statcom = sc) for a in ALLOCS]
        barplot!(ax, X .+ (k - 2.5) * (w + 0.02), vals; width = w, color = col)
    end
    seg = Point2f[]
    for i in X
        push!(seg, Point2f(i - 0.46, HC1[i])); push!(seg, Point2f(i + 0.46, HC1[i]))
    end
    linesegments!(ax, seg; color = INK, linewidth = 2.5)
    hlines!(ax, [HCBAL]; color = INK2, linestyle = :dash, linewidth = 1.2)
    ylims!(ax, 0, 16)
    bottom_legend(fig[2, 1],
        [[pelem(m[4]) for m in MODES]..., lelem(INK), lelem(INK2; ls = :dash)],
        Any[["$(m[2])  $(m[3])" for m in MODES]..., rich(HCn(1), " (no device)"),
             rich("HC", subscript("bal"), @sprintf(" (%.2f)", HCBAL))]; nbanks = 2)
    save_fig(fig, "fig2_hc2_by_mode")
end

# =======================================================================
# Fig 3 -- limit-aware vs minimum-grid-import control, one panel per device
# =======================================================================
let
    fig = mkfig(980, 360)
    axs = Axis[]
    for (k, (sc, code, lab, col)) in enumerate(MODES)
        ax = Axis(fig[1, k]; xticks = (X, ALLOCS), xticklabelrotation = pi / 4,
                  title = "$code  $lab", ylabel = k == 1 ? "hosting capacity (kVA per customer)" : "")
        ideal = [hcval("2",  a; pv_mode = "vvw", statcom = sc) for a in ALLOCS]
        cost  = [hcval("2b", a; pv_mode = "vvw", statcom = sc) for a in ALLOCS]
        lines!(ax, X, HC0; color = C_HC0, linewidth = 2, linestyle = :dash)
        lines!(ax, X, HC1; color = C_HC1, linewidth = 2)
        for (vals, c) in ((ideal, C_IDEAL), (cost, C_COST))
            lines!(ax, X, vals; color = c, linewidth = 2)
            scatter!(ax, X, vals; color = c, markersize = 8, strokecolor = :white, strokewidth = 1.2)
        end
        text!(ax, [X[end]], [cost[end]]; text = [@sprintf("%.2f", cost[end])],
              align = (:left, :center), offset = (6, 0), fontsize = 8, color = INK2)
        xlims!(ax, 0.6, 5.9); ylims!(ax, 0, 15.5)
        k > 1 && hideydecorations!(ax; grid = false)
        push!(axs, ax)
    end
    linkyaxes!(axs...)
    bottom_legend(fig[2, 1:4],
        [lelem(C_HC0; ls = :dash), lelem(C_HC1), lelem(C_IDEAL), lelem(C_COST)],
        Any[rich(HCn(0), " unity PF"), rich(HCn(1), " no device"), "limit-aware dispatch", "minimum grid import"])
    WITH_TITLES && Label(fig[0, 1:4], "Under minimum-grid-import control, reactive capability lowers hosting capacity";
                         halign = :left, fontsize = 13)
    save_fig(fig, "fig3_control_objective")
end

# =======================================================================
# Fig 4 -- delivered and curtailed PV at the HC point
# =======================================================================
let
    cfgs = [(rich(HCn(0), " delivered"), C_HC0, a -> sel(pts; block = "1", allocation = a, pv_mode = "unity")),
            (rich(HCn(1), " delivered"), C_HC1, a -> sel(pts; block = "1", allocation = a, pv_mode = "vvw")),
            (rich(HCn(2), " C delivered"), MODES[2][4], a -> sel(pts; block = "2", allocation = a, statcom = "ponly")),
            (rich(HCn(2), " D delivered"), MODES[1][4], a -> sel(pts; block = "2", allocation = a, statcom = "pq"))]
    fig = mkfig(700, 400)
    ax = Axis(fig[1, 1]; xticks = (X, ALLOCS), xlabel = XLAB,
              ylabel = "PV active power at the HC point (kW, feeder total)",
              title = ttl("At hosting capacity, volt-var curtails 1.6-5% of available PV through the capability circle"))
    w = 0.19
    for (k, (lab, col, rowsof)) in enumerate(cfgs)
        r = [at_hc(rowsof(a)) for a in ALLOCS]
        deliv = [num(q, "pv_p_kw") for q in r]
        curt  = [max(num(q, "pv_curt_kw"), 0.0) for q in r]
        xs = X .+ (k - 2.5) * (w + 0.02)
        barplot!(ax, xs, deliv; width = w, color = col)
        barplot!(ax, xs, deliv .+ curt; fillto = deliv, width = w, color = fade(col, 0.35))
        big = [i for i in X if curt[i] > 0.5]
        isempty(big) || text!(ax, xs[big], deliv[big] .+ curt[big] .+ 3;
                              text = [@sprintf("%.1f%%", 100 * curt[i] / (deliv[i] + curt[i])) for i in big],
                              align = (:left, :center), rotation = pi / 2, fontsize = 7, color = INK2)
    end
    ylims!(ax, 0, 345)
    bottom_legend(fig[2, 1],
        [[pelem(c[2]) for c in cfgs]..., pelem(fade(INK2, 0.35))],
        Any[[c[1] for c in cfgs]..., "curtailed (% of available)"]; nbanks = 2)
    save_fig(fig, "fig4_delivered_curtailed")
end

# =======================================================================
# Fig 5 -- substitution (E7): PV fleet reactive absorption at x = HC1
# =======================================================================
let
    fig = mkfig(900, 380)
    panels = [("cost", "device objective: minimum grid import"), ("VUF", "device objective: minimum VUF")]
    axs = Axis[]
    for (j, (obj, title)) in enumerate(panels)
        allocs = [a for a in ALLOCS if !isempty(sel(sub; allocation = a, objective = obj))]
        xx = collect(1:length(allocs))
        ax = Axis(fig[1, j]; xticks = (xx, allocs), xlabel = rich("endbus phase allocation, at x = ", HCn(1)),
                  ylabel = j == 1 ? "PV fleet reactive absorption (kvar)" : "", title = title)
        w = 0.16
        none = [-num(first(sel(sub; allocation = a)), "q_pv_none_kvar") for a in allocs]
        barplot!(ax, xx .- 2 * (w + 0.02), none; width = w, color = C_HC1)
        for (k, (sc, code, lab, col)) in enumerate(MODES)
            v = [-num(only_row(sub; allocation = a, scenario = code, objective = obj), "q_pv_with_kvar") for a in allocs]
            barplot!(ax, xx .+ (k - 2) * (w + 0.02), v; width = w, color = col)
        end
        j > 1 && hideydecorations!(ax; grid = false)
        push!(axs, ax)
    end
    linkyaxes!(axs...)
    bottom_legend(fig[2, 1:2], [pelem(C_HC1), [pelem(m[4]) for m in MODES]...],
                  ["no device", ["$(m[2])  $(m[3])" for m in MODES]...])
    save_fig(fig, "fig5_substitution")
end

# =======================================================================
# Figs 6-8 -- need pv_sweep.csv (Block 4 of Small_Network_PV.jl)
# =======================================================================
if isnothing(sweep)
    println("  pv_sweep.csv not found in $INDIR -- Figs 6-8 skipped (run Block 4 of Small_Network_PV.jl)")
else
    ok = filter(r -> r["status"] == "SOLVED", sweep)
    sweep_series(alloc, mode, xs) = sort(filter(r -> r["allocation"] == alloc && r["pv_mode"] == mode &&
                                               num(r, "x_scale") == xs, ok);
                                   by = r -> num(r, "x_kva_per_lot"))
    col_of(mode) = mode == "unity" ? C_HC0 : C_HC1
    sweep_allocs = [a for a in ALLOCS if a != "3/3/3" && !isempty(sweep_series(a, "vvw", 1.0))]
    xscales = sort(unique(num(r, "x_scale") for r in ok))
    isempty(sweep_allocs) && error("pv_sweep.csv has no solved unbalanced allocation to plot")
    mode_legend(pos) = bottom_legend(pos, [lelem(C_HC0), lelem(C_HC1), lelem(C_LIMIT; ls = :dash)],
                                     ["unity PF", "volt-var + volt-watt", "limit"])

    function curve!(ax, alloc, mode, xs, col_; ls = :solid, c = col_of(mode))
        s = sweep_series(alloc, mode, xs)
        isempty(s) && return
        lines!(ax, [num(r, "x_kva_per_lot") for r in s], [num(r, col_) for r in s];
               color = c, linewidth = 2, linestyle = ls)
    end

    # ── Fig 6: worst voltage and VUF against size, one row per allocation ──
    let
        fig = mkfig(820, 300 * length(sweep_allocs) + 40)
        for (i, a) in enumerate(sweep_allocs)
            axv = Axis(fig[i, 1]; ylabel = "worst phase-to-neutral voltage (V)",
                       xlabel = i == length(sweep_allocs) ? "PV size (kVA per customer)" : "",
                       title = "$a   volt-var holds voltage down...")
            axu = Axis(fig[i, 2]; ylabel = "worst VUF (%)",
                       xlabel = i == length(sweep_allocs) ? "PV size (kVA per customer)" : "",
                       title = "...but raises unbalance at the same size")
            for mode in ("unity", "vvw")
                curve!(axv, a, mode, 1.0, "vmax_v")
                curve!(axu, a, mode, 1.0, "vuf_pct")
            end
            hlines!(axv, [253.0]; color = C_LIMIT, linestyle = :dash, linewidth = 1.2)
            hlines!(axu, [2.0];   color = C_LIMIT, linestyle = :dash, linewidth = 1.2)
        end
        mode_legend(fig[length(sweep_allocs) + 1, 1:2])
        save_fig(fig, "fig6_voltage_vuf_vs_size")
    end

    # ── Fig 7: what volt-var does to the heavy phase: magnitude, angle, |V2| ──
    # Phase a is the heavy phase in every swept allocation. Magnitude is
    # phase-to-neutral (the 253 V limit's quantity); the angle shift is
    # phase-to-ground (VUF's quantity); |V2| is what VUF's numerator measures.
    let
        a = last(sweep_allocs)
        fig = mkfig(980, 340)
        specs = [("vpn_mag_a", "heavy-phase voltage magnitude (V)",           "magnitude: lower under volt-var"),
                 ("vpg_ang_a", "heavy-phase angle shift from nominal (°)",    "angle: shifted further"),
                 ("v2_v",      rich("negative-sequence voltage |V", subscript("2"), "| (V)"), "net effect on unbalance")]
        for (k, (col_, ylab, title)) in enumerate(specs)
            ax = Axis(fig[1, k]; xlabel = "PV size (kVA per customer)", ylabel = ylab, title = title)
            for mode in ("unity", "vvw")
                curve!(ax, a, mode, 1.0, col_)
            end
            k == 1 && hlines!(ax, [253.0]; color = C_LIMIT, linestyle = :dash, linewidth = 1.2)
        end
        mode_legend(fig[2, 1:3])
        WITH_TITLES && Label(fig[0, 1:3], "$a at endbus: volt-var trades voltage magnitude for phase angle";
                             halign = :left, fontsize = 13)
        save_fig(fig, "fig7_heavy_phase_phasor")
    end

    # ── Fig 8: R/X test. VUF against size at each reactance scale, and the
    #    VUF penalty of volt-var (vvw minus unity) for every scale together ──
    let
        a = last(sweep_allocs)
        n = length(xscales)
        fig = mkfig(330 * (n + 1), 340)
        axs = Axis[]
        for (k, xs) in enumerate(xscales)
            ax = Axis(fig[1, k]; xlabel = "PV size (kVA per customer)", ylabel = k == 1 ? "worst VUF (%)" : "",
                      title = @sprintf("line reactance × %.0f", xs))
            for mode in ("unity", "vvw")
                curve!(ax, a, mode, xs, "vuf_pct")
            end
            hlines!(ax, [2.0]; color = C_LIMIT, linestyle = :dash, linewidth = 1.2)
            push!(axs, ax)
        end
        linkyaxes!(axs...)
        axd = Axis(fig[1, n + 1]; xlabel = "PV size (kVA per customer)",
                   ylabel = "VUF penalty of volt-var (percentage points)", title = "volt-var minus unity PF")
        sc_cols = [C_IDEAL, C_COST, MODES[3][4], MODES[4][4]]
        elems, labs = LineElement[], String[]
        for (k, xs) in enumerate(xscales)
            u = Dict(num(r, "x_kva_per_lot") => num(r, "vuf_pct") for r in sweep_series(a, "unity", xs))
            v = sweep_series(a, "vvw", xs)
            xv = [num(r, "x_kva_per_lot") for r in v if haskey(u, num(r, "x_kva_per_lot"))]
            dv = [num(r, "vuf_pct") - u[num(r, "x_kva_per_lot")] for r in v if haskey(u, num(r, "x_kva_per_lot"))]
            c = sc_cols[min(k, length(sc_cols))]
            lines!(axd, xv, dv; color = c, linewidth = 2)
            push!(elems, lelem(c)); push!(labs, @sprintf("reactance × %.0f", xs))
        end
        hlines!(axd, [0.0]; color = INK2, linewidth = 1)
        bottom_legend(fig[2, 1:n], [lelem(C_HC0), lelem(C_HC1), lelem(C_LIMIT; ls = :dash)],
                      ["unity PF", "volt-var + volt-watt", "2% limit"])
        bottom_legend(fig[2, n + 1], elems, labs)
        WITH_TITLES && Label(fig[0, 1:(n + 1)], "$a: the unbalance penalty of volt-var depends on the feeder's R/X";
                             halign = :left, fontsize = 13)
        save_fig(fig, "fig8_rx_test")
    end
end

println("\nDone. Figures in $OUTDIR")