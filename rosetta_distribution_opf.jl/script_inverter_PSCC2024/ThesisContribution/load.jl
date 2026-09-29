#=
==============================================================================
DIAGNOSTIC: terminals/grounded at single-phase household buses
==============================================================================
Question being answered: does the hard neg-seq VUF constraint in
constraints_PBalance.jl -- which runs unconditionally for EVERY bus, using
vr[1:3,i]/vi[1:3,i] hardcoded to rows 1,2,3 -- get corrupted by zero-padding
at single-phase buses (where variables.jl fills 0.0 for any terminal not in
that bus's own terminal list)?

If a single-phase bus genuinely has bus["terminals"] == [1,4] (just its one
energized phase + neutral), then vr[2,i], vr[3,i] etc. are literal zero
placeholders, not real physics -- and the neg-seq formula would reduce to
roughly |V1|^2/9, wildly exceeding the 0.02^2 cap at every such bus. Since
the sweep has been solving cleanly, EITHER these buses carry a full
[1,2,3,4] terminal set regardless of which phase is actually loaded (no
bug), OR something else is going on that needs a closer look.

This script settles it directly: prints terminals/grounded for the specific
single-phase household buses named in loads.txt (861, 813, 817, 835, 860,
896, 898, 900, 906), plus a full histogram of terminal-set sizes across
every bus in the network.

No solve required -- data_math carries "terminals"/"grounded" directly.

Run with:
    julia --project=. check_terminals.jl
==============================================================================
=#

using Pkg
Pkg.activate("./")
import PowerModelsDistribution
const PMD = PowerModelsDistribution
PMD.silence!()

data_path = "./rosetta_distribution_opf.jl/data/ENWL_4w_Network1_Feeder1/Master.dss"

data_eng  = PMD.parse_file(data_path, transformations=[PMD.transform_loops!])
data_math = PMD.transform_data_model(
    data_eng, multinetwork=false, kron_reduce=false, phase_project=false
)

println("="^78)
println(" Single-phase household buses (from loads_v2.txt) -- terminals/grounded")
println("="^78)

# Bus numbers taken directly from the loads.txt excerpt you shared
# (Bus1=861.1.4, Bus1=813.2.4, Bus1=817.1.4, Bus1=835.3.4, Bus1=860.1.4,
#  Bus1=896.1.4, Bus1=898.1.4, Bus1=900.1.4, Bus1=906.1.4)
single_phase_bus_ids = ["861", "813", "817", "835", "860", "896", "898", "900", "906"]

for bid in single_phase_bus_ids
    if haskey(data_math["bus"], bid)
        bus = data_math["bus"][bid]
        println("  bus $bid:  terminals = $(bus["terminals"])   grounded = $(bus["grounded"])")
    else
        println("  bus $bid:  NOT FOUND in data_math[\"bus\"] (key may differ -- check indexing)")
    end
end

println()
println("="^78)
println(" Full histogram: terminal-set size across ALL buses")
println("="^78)

histogram = Dict{Int,Int}()
examples  = Dict{Int,Vector{String}}()
for (bid, bus) in data_math["bus"]
    n = length(bus["terminals"])
    histogram[n] = get(histogram, n, 0) + 1
    push!(get!(examples, n, String[]), bid)
end

for n in sort(collect(keys(histogram)))
    ex = examples[n][1:min(5, length(examples[n]))]
    println("  $(histogram[n]) buses with $n terminal(s)   e.g. $(join(ex, ", "))")
end

println()
println("="^78)
println(" A few examples FROM EACH terminal-count group, with full detail")
println("="^78)
for n in sort(collect(keys(histogram)))
    bid = examples[n][1]
    bus = data_math["bus"][bid]
    println("  bus $bid ($n terminals): terminals=$(bus["terminals"])  grounded=$(bus["grounded"])  bus_type=$(get(bus,"bus_type","n/a"))")
end

println()
println("Interpretation:")
println("- If ALL buses show 4 terminals (including the single-phase household")
println("  buses above), the neg-seq constraint is evaluating real terminals")
println("  (some possibly grounded, i.e. tied to 0V by the grounded flag rather")
println("  than by variables.jl zero-padding) -- likely NOT a bug.")
println("- If single-phase household buses show FEWER than 4 terminals (e.g. 2:")
println("  just their one live phase + neutral), the neg-seq constraint at")
println("  those buses IS being computed from variables.jl's zero-padding, not")
println("  real physics -- this would need a real fix (e.g. skip the neg-seq")
println("  constraint at buses without all 3 phases present, matching the same")
println("  filter worst_case_vuf() already uses for reporting).")