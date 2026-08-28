using PowerModelsDistribution
const PMD = PowerModelsDistribution
PMD.silence!()

data_path = "./rosetta_distribution_opf.jl/data/ENWL_4w_Network1_Feeder1/Master.dss"

data_eng  = PMD.parse_file(data_path, transformations=[PMD.transform_loops!])
data_math = PMD.transform_data_model(
    data_eng, multinetwork=false, kron_reduce=false, phase_project=false
)

println("sbase             : ", data_math["settings"]["sbase"])
println("power_scale_factor: ", get(data_math["settings"], "power_scale_factor", "n/a"))
println()

# one example load, in pu
first_load_id = first(keys(data_math["load"]))
load = data_math["load"][first_load_id]
println("Example load ($first_load_id):")
println("  pd (pu) = ", load["pd"])
println("  qd (pu) = ", load["qd"])
println()

# total feeder demand, in pu
total_pd = sum(sum(l["pd"]) for (i,l) in data_math["load"])
total_qd = sum(sum(l["qd"]) for (i,l) in data_math["load"])
println("Total feeder demand (pu): pd = $total_pd   qd = $total_qd")
println("Number of loads: ", length(data_math["load"]))