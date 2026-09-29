using Pkg
Pkg.activate("./")
using rosetta_distribution_opf
import PowerModelsDistribution
import InfrastructureModels
using Ipopt
using JuMP
import LinearAlgebra: diag, diagm
const PMD = PowerModelsDistribution
const RPMD = rosetta_distribution_opf
const IM = InfrastructureModels

## Data Loading and Pre-processing
data_path = "./data/ENWL_4w_Network1_Feeder1/Master.dss"
ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes","warm_start_init_point"=>"yes")

data_eng = PMD.parse_file(data_path, transformations=[PMD.transform_loops!])
data_eng["settings"]["sbase_default"] = 1

# Keep generator impedance at 0
data_eng["voltage_source"]["source"]["rs"] *= 0
data_eng["voltage_source"]["source"]["xs"] *= 0 

data_math = PMD.transform_data_model(data_eng, multinetwork=false, kron_reduce=false, phase_project=false)

# Set Voltage Bounds
for (i, bus) in data_math["bus"]
    bus["vmin"] = [0.9 * ones(3) ; 0 ]
    bus["vmax"] = [1.1 * ones(3) ; Inf]
end

# Set Load Scaling
for (i, load) in data_math["load"]
    load["pd"] *= 1
    load["qd"] *= 3
end


# Slack Bus Cost
data_math["gen"]["1"]["cost"] = [1000 0]

## Add STATCOMs to the Network
# Every 10th load gets a STATCOM
for i = 1:10:length(data_math["load"])
    load = data_math["load"]["$i"]
    pd = load["pd"][1]

    gen_id = length(data_math["gen"]) + 1
    data_math["gen"]["$gen_id"] = deepcopy(data_math["gen"]["1"])
    
    Smax = 10 * ceil(pd) * ones(3)
    gen = data_math["gen"]["$gen_id"]
    gen["gen_bus"] = copy(load["load_bus"])
    gen["pmax"] = 0.0 * Smax
    gen["pmin"] = 0.0 * Smax
    gen["qmax"] = sqrt.(Smax.^2 - gen["pmax"].^2)
    gen["qmin"] = -gen["qmax"]
    gen["cost"] = [10 0]
    gen["type"] = "Gen" # Label for identification
end

ref = IM.build_ref(data_math, PMD.ref_add_core!, PMD._pmd_global_keys, PMD.pmd_it_name)[:it][:pmd][:nw][0]



## Model Construction
objective = "cost"
model = JuMP.Model(ipopt_solver)
include("./core/variables.jl")
include("./core/constraints.jl")
include("./core/objectives.jl")

## Optimization
JuMP.optimize!(model)
status = JuMP.termination_status(model)
@assert status in [LOCALLY_SOLVED, ALMOST_LOCALLY_SOLVED] "Solver failed with status: $status"

obj_val_GFM = JuMP.objective_value(model)
solve_time_GFM = JuMP.solve_time(model)

## 1. Voltage Performance Report
v = value.(vr) .+ im * value.(vi)
v012 = T * Array(v[1:3,:]) 

v1 = abs.(v012[2,:]) # Positive sequence
v2 = abs.(v012[3,:]) # Negative sequence
vuf = (v2 ./ v1) .* 100 # Unbalance Factor

println("--- STATCOM PERFORMANCE REPORT ---")
println("Max Negative Sequence Voltage (V2): ", maximum(v2))
println("Max Voltage Unbalance Factor (VUF): ", maximum(vuf), "%")
println("Min Voltage Magnitude (Phase A-C): ", minimum(abs.(v[1:3,:])))

## 2. STATCOM Utilization Report
# Map keys to ensure index 1 is Slack and 2+ are STATCOMs
gen_keys = sort(collect(keys(data_math["gen"])), by=x->parse(Int, x))
statcom_keys = [k for k in gen_keys if get(data_math["gen"][k], "type", "") == "Gen"]

# Extract numeric dispatch and force to standard Array
qg_val_matrix = Array(value.(qg)) 

# Theoretical Capacity (Excluding Slack)
total_qmax = sum((sum(data_math["gen"][k]["qmax"]) for k in statcom_keys), init=0.0)

# Total Dispatch (Columns 2 to end are the STATCOMs)
qg_statcom_total = sum(abs.(qg_val_matrix[:, 2:end]), init=0.0) ##abs

println("--- STATCOM UTILIZATION ---")
println("Total STATCOM Support: ", round(qg_statcom_total, digits=4), " kVAR")
println("Total STATCOM Capacity: ", round(total_qmax, digits=4), " kVAR")

if total_qmax > 1e-6
    utilization = (qg_statcom_total / total_qmax) * 100
    println("STATCOM Utilization: ", round(utilization, digits=2), "%")
end

## 3. Identification of "Hardest Working" STATCOM
if size(qg_val_matrix, 2) > 1
    # Sum absolute reactive power across phases for each STATCOM column
    q_per_unit = [sum(abs.(qg_val_matrix[:, i]), init=0.0) for i in 2:size(qg_val_matrix, 2)]
    max_q, local_idx = findmax(q_per_unit)
    
    # Map back to the specific Key and Bus ID
    # gen_keys[local_idx + 1] accounts for the offset of skipping the slack bus
    top_statcom_key = gen_keys[local_idx + 1] 
    top_statcom_bus = data_math["gen"][top_statcom_key]["gen_bus"]
    
    println("Max Dispatch (Single Unit): ", round(max_q, digits=4), " kVAR")
    println("Hardest Working STATCOM ID: ", top_statcom_key)
    println("Connected to Bus: ", top_statcom_bus)
else
    println("No STATCOMs detected in network.")
end
#=
v = value.(vr) .+ im * value.(vi)
v_axes2 = v.axes[2]
v012 = T * Array(v[1:3,:])
v2 = v012[3,:]
v2m_GEN = abs.(v2)
v0m_GEN = abs.(v012[1,:])
(v2m_max, idx) = findmax(v2m_GEN)
idx_bus = v_axes2[idx]
round.(abs.(v[:,idx_bus]), digits=4)
round.(angle.(v[:,idx_bus])*180/pi, digits=2)

c = value.(cr) .+ im * value.(ci)
c_axes2 = c.axes[2]
c012 = T * Array(c[1:3,:])
c2 = c012[3,:]
c2m_GEN = abs.(c2)


[(i, gen["pmax"][1]) for (i, gen) in data_math["gen"]]
pg_vals = Array(value.(pg)[:,1:7])
qg_vals = Array(value.(qg)[:,1:7])


gen_buses = [(parse(Int,i), gen["gen_bus"]) for (i, gen) in data_math["gen"]]
vm_gens = abs.(v[:,last.(gen_buses)])
round.(vm_gens[1:3,:] .- 1, digits=3) ./ 2
=#