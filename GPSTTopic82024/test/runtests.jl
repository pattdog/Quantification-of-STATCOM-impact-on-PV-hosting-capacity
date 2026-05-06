using Test
using GPSTTopic82024
using PowerModelsDistribution
using Ipopt

PowerModelsDistribution.silence!()

pmd_path = joinpath(dirname(pathof(PowerModelsDistribution)), "..")

ipopt = Ipopt.Optimizer

# include("./common.jl")
include("./testcases.jl")

@testset "GPSTTopic82024" begin

    include("doeobjective.jl")

    include("vvvwdoe.jl")
end