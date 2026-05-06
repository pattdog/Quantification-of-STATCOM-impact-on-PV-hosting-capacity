using Pkg
Pkg.activate(".") # This creates a Project.toml in your folder
Pkg.add(["PowerModelsDistribution", "OpenDSSDirect", "Ipopt", "JuMP", "DataFrames"])