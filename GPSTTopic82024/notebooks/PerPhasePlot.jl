function plot_3ph_load(; ymax::Union{Float64, Nothing}=nothing, ymin::Union{Float64, Nothing}=nothing)
	
	# Bus connections
	local line_nr = OpenDSSDirect.Lines.First()
    local vmax = -Inf
    local vmin = Inf

	if ymax !== nothing || ymin !== nothing
		local vmax = ymax
		local vmin = ymin
	end

    p1 = Plots.plot(size=(600,300), xlims=[-Inf,0.4], widen=true)
	Plots.xlabel!("Distance from source (km)")
	Plots.ylabel!("Average Voltage Magnitude (V)")

	while line_nr > 0
		voltage = [[],[],[]]
		dist = []
		buses = [OpenDSSDirect.Lines.Bus1(), OpenDSSDirect.Lines.Bus2()]
		for b in buses
			OpenDSSDirect.Circuit.SetActiveBus(b)
            [append!(voltage[i], hypot.(OpenDSSDirect.Bus.Voltages())[i]) for i in 1:3]
			append!(dist, OpenDSSDirect.Bus.Distance())
		end
		Plots.plot!(dist, voltage[1], linecolor=:blue, label=false)
		Plots.plot!(dist, voltage[2], linecolor=:purple, label=false)
		Plots.plot!(dist, voltage[3], linecolor=:green, label=false)

		line_nr = OpenDSSDirect.Lines.Next()

		if ymax !== nothing || ymin !== nothing
			continue
		end
        for v in voltage
            if v === missing || v === nothing || length(v) <= 0
                continue
            end
            minV = minimum(v)
            if minV < vmin 
                vmin = minV 
            end
            maxV = maximum(v)
            if maxV > vmax
                vmax = maxV
            end
        end
	end

    Plots.ylims!(vmin, vmax)
    Plots.hline!([230 * 1.1], label="V Max p.u.")
    Plots.hline!([230 * 0.94], label="V Min p.u.")

	# Load buses
	local load_nr = OpenDSSDirect.Loads.First()
	local voltage = [[],[],[]]
	local dist = []
	while load_nr > 0
		name = OpenDSSDirect.Loads.Name()
		OpenDSSDirect.Circuit.SetActiveElement("Load.$name")
		bus1 = OpenDSSDirect.Properties.Value("bus1")
		OpenDSSDirect.Circuit.SetActiveBus(bus1)
		[append!(voltage[i], hypot.(OpenDSSDirect.Bus.Voltages())[i]) for i in 1:3]
		append!(dist, OpenDSSDirect.Bus.Distance())
		load_nr = OpenDSSDirect.Loads.Next()
	end
	labels = ["Load Buses", false, false, false]
	[Plots.scatter!(dist, voltage[i], label=labels[i], color=:red) for i in 1:length(voltage)]
	Plots.plot!(legend=true)
	return p1, vmax, vmin
end