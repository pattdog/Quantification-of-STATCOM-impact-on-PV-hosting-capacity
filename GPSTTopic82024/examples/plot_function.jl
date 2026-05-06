function plot_3ph_load()
	local p = Plots.plot(size=(600,300), xlims=[-Inf,0.4], widen=true)
	Plots.xlabel!("Distance from source (km)")
	Plots.ylabel!("Voltage Magnitude (V)")
	# Plots.title!("Voltage at LoadBus by Distance from Source")
	# Bus connections
	local line_nr = ODSS.Lines.First()
	local vlims = []
	while line_nr > 0
		local voltage = [[],[],[]]
		local dist = []
		local buses = [ODSS.Lines.Bus1(), ODSS.Lines.Bus2()]
		for b in buses
			ODSS.Circuit.SetActiveBus(b)
			[append!(voltage[i], hypot.(ODSS.Bus.Voltages())[i]) for i in 1:3]
			append!(dist, ODSS.Bus.Distance())
		end
		Plots.plot!(dist, voltage[1], linecolor=:blue, label=false)
		Plots.plot!(dist, voltage[2], linecolor=:purple, label=false)
		Plots.plot!(dist, voltage[3], linecolor=:green, label=false)
		# Plots.plot!(dist, voltage[4], linecolor=:brown, label=false)
		[append!(vlims, v) for v in voltage]
		line_nr = ODSS.Lines.Next()
	end
	# Fix ylims
	local ymin = minimum(vlims)
	local ymax = maximum(vlims)
	if ymin < minimum(global_ymin)
		append!(global_ymin, ymin)
	end
	if ymax > maximum(global_ymax)
		append!(global_ymax, ymax)
	end
	Plots.ylims!(minimum(global_ymin), maximum(global_ymax))
	# Load buses
	local load_nr = ODSS.Loads.First()
	local voltage = [[],[],[]]
	local dist = []
	while load_nr > 0
		local name = ODSS.Loads.Name()
		ODSS.Circuit.SetActiveElement("Load.$name")
		local bus1 = ODSS.Properties.Value("bus1")
		ODSS.Circuit.SetActiveBus(bus1)
		[append!(voltage[i], hypot.(ODSS.Bus.Voltages())[i]) for i in 1:3]
		append!(dist, ODSS.Bus.Distance())
		load_nr = ODSS.Loads.Next()
	end
	local labels = ["Load Buses", false, false, false]
	[Plots.scatter!(dist, voltage[i], label=labels[i], color=:red) for i in 1:length(voltage)]
	Plots.plot!(legend=true)
	p
end