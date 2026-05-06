@info "Running vvwdoe objective tests"

function add_gens!(math4w)
    gen_counter = 2
    for (d, load) in math4w["load"]
        if mod(load["index"], 2) == 1
            # phases = 3
            phases = length(load["connections"])-1
            math4w["gen"]["$gen_counter"] = deepcopy(math4w["gen"]["1"])
            math4w["gen"]["$gen_counter"]["name"] = "$gen_counter"
            math4w["gen"]["$gen_counter"]["index"] = gen_counter
            math4w["gen"]["$gen_counter"]["cost"] = 1.0 #*math4w["gen"]["1"]["cost"]
            math4w["gen"]["$gen_counter"]["gen_bus"] = load["load_bus"]
            math4w["gen"]["$gen_counter"]["pmax"] = 5.0*ones(phases)
            math4w["gen"]["$gen_counter"]["pmin"] = 0.0*ones(phases)
            math4w["gen"]["$gen_counter"]["qmax"] = 5.0*ones(phases)
            math4w["gen"]["$gen_counter"]["qmin"] = -5.0*ones(phases)
            math4w["gen"]["$gen_counter"]["connections"] = load["connections"]
            gen_counter = gen_counter + 1
        end
    end
end

@testset "vvvwdoe" begin
    
    @testset "solve_mc_vvvw_doe_competitive" begin
        case2_math4w = initialise_case2()
        for (i,bus) in case2_math4w["bus"]
            if bus["bus_type"] != 3 && !startswith(bus["source_id"], "transformer")
                bus["vm_pair_lb"] = [(1, 4, 0.9);(2, 4, 0.9);(3, 4, 0.9)]
                bus["vm_pair_ub"] = [(1, 4, 1.1);(2, 4, 1.1);(3, 4, 1.1)]
                # bus["grounded"] .=  0
            else
                @show bus
            end

        end

        for (g,gen) in case2_math4w["gen"]
            gen["cost"] = 0.0
        end

        for (d,load) in case2_math4w["load"]
            load["pd"] .*= 30.0
            load["qd"] .*= 30.0
        end
        add_gens!(case2_math4w)
        case2_math4w["gen"]["5"]["pmax"]*=0.2
        res_comp = solve_mc_vvvw_doe_competitive(case2_math4w, ipopt)
        pg_cost = [gen["pg_cost"] for (g,gen) in res_comp["solution"]["gen"] if g!="1"]

        @test isapprox(res_comp["objective"], 51.00000048242378)
        @test res_comp["termination_status"] == LOCALLY_SOLVED
        @test isapprox(pg_cost, [5.000000047493472, 5.000000047493286, 5.000000047493474, 5.000000047493147, 5.000000047493147, 1.0000000074904265, 5.0000000474932955, 5.0000000474932955, 5.000000047493472, 5.000000047493473, 5.000000047493292]; atol=0.001)
    end

    @testset "solve_mc_vvvw_doe_mse" begin
        case2_math4w = initialise_case2()
        for (i,bus) in case2_math4w["bus"]
            if bus["bus_type"] != 3 && !startswith(bus["source_id"], "transformer")
                bus["vm_pair_lb"] = [(1, 4, 0.9);(2, 4, 0.9);(3, 4, 0.9)]
                bus["vm_pair_ub"] = [(1, 4, 1.1);(2, 4, 1.1);(3, 4, 1.1)]
                # bus["grounded"] .=  0
            else
                @show bus
            end

        end

        for (g,gen) in case2_math4w["gen"]
            gen["cost"] = 0.0
        end

        for (d,load) in case2_math4w["load"]
            load["pd"] .*= 30.0
            load["qd"] .*= 30.0
        end
        add_gens!(case2_math4w)
        case2_math4w["gen"]["5"]["pmax"]*=0.2
        res_ms = solve_mc_vvvw_doe_mse(case2_math4w, ipopt)
        pg_cost = [gen["pg_cost"] for (g,gen) in res_ms["solution"]["gen"] if g!="1"]

        @test isapprox(res_ms["objective"], -20.50000004749447)
        @test res_ms["termination_status"] == LOCALLY_SOLVED
        @test isapprox(pg_cost, [3.5000000020102293, 3.500000002080621, 3.500000002010128, 3.500000002129766, 3.500000002129661, 1.000000009498892, 3.5000000020768756, 3.500000002076897, 3.5000000020101267, 3.500000002010125, 3.500000002077404]; atol=0.001)
    end

    @testset "solve_mc_vvvw_doe_abs" begin
        case2_math4w = initialise_case2()
        for (i,bus) in case2_math4w["bus"]
            if bus["bus_type"] != 3 && !startswith(bus["source_id"], "transformer")
                bus["vm_pair_lb"] = [(1, 4, 0.9);(2, 4, 0.9);(3, 4, 0.9)]
                bus["vm_pair_ub"] = [(1, 4, 1.1);(2, 4, 1.1);(3, 4, 1.1)]
                # bus["grounded"] .=  0
            else
                @show bus
            end
        
        end

        for (g,gen) in case2_math4w["gen"]
            gen["cost"] = 0.0
        end
        
        for (d,load) in case2_math4w["load"]
            load["pd"] .*= 30.0
            load["qd"] .*= 30.0
        end

        add_gens!(case2_math4w)
        case2_math4w["gen"]["5"]["pmax"]*=0.2

        res_abs = solve_mc_vvvw_doe_abs(case2_math4w, ipopt)
        pg_cost = [gen["pg_cost"] for (g,gen) in res_abs["solution"]["gen"] if g!="1"]

        @test isapprox(res_abs["objective"], -36.00000042236287)
        @test res_abs["termination_status"] == LOCALLY_SOLVED
        @test isapprox(pg_cost, [5.000000046237701, 5.000000046237235, 5.000000046237705, 5.000000046236949, 5.000000046236949, 1.0000000090604235, 5.000000046237269, 5.000000046237269, 5.000000046237704, 5.000000046237704, 5.000000046237263]; atol=0.001)
    end

    @testset "solve_mc_vvvw_doe_equal" begin
        case2_math4w = initialise_case2()
        for (i,bus) in case2_math4w["bus"]
            if bus["bus_type"] != 3 && !startswith(bus["source_id"], "transformer")
                bus["vm_pair_lb"] = [(1, 4, 0.9);(2, 4, 0.9);(3, 4, 0.9)]
                bus["vm_pair_ub"] = [(1, 4, 1.1);(2, 4, 1.1);(3, 4, 1.1)]
                # bus["grounded"] .=  0
            else
                @show bus
            end
        
        end

        for (g,gen) in case2_math4w["gen"]
            gen["cost"] = 0.0
        end
        
        for (d,load) in case2_math4w["load"]
            load["pd"] .*= 30.0
            load["qd"] .*= 30.0
        end

        add_gens!(case2_math4w)
        case2_math4w["gen"]["5"]["pmax"]*=0.2

        res_eq = solve_mc_vvvw_doe_equal(case2_math4w, ipopt)
        pg_cost = [gen["pg_cost"] for (g,gen) in res_eq["solution"]["gen"] if g!="1"]

        @test isapprox(res_eq["objective"], 1.0000000074640791)
        @test res_eq["termination_status"] == LOCALLY_SOLVED
        @test isapprox(pg_cost, [1.0000000074640765, 1.0000000074640765, 1.0000000074640765, 1.0000000074640765, 1.0000000074640765, 1.0000000074640765, 1.0000000074640765, 1.0000000074640765, 1.0000000074640765, 1.0000000074640765, 1.0000000074640765]; atol=0.001)
    end
end