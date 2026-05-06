#=
    File containing the network test environments for use throughout testing.
=# 
case1_eg4w = parse_file("../test/data/opendss/ENWLNW9F6/Master.dss", transformations=[transform_loops!,remove_all_bounds!])
case1_eg4w["settings"]["sbase_default"] = 1
case1_math4w = transform_data_model(case1_eg4w, kron_reduce=false, phase_project=false)
add_start_vrvi!(case1_math4w)

function initialise_case2()
    case2_file = "../test/data/opendss/ENWLNW9F6/Master.dss"
    case2_eng4w = parse_file(case2_file, transformations=[transform_loops!,remove_all_bounds!])
    case2_eng4w["settings"]["sbase_default"] = 1
    reduce_line_series!(case2_eng4w)
    case2_math4w = transform_data_model(case2_eng4w, kron_reduce=false, phase_project=false)
    add_start_vrvi!(case2_math4w)
    return case2_math4w
end