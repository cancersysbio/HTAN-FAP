include("init.jl")
include("expandpop.jl")
include("sequence_cells.jl")

#~~~~~~~~~ Imports ~~~~~~~~~~#
# Pkg.add(["FileIO", "CSV", "DataFrames", "HDF5", "Random", "Distributions", "ArgParse", "PoissonRandom", "StatsBase"])
using Pkg
using FileIO
using CSV
using DataFrames
using HDF5
using Random
using Distributions
using ArgParse

#~~~~~~~~~ Outputs ~~~~~~~~~~#
function save_model_outputs(output_dir::String, outputs::Dict{String, Any})
    # Create the directory if it doesn't exist
    try
        if !isdir(output_dir)
            mkpath(output_dir)
        end
    catch e
        println("Error creating directory: $e")
        return
    end

    # Save each output to a separate file
    for (filename, content) in outputs
        if content != "empty"
            try
                file_path = joinpath(output_dir, filename)
                open(file_path, "w") do f
                    if isa(content, AbstractVector) && length(size(content)) == 1
                        for value in content
                            write(f, string(value) * "\n")
                        end
                    elseif isa(content, DataFrame)
                        CSV.write(file_path, content, delim='\t', writeheader=true, rownumber=false)
                    else
                        write(f, string(content))
                    end
                end
            catch e
                println("Error saving file $filename: $e")
            end
        end
    end
end

function write_arrays_to_file(arrays::Vector{Vector{T}}, filename::String) where T
    open(filename, "w") do file
        for array in arrays
            println(file, join(array, " "))
        end
    end
end

#~~~~~~~~~ Parameters ~~~~~~~~~~#

function main()
    s = ArgParseSettings()
    
    @add_arg_table s begin
        "--runtype"
        help = "Type of run"
        arg_type = String
        default = "full_model"
        
        "--initSize"
        help = "Initial size"
        arg_type = Int
        default = 10^4
        
        "--birth_rate"
        help = "Birth rate"
        arg_type = Float64
        default = 0.05548397814
        
        "--death_rate"
        help = "Death rate"
        arg_type = Float64
        default = 0.00548397814
        
        "--mut_rate"
        help = "Mutation rate"
        arg_type = Float64
        default = 12.0
        
        "--adv_mut_rate"
        help = "Advantageous mutation rate"
        arg_type = Float64
        default = 0.0
        
        "--s_coef"
        help = "Selection coefficient"
        arg_type = Float64
        default = 0.0
        
        "--num_seeds"
        help = "Number of seeds"
        arg_type = Int
        default = 1
        
        "--final_pop_size"
        help = "Final population size"
        arg_type = Int
        default = 10000
        
        "--s_coef_polyp"
        help = "Selection coefficient for polyp"
        arg_type = Float64
        default = 0.0
        
        "--adv_mut_rate_polyp"
        help = "Advantageous mutation rate for polyp"
        arg_type = Float64
        default = 0.0
        
        "--polyp_birth_rate"
        help = "Birth rate for polyp"
        arg_type = Float64
        default = 0.9

        "--polyp_death_rate"
        help = "Death rate for polyp"
        arg_type = Float64
        default = 0.05
        
        "--mut_rate_polyp"
        help = "Mutation rate for polyp"
        arg_type = Float64
        default = 12.0

        "--polyp_init_time"
        help = "Initiation time of polyp (also when the normal stops)"
        arg_type = Float64
        default = 5
        
        "--outfile"
        help = "Output file"
        default = "/home/rschenck/oak_dir/cluster/outputs/inferences/inferences_set1/"
        
        "--seed"
        help = "Random seed for simulations"
        arg_type = Int128
        default = 1234
        
        "--max_t"
        help = "Maximum time for the total simulation. All parameters are in 6 month ticks and adjusted within to monthly (inutero) or 6 months in other phases"
        arg_type = Float64
        default = 10

        "--verbose"
        help = "Verbose output"
        arg_type = Bool
        default = false

        "--mean_depth"
        help = "Mean sequencing depth"
        arg_type = Int
        default = 100

        "--sd_depth"
        help = "Standard deviation of sequencing depth"
        arg_type = Int
        default = 10

        "--sample_size"
        help = "Sample size of cells to take"
        arg_type = Int
        default = 10000

        "--purity"
        help = "Purity of the sample"
        arg_type = Float64
        default = 1.0

        "--runexpandmuts"
        help = "Whether to expand mutations in the polyp only to recover subclonal mutations and adjust VAF based on purity."
        arg_type = Bool
        default = false
    end
    
    parsed_args = parse_args(s)
    
    runtype = parsed_args["runtype"]
    run_normal = runtype == "full_model"
    
    initSize = parsed_args["initSize"]
    birth_rate = parsed_args["birth_rate"]
    death_rate = parsed_args["death_rate"]
    mut_rate = parsed_args["mut_rate"]
    adv_mut_rate = parsed_args["adv_mut_rate"]
    s_coef = parsed_args["s_coef"]
    
    num_seeds = parsed_args["num_seeds"]
    final_pop_size = parsed_args["final_pop_size"]
    s_coef_polyp = parsed_args["s_coef_polyp"]
    adv_mut_rate_polyp = parsed_args["adv_mut_rate_polyp"]
    polyp_birth_rate = parsed_args["polyp_birth_rate"]
    plyp_death_rate = parsed_args["polyp_death_rate"]
    mut_rate_polyp = parsed_args["mut_rate_polyp"]
    polyp_init_time = parsed_args["polyp_init_time"]
    
    outfile = parsed_args["outfile"]
    max_t = parsed_args["max_t"]

    verbose = parsed_args["verbose"]
    mean_depth = parsed_args["mean_depth"]
    sd_depth = parsed_args["sd_depth"]
    sample_size = parsed_args["sample_size"]
    purity = parsed_args["purity"]
    runexpandmuts = parsed_args["runexpandmuts"]

    Random.seed!(parsed_args["seed"])

    if runtype != "full_model"
        run_normal = false
    else
        run_normal = true
    end

    if birth_rate == 0.0
        # Exit 
        println("Birth rate cannot be 0.0")
        exit()
    end

    #~~~~~~~~~ Run Normal Population ~~~~~~~~~~#
    # Setup arrays and variables
    mutId::Int64 = 1
    cloneId::Int64 = 1
    cellArr = zeros(Int64, initSize)
    parentCellArr = zeros(Int64, initSize)
    cellGenomes = Vector{Vector{Int64}}(undef, 1)
    cellGenomes[1] = Int64[]
    adv_clones_arr = [0]
    cellArr = [0]
    parentCellArr = [1]
    cellArr[1] = cloneId
    cellGenomes[1] = [0]

    mutInductionTimes = Dict(mutId => 0.0)


    #~~~~~~~~~ STAGE 1 IN UTERO: Timescale == Months
    init_birth = log(initSize) / 9.0

    m1 = Model_BD(Int64(initSize), init_birth , death_rate / 12.0, mut_rate / 12.0, adv_mut_rate / 12.0, s_coef , cellArr, adv_clones_arr, cellGenomes, parentCellArr, mutId, cloneId, mutInductionTimes)

    cellArr, cellGenomes, adv_clones_arr, mutId, driver_muts, model, popsize1 = expandPop(m1, group="inutero", verbose=verbose, init_birth=init_birth, polyp_init_time=Int(polyp_init_time), max_time=9)

    mutId_delineator = [mutId]

    origindict = Dict{Int64, Int64}()
    for i in 1:length(cellArr)
        origindict[cellArr[i]] = 0
    end

    df1, vafs1 = build_output_table(cellGenomes, m1.mutInductionTimes, driver_muts, "inutero", origindict, mean_depth, sd_depth, sample_size=15000, mut_rate=Int(floor(mut_rate/2.)), purity=purity, samplehist=false, overdispersion=.001)
    print(length(vafs1))

   #~~~~~~~~~~ STAGE 2 NORMAL: Timescale == 6 months
    m = Model_BD(initSize, birth_rate/2., (death_rate)/2., mut_rate/2., adv_mut_rate/2., s_coef, cellArr, adv_clones_arr, cellGenomes, parentCellArr, mutId, cloneId, mutInductionTimes)

    cellArr, cellGenomes, adv_clones_arr, mutId, driver_muts_fission, model, popsize2 = expandPop(m, group="fission", verbose=verbose, polyp_init_time=Int(polyp_init_time), max_time=Int(polyp_init_time))

    push!(mutId_delineator, mutId)

    origindict = Dict{Int64, Int64}()
    for i in 1:length(cellArr)
        origindict[cellArr[i]] = 0
    end

    df2, vafs2 = build_output_table(cellGenomes, m.mutInductionTimes, driver_muts, "Normal", origindict, mean_depth, sd_depth, sample_size=15000, mut_rate=Int(floor(mut_rate/2.)), purity=purity, samplehist=false, overdispersion=.001)

    df = vcat(df1, df2)

    if verbose
        println("Polyp Growth Start")
    end

    #~~~~~~~~~ Polyp growth
    too_pick = deepcopy(adv_clones_arr)
    # add one to all of too_pick
    too_pick = too_pick .+ 1
    too_pick = too_pick ./ sum(too_pick)
    seed_idxs = sample(1:length(cellArr), Weights(too_pick), num_seeds, replace=false)
    global too_pick # Remove this from environment
    cellArr = deepcopy(cellArr[seed_idxs])
    parentCellArr = zeros(Int64, final_pop_size)
    cellGenomes = deepcopy(cellGenomes[seed_idxs])
    adv_clones_arr = deepcopy(adv_clones_arr[seed_idxs])
    parentCellArr[1:num_seeds] .= 1
    cloneId = num_seeds

    seed_genomes = deepcopy(cellGenomes)

    m3 = Model_BD(final_pop_size, polyp_birth_rate/2., plyp_death_rate/2., mut_rate_polyp/2., adv_mut_rate_polyp/2., s_coef_polyp, cellArr, adv_clones_arr, cellGenomes, parentCellArr, mutId, cloneId, mutInductionTimes)
    
    #~~~~~~~~ Expand population
    cellArr, polyp_cellGenomes, adv_clones_arr, mutId, polyp_driver_muts, model, popsize3 = expandPop(m3, group="polypfission", verbose=verbose, polyp_init_time=Int(polyp_init_time), max_time=Int(max_t) - Int(polyp_init_time))

    origindict = assign_seed_origin(polyp_cellGenomes, seed_genomes)

    df_polyp, vaf_polyp = build_output_table(polyp_cellGenomes, m3.mutInductionTimes, polyp_driver_muts, "Polyp", origindict, mean_depth, sd_depth, sample_size=15000, mut_rate=Int(floor(mut_rate/2.)), purity=purity, samplehist=true, overdispersion=.001)
    print(length(vaf_polyp))
    popsize = vcat(popsize1, popsize2)

    df = vcat(df, df_polyp)

    println("Output Complete ", outfile)

    if verbose
        println("Polyp Growth Complete")
    end
    
    h5open(joinpath(outfile, "out.hdf5"), "w") do file
        g = create_group(file, "params")
        param_names::Vector{String} = []
        values::Vector{String} = []
        for (key, value) in parsed_args
            push!(param_names, String(key))
            push!(values, string(value))
        end
        g["param_names"] = param_names
        g["values"] = values
        attributes(g)["Description"] = "Contains parameters used received in this simulation"

        # Population Size Information
        g = create_group(file, "pops") # create a group
        g["inutero"] = popsize1
        g["fission"] = popsize2
        g["polypfission"] = popsize3
        attributes(g)["Description"] = "Contains population sizes of each growth phase: inutero, fission, polypfission"

        # Mutation delineators for each growth phase
        g = create_group(file, "mutId_delineator") # create a group
        g["delineator_muts"] = mutId_delineator
        attributes(g)["Description"] = "Contains the mutation delineators between each growth phase (i.e. mutID after inutero is mutId_delineator[0])"

        # Seed Genomes
        g = create_group(file, "seed_genomes") # create a group
        for i in 1:length(seed_genomes)
            g["seed_genome_$i"] = seed_genomes[i]
        end
        attributes(g)["Description"] = "Contains the seed genomes used for the polyp growth phase initialization"

        # Driver Mutations
        g = create_group(file, "driver_muts") # create a group
        g["inutero"] = driver_muts
        g["fission"] = driver_muts_fission
        g["polypfission"] = polyp_driver_muts

        # VAFs
        g = create_group(file, "samplehist") # create a group
        g["vafs_inutero"] = vafs1
        g["vafs_fission"] = vafs2
        g["vafs_polyp"] = vaf_polyp
    end

    output_dict = Dict("muts.tsv" => df, "empty" => "empty")

    save_model_outputs(outfile, output_dict)
    
    if verbose
        println("Output Complete ", outfile)
    end
end

main()
