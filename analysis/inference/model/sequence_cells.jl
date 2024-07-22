using Pkg
# Pkg.add("DataFrames")
# Pkg.add("Distributions")

using DataFrames
using Random
using Distributions

function mutationFrequency(cellGenomeFull, pop_size, seed_origins, raw_threshold=0.01)
    # Flatten mutations
    muts = vcat(cellGenomeFull...)

    # Count mutations
    mutation_counts = countmap(muts)

    # Mutation frequency
    mutation_frequency = Dict{Int, Float64}()
    for (mutation, count) in pairs(mutation_counts)
        mutation_frequency[mutation] = count / (pop_size * 2.0) # diploid
    end

    # Remove mutations with frequency less than threshold
    mutation_frequency = Dict(filter(x -> x[2] >= raw_threshold, mutation_frequency))

    # Rebuild seed_origins dict
    seed_origins_return = Dict{Int, Int}()
    # loop over mutations in the frequency dict
    for (mutation, frequency) in mutation_frequency
        seed_origins_return[mutation] = get(seed_origins, mutation, -1)
    end

    # println("Counts complete.")
    return mutation_frequency, seed_origins_return
end

# function mutationFrequency(cellGenomeFull, pop_size, raw_threshold=0.01)
#     # Flatten mutations
#     muts = vcat(cellGenomeFull...)

#     # Unique mutations
#     unique_muts = unique(muts)

#     # count mutations
#     mutation_counts = Dict{Int, Int}()
#     for mut in unique_muts
#         count = count_instances(mut, cellGenomeFull)
#         # println("Mutation: ", mut, " Count: ", count)
#         mutation_counts[mut] = count
#     end

#     # mutation_counts = Dict{Int, Int}()
#     # for genome in cellGenomeFull
#     #     for mutation in genome
#     #         mutation_counts[mutation] = get(mutation_counts, mutation, 0) + 1
#     #     end
#     # end

#     mutation_frequency = Dict{Int, Float64}()
#     for (mutation, count) in mutation_counts
#         mutation_frequency[mutation] = count / (pop_size * 2.) # diploid
#     end

#     # Remove all mutations with frequency less than threshold
#     mutation_frequency = Dict(filter(x -> x[2] >= raw_threshold, mutation_frequency))

#     println("Counts complete.")
#     return mutation_frequency
# end

function count_instances(value, list_of_lists)
    count = 0
    for inner_list in list_of_lists
        count += sum(x -> x == value, inner_list)
    end
    return count
end

function build_output_table(cellGenomes, mutInductionTimes, driver_muts, group, seed_origin, mean_depth, sd_depth; sample_size=15000, mut_rate=1, purity=1.0, samplehist=true, overdispersion=.001)
    # Get population size
    pop_size = length(cellGenomes)

    # Sample the population
    if sample_size != 0
        sample_size = min(sample_size, pop_size)
        sample_idx = sample(1:pop_size, sample_size, replace=false)
        cellGenomes = cellGenomes[sample_idx]
        pop_size = length(cellGenomes)
    end
    
    raw_vafs, seed_origins = mutationFrequency(cellGenomes, pop_size, seed_origin)
    
    output_table = DataFrame(mutation = Int[], raw_vaf = Float64[], vaf = Float64[], t= Float64[], is_driver = Int[], seed_parent = Int[], group = String[])
    
    for (mutation, frequency) in raw_vafs
        t = get(mutInductionTimes, mutation, -1) # -1 if can't find the mutation induction time in dict
        is_driver = 0
        for d in 1:length(driver_muts)
            if driver_muts[d] == mutation
                is_driver = 1
            else
                is_driver = 0
                # println("Driver mutation found: ", mutation)
            end
        end

        # Get seq_vafs
        vaf = sequence_mutations(frequency, mean_depth, sd_depth)
        seed_num = get(seed_origins, mutation, -1)

        push!(output_table, (mutation, frequency, vaf, t, is_driver, seed_num, group))
    end

    if samplehist
        # Beta binomial and sampled histogram for subclonal mutations
        AF = allelefreq(vcat([cellGenomes[c] for c in 1:length(cellGenomes)]...))
        AF, cmut = allelefreqexpand(AF, mut_rate, [])
        ret = sampledhist(AF, length(cellGenomes), overdispersion, read_depth=mean_depth, cellularity=purity)
        VAF = ret.VAF
    else
        VAF = [0.0]
    end

    return output_table, VAF
end

function sequence_mutations(vaf, mean_depth, sd_depth)
    # Simulate sequencing mutations getting reads, and depth
    
    # Reads from a binomial distribution with vaf as probability and depth as number of trials
    # Depth is normal distribution with mean of 'depth' and standard deviation 10
    # Return the number of reads
    total_reads = Int(floor(rand(Normal(mean_depth, sd_depth))))

    if total_reads < 0 || vaf > 1 || vaf < 0
        println("Total Reads: ", total_reads)
        println("VAF: ", vaf)
        println("ERROR!")
        exit()
    end

    var_reads_dist = Binomial( total_reads, vaf )
    var_reads = rand(var_reads_dist)

    # var_reads = rand(var_reads_dist, rand_total_reads)

    # Divide each var_reads by the total reads
    return var_reads / total_reads
end

function assign_seed_origin(polyp_cellGenomes::Vector{Vector{Int}}, seed_genomes::Vector{Vector{Int}})
    # Create a dictionary to store the origin of each mutation
    origindict = Dict{Int, Int}()

    # Convert seed genomes to sets for faster operations
    seed_sets = [Set(seed) for seed in seed_genomes]

    # Loop through each cell genome
    for (i, cell_genome) in enumerate(polyp_cellGenomes)
        cell_set = Set(cell_genome)
        # Loop through each seed genome set
        for (j, seed_set) in enumerate(seed_sets)
            # Check if the seed genome is a subset of the cell genome
            if issubset(seed_set, cell_set)
                # Assign the origin for each mutation in the cell genome
                for m in cell_genome
                    origindict[m] = j
                end
                break # Once matched, no need to check other seed genomes
            end
        end
    end

    

    return origindict
end

mutable struct SampledData
    DF::DataFrame
    VAF::Array{Float64,1}
    counts::Array{Int64,1}
    depth::Array{Int64,1}
end

function betabinom(p, n, ρ)
    μ = p * n
    shape1 = (μ / n) * ((1 / ρ) - 1)
    shape2 = n * shape1/μ - shape1
    return rand(Binomial(n, rand(Beta(shape1, shape2))))
end

function sampledhist(AF, cellnum, ρ; detectionlimit = 0.1, ploidy = 2.0, read_depth = 100.0, cellularity = 1.0)

    AF = AF./ploidy
    AF = AF .* cellularity
    filter!(x -> x > detectionlimit * cellnum, AF)
    samp_percent = read_depth/cellnum
    #depth = rand(Binomial(cellnum, samp_percent), length(AF))
    depth = rand(Poisson(read_depth), length(AF))
    samp_alleles = map((x, y) -> betabinom(x, y, ρ), AF/cellnum, depth)
    VAF = samp_alleles./depth

    #data for histogram
    x = 0.005:0.01:1.005
    y = fit(Histogram, VAF, x, closed=:right)
    DFhist = DataFrame(VAF = x[1:end-1], freq = y.weights)

    value = SampledData(DFhist, VAF, samp_alleles, depth)
    return(value)
end

function allelefreq(mutations)
    #create dictionary that maps mutation ID to allele frequency
    f = counts(mutations,minimum(mutations):maximum(mutations))
    muts = collect(minimum(mutations):maximum(mutations))
    idx = f .> 0.01
    f = map(Float64, f[idx])
    muts = muts[idx]
    Dict{Int64, Float64}(muts[i]::Int64 => f[i]::Float64 for i in 1:length(f))
end

function allelefreqexpand(AFDict, μ, subclonemutations; fixedmu = true)

    #expand allele frequency given mutation rate and calculate number of mutations in the subclones
    #subclonemutations = convert(Array{Array{Int64,1},1}, subclonemutations)
    if fixedmu == false
        cmuts = zeros(Int64, length(subclonemutations))
        mutfreqs = collect(values(AFDict))
        mutids = collect(keys(AFDict))
        mutations = rand(Poisson(μ), length(mutfreqs))
        AFnew = zeros(Int64, sum(mutations))

        for i in 1:length(cmuts)
        idx = findall((in)(mutids), subclonemutations[i])
        cmuts[i] = sum(mutations[idx])
        end

        j = 0
        for f in 1:length(mutfreqs)
            AFnew[(j + 1): j + mutations[f]] = fill(mutfreqs[f], mutations[f])
            j = j + mutations[f]
        end

    else
        mutfreqs = collect(values(AFDict))
        mutids = collect(keys(AFDict))
        μint = round(Int64, μ)
        mutations = fill(Int64(μ), length(mutfreqs))
        AFnew = zeros(Int64, sum(mutations))
        cmuts = zeros(Int64, length(subclonemutations))

        for i in 1:length(cmuts)
        idx = findall((in)(mutids), subclonemutations[i])
        cmuts[i] = sum(mutations[idx])
        end

        j = 0
        for f in 1:length(mutfreqs)
        AFnew[(j + 1): j + mutations[f]] = fill(mutfreqs[f], mutations[f])
        j = j + mutations[f]
        end
    end

    return AFnew, cmuts
end