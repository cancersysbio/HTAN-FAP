function incurMuts(parentGenomes, nmuts, mutId, adv_arr, adv_mut_rate)
    ret_genomes = Vector{Vector{Int}}(undef, length(parentGenomes))
    driver_muts = Int[]
    if sum(nmuts) > 0
        new_muts_start = mutId+1 # Mutations gained
    else
        new_muts_start = mutId
    end

    for i in 1:length(parentGenomes)
        pGenome = parentGenomes[i]
        nmut = nmuts[i]
        
        if nmut > 0
            # Julia is inclusive of first and last value in range unlike Python!!!!
            mutId = mutId + 1 # Set the first new mutation ID for this cell to never use same ID!
            cUniqueMuts = mutId:(mutId + nmut)
            
            if adv_arr[i] != 1 && rand() < adv_mut_rate
                adv_arr[i] = 1
                # Randomly sample a mutation to be a driver mutation
                d = sample(cUniqueMuts)
                push!(driver_muts, d)
            end
            
            ret_genomes[i] = vcat(pGenome, cUniqueMuts)

            # Finally set the mutId to the last mutation ID for this cell
            mutId += nmut
        else
            ret_genomes[i] = pGenome
        end
    end

    new_muts_end = mutId
    
    # new_muts_start:new_muts_end is the the mutations at this timestep
    return ret_genomes, adv_arr, mutId, driver_muts, new_muts_start:new_muts_end
end

function incurMutsSingleCell(parentGenome, nmuts, mutId, adv_arr, adv_mut_rate)
    ret_genomes = Vector{Vector{Int}}(undef, length(parentGenomes))
    driver_muts = Int[]

    for cell_muts in 1:length(parentGenomes)
        pGenome = parentGenomes[cell_muts]
        cUniqueMuts = Int[]
        
        if nmuts[cell_muts] > 0
            cUniqueMuts = collect(mutId:mutId+nmuts[cell_muts]-1)
            mutId += nmuts[cell_muts]
            
            if adv_arr != 1
                if rand() < adv_mut_rate
                    adv_arr = 1
                    driver_muts = vcat(driver_muts, mutId-1)
                end
            end
        end
        
        cGenome = vcat(pGenome, cUniqueMuts)
        ret_genomes[cell_muts] = cGenome
    end
    return ret_genomes, adv_arr, mutId, driver_muts
end


# function incurMuts(parentGenomes::Vector{Vector{Int}}, nmuts::Array{Int}, mutId::Int, cellGenomeList::Vector{Vector{Int}})
#     for n in 1:size(nmuts, 1)
#         pGenome = parentGenomes[n]
#         if !isempty(pGenome)
#             cGenome = vcat(pGenome, (mutId:mutId+nmuts[n, 1]-1))
#         else
#             cGenome = (mutId:mutId+nmuts[n, 1]-1)
#         end
#         push!(cellGenomeList, cGenome)
#         mutId += nmuts[n, 1] + 1
#     end
#     return cellGenomeList, mutId
# end