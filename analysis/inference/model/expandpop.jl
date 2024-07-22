using Pkg
using StatsBase
using PoissonRandom
using Distributions
using Random

include("init.jl")
include("incurmuts.jl")

function expandPop(model::Model_BD; group::String="Group", verbose::Bool=false, init_birth::Float64=1.2, polyp_init_time::Int64=10, max_time::Int64=10)
    birth_rate = model.birth_rate
    death_rate = model.death_rate
    mut_rate = model.mut_rate
    adv_mut_rate = model.adv_mut_rate
    mut_rate_store = mut_rate
    adv_mut_rate_store = adv_mut_rate
    s_coef = model.s_coef
    cellArr = model.cellArr
    adv_clones_arr = model.adv_clones_arr
    cellGenomes = model.cellGenomes
    mutId = model.mutId
    cloneId = model.cloneId
    max_pop = model.initSize

    fission_poisson = Poisson(birth_rate)
    adv_fission_poisson = Poisson(birth_rate*(1+s_coef))
    death_poisson = Poisson(death_rate)

    currentPop = length(cellArr)
    prevPop = currentPop
    trigger_stop = false
    t = 0

    # Track population size over time
    pop_size = Int[]
    driver_muts = Int[]

    
    for t in 1:max_time
        # Record prev pop
        push!(pop_size, currentPop)
        # Set the starting mutation ID for this tick

        # Neutral population
        neu_pop = length(cellArr) - sum(adv_clones_arr)
        adv_pop = sum(adv_clones_arr)

        adv_divs = 0 # Holder for number of advantageous divisions in case no adv cells
        if group=="inutero"
            t_offset = (t-9.) / 12.
            neu_divs = exponential_births(neu_pop, init_birth)
            total_divs = neu_divs
            if adv_pop > 0
                adv_divs = exponential_births(adv_pop, init_birth*(1+s_coef))
                total_divs += adv_divs
            end
            ndeaths = sum(rand(death_poisson, prevPop))
            replace = true

            if t<2
                mut_rate=0
                adv_mut_rate=0
            else
                mut_rate=mut_rate_store
                adv_mut_rate=adv_mut_rate_store
            end

        elseif group=="fission" || group=="polypfission"
            t_offset = t-1.
            neu_divs = sum(rand(fission_poisson, abs(length(cellArr)-sum(adv_clones_arr))))
            total_divs = neu_divs
            if sum(adv_clones_arr) > 0
                adv_divs = sum(rand(adv_fission_poisson, sum(adv_clones_arr)))
                total_divs += adv_divs
            end

            if group == "polypfission"
                t_offset = (polyp_init_time-1.) + (t-1.)
                ndeaths = sum(rand(death_poisson, prevPop))
                replace = true
                if prevPop <= 20
                    ndeaths = 0
                end
            else
                # ndeaths = Int(ceil(total_divs * .5))
                ndeaths = sum(rand(death_poisson, prevPop))
                replace = false
            end
            
            if total_divs + prevPop >= max_pop && group != "fission" # Normal turnover not subject to maxpop
                # Stop at max_pop
                total_divs = (total_divs + prevPop) - max_pop
                neu_divs = Int(floor(total_divs * (neu_divs / total_divs)))
                if adv_divs > 0
                    adv_divs = Int(floor(total_divs * (adv_divs / total_divs)))
                end
                trigger_stop = true
                if verbose
                    println("Trigger stop on pop exceeded", max_pop)
                end
            end
        else
            exit("Error: Group not recognized")
        end

        if verbose
            println("Group: ", group, " Birth rate: ", birth_rate, " Fissions: ", total_divs, " Population ", length(cellArr), " Deaths ", ndeaths)
        end

        # Update cells and genomes
        if total_divs > 0 || ndeaths > 0 # If not empty
            cellArr, cellGenomes, adv_clones_arr, mutId, new_driver_muts, new_muts_incl_range, nmuts = update_growth_arr_clean(neu_divs, adv_divs, ndeaths, cellArr, cellGenomes, adv_clones_arr, mut_rate, mutId, adv_mut_rate, replace)
            driver_muts = vcat(driver_muts, new_driver_muts)
            # Record the time that each mutation occurred
            if nmuts > 0
                for m in new_muts_incl_range
                    model.mutInductionTimes[m] = Float64(t_offset)
                end
            end
        end

        # Check array lengths are all the same
        if length(cellArr) != length(cellGenomes) || length(cellArr) != length(adv_clones_arr)
            println("Error: Array lengths are not the same")
            println(length(cellArr), " ", length(cellGenomes), " ", length(adv_clones_arr))
            exit("Exit after birth and death")
        end

        currentPop = length(cellArr)
        
        prevPop = currentPop
        cloneId = currentPop

        if trigger_stop
            break
        end

        t += 1
    end

    # Record final time pop
    push!(pop_size, currentPop)

    model.cloneId = cloneId
    model.mutId = mutId
    return cellArr, cellGenomes, adv_clones_arr, mutId, driver_muts, model, pop_size
end

function sample_poisson(rate::Float64, n::Int)
    return [pois_rand(rate) for i in 1:n]
end

function exponential_births(prevPop, r)
    births = (r+1) * prevPop
    ndivs = Int(ceil(births))
    return(ndivs)
end

function get_cell_actions(ndivs, ndeaths, prevPop, adv_clones_arr, s_coef)
    adv_cells = sum(adv_clones_arr)
    neu_p_birth = ndivs / prevPop - adv_cells
    adv_p_birth = adv_ndivs / adv_cells
    p_death = ndeaths / prevPop

    n = length(adv_clones_arr)
    birth_weights = fill(p_birth, n)
    death_weights = fill(ndeaths / prevPop, n)

    # Get array of 1 for divide and 0 for 
    neu_divs = [rand() < neu_p_birth ? 1 : 0 for i in 1:ndivs]

    return birth_weights, death_weights
end

function update_growth_arr_clean(neu_ndivs, adv_ndivs, ndeaths, cellArr, cellGenomes, adv_clones_arr, mut_rate, mutId, adv_mut_rate, replace=false)
    # Sample the cells that will die
    death_weights = fill(ndeaths / length(cellArr), length(cellArr))
    death_idx = sample(1:length(cellArr), Weights(death_weights), ndeaths, replace=replace)
    death_set = Set(death_idx)
    
    # Remove parent cells that died
    keep_idx = setdiff(1:length(cellArr), death_set)

    new_cellArr = cellArr[keep_idx]
    new_cellGenomes = cellGenomes[keep_idx]
    new_adv_clones_arr = adv_clones_arr[keep_idx]

    # birth_wieghts = birth_wieghts[keep_idx]
    # this_parentCellIdx = sample(1:length(new_cellArr), Weights(birth_wieghts), ndivs, replace=replace)

    # Get idx where new_adv_clones_arr == 0
    neu_cell_idx = findall(x -> x == 0, new_adv_clones_arr)
    neu_div_idx = sample(neu_cell_idx, neu_ndivs, replace=replace)

    # Get idx where new_adv_clones_arr == 1
    adv_cell_idx = findall(x -> x == 1, new_adv_clones_arr)
    adv_div_idx = sample(adv_cell_idx, adv_ndivs, replace=replace)

    # Combine the idx
    this_parentCellIdx = vcat(neu_div_idx, adv_div_idx)

    # Add daughter cells
    daughter_cells = this_parentCellIdx
    daughter_cellGenomes = cellGenomes[this_parentCellIdx]
    daughter_adv_clones_arr = adv_clones_arr[this_parentCellIdx]

    new_cellArr = vcat(new_cellArr, daughter_cells)
    new_cellGenomes = vcat(new_cellGenomes, daughter_cellGenomes)
    new_adv_clones_arr = vcat(new_adv_clones_arr, daughter_adv_clones_arr)

    driver_muts = Int[]
    # # Incur mutations
    
    nmuts = sample_poisson(Float64(mut_rate), length(new_cellGenomes))

    # Check if any cells are going to incur mutations
    if sum(nmuts) > 0
        new_cellGenomes, new_adv_clones_arr, mutId, driver_muts, new_muts_incl_range = incurMuts(new_cellGenomes, nmuts, mutId, new_adv_clones_arr, adv_mut_rate)
    else
        new_muts_incl_range = 0:0
    end

    # Check array lengths are all the same
    if length(new_cellArr) != length(new_cellGenomes) || length(new_cellArr) != length(new_adv_clones_arr)
        println("Error births: Array lengths are not the same")
        println(length(new_cellArr), " ", length(new_cellGenomes), " ", length(new_adv_clones_arr))
        exit()
    end

    return new_cellArr, new_cellGenomes, new_adv_clones_arr, mutId, driver_muts, new_muts_incl_range, sum(nmuts)
end



