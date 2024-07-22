using Random

mutable struct Model
    initSize::Int64
    birth_rate::Float64
    mut_rate::Float64
    adv_mut_rate::Float64
    s_coef::Float64
    cellArr::Vector{Int}
    adv_clones_arr::Vector{Int}
    cellGenomes::Vector{Vector{Int}}
    parentCellArr::Array{Int64, 1}
    mutId::Int
    cloneId::Int
    mutInductionTimes::Dict{Int, Float64}

    # Constructor
    function Model(initSize::Int64 = 10^6,
            birth_rate::Float64 = 0.9,
            mut_rate::Float64 = 3.2,
            adv_mut_rate::Float64 = 10^-2,
            s_coef::Float64 = 0.0,
            cellArr::Vector{Int} = [0],
            adv_clones_arr::Vector{Int} = [0],
            cellGenomes::Vector{Vector{Int}} = [[]],
            parentCellArr::Array{Int64} = [],
            mutId::Int = 0,
            cloneId::Int = 1,
            mutInductionTimes::Dict{Int, Float64} = Dict())
        new(initSize, birth_rate, mut_rate, adv_mut_rate, s_coef, cellArr, adv_clones_arr, cellGenomes, parentCellArr, mutId, cloneId, mutInductionTimes)
    end
end

mutable struct Model_BD
    initSize::Int64
    birth_rate::Float64
    death_rate::Float64
    mut_rate::Float64
    adv_mut_rate::Float64
    s_coef::Float64
    cellArr::Vector{Int}
    adv_clones_arr::Vector{Int}
    cellGenomes::Vector{Vector{Int}}
    parentCellArr::Array{Int64, 1}
    mutId::Int
    cloneId::Int
    mutInductionTimes::Dict{Int, Float64}
    

    # Constructor
    function Model_BD(initSize::Int64 = 10^6,
            birth_rate::Float64 = 0.9,
            death_rate::Float64 = 0.1,
            mut_rate::Float64 = 3.2,
            adv_mut_rate::Float64 = 10^-2,
            s_coef::Float64 = 0.0,
            cellArr::Vector{Int} = [0],
            adv_clones_arr::Vector{Int} = [0],
            cellGenomes::Vector{Vector{Int}} = [[]],
            parentCellArr::Array{Int64} = [],
            mutId::Int = 0,
            cloneId::Int = 1,
            mutInductionTimes::Dict{Int, Float64} = Dict())
        new(initSize, birth_rate, death_rate, mut_rate, adv_mut_rate, s_coef, cellArr, adv_clones_arr, cellGenomes, parentCellArr, mutId, cloneId, mutInductionTimes)
    end
end