using SurrealNumbers

using Random
using Distributions
using StableRNGs

using SpecialFunctions

using DataFrames

using OffsetArrays

import Statistics.mean
import Base.maximum
import Base.minimum
 
include("utilities.jl")

# P( g(x) ≤ k )
function generation_cdf(k; α=0.5, λ=3.0)
    # cx2 = 0 : 10*maximum(cx)
    # zn = cumsum( OffsetArray(α .^ cx2, cx2 ) )
    # zn = zn ./ zn[end]
    # zn = collect(zn)
    if k==0
        return exp(-λ)
    elseif k>0
        return exp( -λ * α ^ k )    
    else    
        throw(DomainError(k, "k must be >= 0"))
    end 
end 

# P( g(x) = k )
function generation_pmf(k; α=0.5, λ=3.0)
    if k==0
        return exp(-λ)
    elseif k>0
        return generation_cdf(k; α=α, λ=λ) - generation_cdf(k-1; α=α, λ=λ) 
    else    
        throw(DomainError(k, "k must be >= 0"))
    end 
end 

function generation_exp( α, λ;  ϵ = 1.0e-12)  
    if α <= 0.0 || α >= 1.0
        throw(DomainError(α, "Invalid input parameterrs"))
    end
    if λ <= 0.0 || λ > 5
        throw(DomainError(λ, "Invalid input parameterrs"))
    end
    not_end = true
    expectation = 0.0
    k = 0
    while not_end
        tmp = 1.0 - generation_cdf(k; α=α, λ=λ)
        expectation += tmp
        k += 1
        if tmp < ϵ
            not_end = false
        end
    end
    return (expectation, k)
end

function generation_exp_approx1( α, λ )
    if α <= 0.0 || α >= 1.0 || λ <= 0.0 || λ > 5
        throw(DomainError("Invalid input parameterrs"))
    end
    expectation = -( Base.MathConstants.γ + log(λ) - expint(λ)) / log(α)
    return (expectation, 1)
end

function generation_exp_approx3( α, λ; ϵ = 1.0e-2 )
    if α <= 0.0 || α >= 1.0 || λ <= 0.0 || λ > 5
        throw(DomainError("Invalid input parameterrs"))
    end
    (tmp, k) = generation_exp( α, λ;  ϵ = ϵ)
    expectation =  tmp +  λ * α^(k-1)
    return (expectation, k)
end

function generate_canonicals(generation::Integer)
    # create all of the dyadics up to given generation
    list = Vector{SurrealFinite}(undef, 2^(generation+1) - 1)
    list[1] = dali(0)

    m = 2
    for k=1:generation
        cn = mod(k,2) + 1
        list[m] = convert(SurrealFinite, -k)
        m += 1

        list[m] = convert(SurrealFinite, k)
        m += 1

        for f=2:k
            for j=(f-2)*2^(k-f+1) .+ (1:2:2^(k-f+1))

                v = j // 2^(k-f+1)
                list[m] = convert(SurrealFinite, v)
                m += 1

                v = -j // 2^(k-f+1)
                list[m] = convert(SurrealFinite, v)
                m += 1
            end
        end 

    end
    return list
end

function one_step_generation( input::Vector{SurrealFinite}, n::Integer, 
                                rng::AbstractRNG, 
                                Dn::Distributions.DiscreteDistribution,
                                Ds::Vector,
                                α::AbstractFloat)
    input_len = length(input)
    stats = dag_stats.(input; V=false)
    gen = [ s.generation for s in stats ]
    co, cx = count_n(gen)
    g_n =  co ./ input_len
    weights = OffsetArray(α .^ cx, cx ) ./ g_n
    prob = weights[ gen ]    # give a sampling probability to each element of the population
    prob = prob ./ sum(prob) # normalise
    D2 = Categorical( prob )
    output = Vector{SurrealFinite}(undef, n)
    n_parents = Vector{Int64}(undef,n)
    split_point = Vector{Int64}(undef,n)
    for i=1:n 
        n_parents[i] = rand( rng, Dn )
        # n_parents[i] = 3
        if n_parents[i] == 0
            output[i] = SurrealZero
        else
            not_valid = true
            count = 0
            while not_valid && count < 30
                count += 1
                P = Vector{SurrealFinite}(undef, n_parents[i])
                for j=1:n_parents[i]
                    k = rand( rng, D2 ) 
                    P[j] = input[ k ]
                end
                # if we try to get n_parents unique samples, we run into
                # problems if the intial space is small, so instead accept that
                # there will be some non-uniqueness, and drop the duplicates
                P = sort( unique(P) )
                n_parents[i] = length(P)        

                split_point[i] = rand( rng, Ds[n_parents[i]] )
                if split_point[i] > 0
                    L = P[ 1:split_point[i] ]
                else 
                    L = ϕ
                end
                if split_point[i] < n_parents[i]
                    R = P[ split_point[i]+1 : n_parents[i] ]           
                else 
                    R = ϕ
                end
                if length(L)==0 || length(R)==0 || L[end] < R[1] # could get the case that L[end] ≡ R[1]
                    output[i] = SurrealFinite("", L, R)
                    not_valid = false
                end
            end
            if not_valid
                error("taking lots of iterations")
            end
        end 
    end
    return output, n_parents, split_point
end 

mutable struct MeanDAGstats <: SurrealStats
    nodes::Float64 # nodes in its DAG
    tree_nodes::Float64 # nodes in a tree-based representation of the graph
    edges::Float64 # edges in its DAG
    generation::Float64 # generation/birthday of the surreal, i.e., shortest path to root at zero
    # longest_path::Array{SurrealFinite,1} 
    # paths::Float64 # number of paths from source to sink
    value::Float64
    minval::Float64 
    maxval::Float64
    n_zeros::Float64 # number of nodes with value zero: only calculate if V switch is true
    max_parents::Float64
    gap::Float64
end

function compile_stats( X::Vector{SurrealDAGstats}, fn::Function)
    # fn should map a Vector of Reals to a Float, e.g., mean, maximum, minimum, ...
    result = MeanDAGstats(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    for f in fieldnames( MeanDAGstats )
        tmp = fn( [ Float64(getfield(x,f)) for x in X] ) 
        setfield!( result, f, tmp)
    end
    return result
end
mean( X::Vector{SurrealDAGstats} ) = compile_stats( X, mean )
maximum( X::Vector{SurrealDAGstats} ) = compile_stats( X, maximum )
minimum( X::Vector{SurrealDAGstats} ) = compile_stats( X, minimum )

# function mean( X::Vector{SurrealDAGstats} )
#     result = MeanDAGstats(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
#     for f in fieldnames( MeanDAGstats )
#         tmp = mean( [ Float64(getfield(x,f)) for x in X] ) 
#         setfield!( result, f, tmp)
#     end
#     return result
# end

# function maximum( X::Vector{SurrealDAGstats} )
#     result = MeanDAGstats(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
#     for f in fieldnames( MeanDAGstats )
#         tmp = maximum( [ Float64(getfield(x,f)) for x in X] ) 
#         setfield!( result, f, tmp)
#     end
#     return result
# end

# function minimum( X::Vector{SurrealDAGstats} )
#     result = MeanDAGstats(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
#     for f in fieldnames( MeanDAGstats )
#         tmp = minimum( [ Float64(getfield(x,f)) for x in X] ) 
#         setfield!( result, f, tmp)
#     end
#     return result
# end


function m_step_generation( n::Integer, m::Integer, gen::Integer, 
                            α::Real, λ::Real, 
                            seed::Integer;
                            parent_dist = "P", # Poisson distribution is the only case so far
                            split_dist = "U", # default uniform distribution
                            V=true)

    rng = StableRNG(seed)
    if parent_dist == "P"
        Dn = Poisson(λ) # no of parents distriubution
    else
        throw(DomainError(parent_dist))
    end

    maxp = 100
    Ds = Vector{Distributions.DiscreteUnivariateDistribution}(undef, maxp)
    if split_dist == "U"
        for i=1:maxp
            Ds[i] = Distributions.DiscreteUniform(0,i)
        end
    elseif split_dist == "B"
        for i=1:maxp
            Ds[i] = Distributions.Binomial(i,0.5)
        end
    else
        throw(DomainError(split_dist))
    end

    init_list = generate_canonicals(gen)
    input = init_list
    p = []
    s = []
    out_g = []
    output = []
    mean_stats = Vector{MeanDAGstats}(undef, m+1)
    max_stats = Vector{MeanDAGstats}(undef, m+1)
    min_stats = Vector{MeanDAGstats}(undef, m+1)
    stats = dag_stats.( input ) 
    mean_stats[1] = mean( stats )
    max_stats[1] = maximum( stats )
    min_stats[1] = minimum( stats )

    for i=1:m
        if mod(i, 10) == 0
            println("     iteration $i of $m")
        end

        output, p, s = one_step_generation( input, n, rng, Dn, Ds, α)
        input = output

        # collect stats
        stats = dag_stats.(output; V=V)
        mean_stats[i+1] = mean( stats )
        max_stats[i+1] = maximum( stats )
        min_stats[i+1] = minimum( stats )
    end
    
    return output, stats, mean_stats, max_stats, min_stats
end

