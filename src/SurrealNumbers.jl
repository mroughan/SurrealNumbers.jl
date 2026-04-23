module SurrealNumbers

using DataStructures

import Base: convert, promote, promote_rule, hash, show, read, delete!, empty!, 
    size, sign, round, floor, ceil, mod, trunc, abs,
    sort,
    isinteger, isinf, isnan, isfinite, isodd, iseven, iszero, isone, isless, isequal,
    <, <=, zero, one, ==, <, -, +, *, /,
    summarysize,
    AbstractFloat

export Surreal, SurrealFinite, SurrealShort, SurrealDyadic, SurrealIndex
export SurrealZero, SurrealOne, SurrealMinusOne, SurrealTwo, SurrealMinusTwo, SurrealThree, ϕ  

# caching and information on such
export Count, CountUncached, RecursionDepth
export clearcache, cache_stats, cache_stats_summary, cache_stats_summary2, cache_stats_LEQcount, cache_hit_percent
export HashList
export ExistingSurreals, ExistingConversions, ExistingCanonicals, ExistingCanonicalsC, 
        ExistingLEQ, ExistingLEQ1, ExistingLEQ2, ExistingLEQ3, ExistingLEQ4, ExistingLEQ5, ExistingEQ,
        ExistingProducts, 
        ExistingSums, ExistingSums2, ExistingSums3, ExistingSums4, ExistingNegations,
        ExistingFloors, ExistingIntegers, ExistingCanonicalIntegers, ExistingSigns,
        ExistingSurrealDAGstats

# export SurrealAlt2  
# export ExistingSurrealsAlt2, ExistingProductsAlt2, ExistingSumsAlt2

export ∅, ϕ, <=, zero, one, ==, ≅, ≇, <, -, +, *, ≀, ↑, ↓, dag_add,
    dali, convert, overbarString, promote, hash, 
    pf, pff, spf, show, surreal2node, surreal2dot, surreal2dag, surreal2tex, read, expand,
    generation, canonicalise, iscanonical, parents, isancestor, ≺, ≻, ⪯, ⪰, 
    unique2!, 
    sign, round, floor, ceil, trunc, mod, abs,
    isinteger, iscanonicalinteger, isinf, isnan, isfinite, isodd, iseven, isdivisible, isless, isequal, leq, 
    iszero, isone, isminusone,
    nodes, edges, paths, tree_nodes, breadth, width, surrealDAG,
    uniqueness, uniqueness_max, uniqueness_failure,
    equivtozero

# for testing    
export sign_0, sign_1, sign_2, sign_uncached, sign_c

export SurrealStats, SurrealDAGstats, SurrealDegreeDist, dag_stats, degree_stats, print_degree_stats, degree_dist

abstract type Surreal <: Real end 
const HashType = UInt64
if UInt != UInt64
    @warn "Possible hash type mismatch"
end
export HashType
# Julia is using UInt as a hash, but I want hashes to be stable across platforms
# They are heavily used to speed up performance, so I need good tests of hashing, and that is a pain if they
# change between cases.
# One potential problem is that Julia doesn't guarantee that hashes will remain stable over versions, so we
# should probably take their advice and conver to the SHA Standard Library or StableHashTraits.jl
# but that requires quite a lot of work that I don't have a lot of time to complete in the latest update

@static if VERSION < v"0.7.0"
    const stdout = STDOUT
end 
  
# include("Dyadic.jl")
include("SurrealFinite.jl")
# include("SurrealTrans.jl") 
# include("SurrealAlt2.jl")


# function __init__()
#     for C in HashList
#         println("export $C")
#         eval("import SurrealNumbers.$C")
#     end
# end

end # module 



 
