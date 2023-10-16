module SurrealNumbers

using DataStructures

import Base: convert, promote, promote_rule, hash, show, read, delete!, empty!, 
    size, sign, round, floor, ceil, mod, trunc,
    sort,
    isinteger, isinf, isnan, isfinite, isodd, iseven, iszero, isless, isequal,
    <, <=, zero, one, ==, <, -, +, *, /,
    summarysize

export Surreal, SurrealFinite, SurrealShort, SurrealDyadic, SurrealIndex
export SurrealZero, SurrealOne, SurrealMinusOne, SurrealTwo, SurrealMinusTwo, SurrealThree, ϕ  

# eventually these will be removed, but they are convenient for testing
export ExistingSurreals, ExistingConversions, ExistingCanonicals, ExistingCanonicalsC, 
        ExistingProducts,
        ExistingSums, ExistingSums2, ExistingSums3,  ExistingSums4,
        ExistingNegations, 
        ExistingEQ, ExistingLEQ, ExistingLEQ1, ExistingLEQ2, ExistingLEQ3, ExistingLEQ4, ExistingLEQ5, # ExistingGT,
        ExistingFloors, ExistingIntegers, ExistingCanonicalIntegers, ExistingSurrealDAGstats,
        Count, CountUncached, RecursionDepth
export clearcache, cache_stats, cache_stats_summary, cache_stats_LEQcount, cache_hit_percent

# export SurrealAlt2  
# export ExistingSurrealsAlt2, ExistingProductsAlt2, ExistingSumsAlt2

export ∅, ϕ, <=, zero, one, ==, ≅, ≇, <, -, +, *, ≀, ↑, ↓, dag_add,
    dali, convert, overbarString, promote, hash, 
    pf, pff, spf, show, surreal2node, surreal2dot, surreal2dag, surreal2tex, read, expand,
    generation, canonicalise, iscanonical, parents, isancestor, ≺, ≻, ⪯, ⪰, 
    unique2!, 
    sign, round, floor, ceil, trunc, mod, 
    isinteger, iscanonicalinteger, isinf, isnan, isfinite, isodd, iseven, isdivisible, isless, isequal, leq, 
    nodes, edges, paths, tree_nodes, breadth, width, surrealDAG,
    uniqueness, uniqueness_max, uniqueness_failure,
    equivtozero

export SurrealStats, SurrealDAGstats, SurrealDegreeDist, dag_stats, degree_stats, print_degree_stats, degree_dist

abstract type Surreal <: Real end 
 
@static if VERSION < v"0.7.0"
    const stdout = STDOUT
end 
  
# include("Dyadic.jl")
include("SurrealFinite.jl")
# include("SurrealTrans.jl") 
# include("SurrealAlt2.jl")

end # module 



 
