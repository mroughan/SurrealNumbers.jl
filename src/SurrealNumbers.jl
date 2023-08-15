module SurrealNumbers

import Base: convert, promote, promote_rule, hash, show, read, delete!, empty!, 
    size, sign, round, floor, ceil, mod, trunc, 
    isinteger, isinf, isnan, isfinite, isodd, iseven, iszero, 
    <, <=, zero, one, ==, <, -, +, *, /,
    summarysize

export Surreal, SurrealFinite, SurrealShort, SurrealDyadic, SurrealDAGstats
export SurrealZero, SurrealOne, SurrealMinusOne, SurrealTwo, SurrealMinusTwo, SurrealThree, ϕ  

# eventually these will be removed, but they are convenient for testing
export ExistingSurreals, ExistingConversions, ExistingCanonicals, ExistingProducts, ExistingSums, 
        ExistingNegations, ExistingLEQ, ExistingGEQ, ExistingEQ, 
        Count, CountUncached
export clearcache, cache_stats, cache_stats_summary

# export SurrealAlt2  
# export ExistingSurrealsAlt2, ExistingProductsAlt2, ExistingSumsAlt2

export ∅, ϕ, <=, zero, one, ==, ≅, ≇, <, -, +, *, ≀, ↑, ↓, dag_add,
    dali, convert, overbarString, promote, hash, 
    pf, pff, spf, show, surreal2node, surreal2dot, surreal2dag, surreal2tex, read, expand,
    generation, canonicalise, iscanonical, parents, isancestor, ≺, ≻, ⪯, ⪰, 
    unique2!, 
    sign, round, floor, ceil, trunc, mod, 
    isinteger, isinf, isnan, isfinite, isodd, iseven, isdivisible,
    dag_stats, nodes, edges, paths, tree_nodes, breadth, width, surrealDAG,
    uniqueness, uniqueness_max, uniqueness_failure,
    equivtozero
  
abstract type Surreal <: Real end 
 
@static if VERSION < v"0.7.0"
    const stdout = STDOUT
end 
  
# include("Dyadic.jl")
include("SurrealFinite.jl")
# include("SurrealTrans.jl") 
# include("SurrealAlt2.jl")

end # module 




 
