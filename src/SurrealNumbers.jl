module SurrealNumbers

import Base: convert, promote, promote_rule, hash, show, read, expand, 
    size, sign, round, floor, ceil, mod, trunc, 
    isinteger, isinf, isnan, isfinite, isodd, iseven, iszero, 
    <, <=, zero, one, ==, <, -, +, *, /

export Surreal, SurrealFinite, SurrealShort, SurrealDyadic, SurrealDAGstats
export ExistingSurreals, ExistingProducts, ExistingSums, Count
  
# export SurrealAlt2  
# export ExistingSurrealsAlt2, ExistingProductsAlt2, ExistingSumsAlt2

export <=, zero, one, ==, ≅, ≇, <, -, +, *, ϕ, ≀, ↑, ↓, dag_add,
    dali, convert, promote, hash,
    pf, pff, spf, show, surreal2node, surreal2dot, surreal2dag, surreal2tex, read, expand,
    generation, canonicalise, iscanonical, parents, isancestor, ≺, ≻, ⪯, ⪰, 
    unique2!, size, size_u, n_zeros, depth_max, depth_av, list_n, count_n, depth, 
    sign, round, floor, ceil, trunc, mod, 
    isinteger, isinf, isnan, isfinite, isodd, iseven, isdivisible,
    dag_stats, surrealDAG, breadth, width, clearcache

abstract type Surreal <: Real end 
 
 
# include("Dyadic.jl")
include("SurrealFinite.jl")
# include("SurrealTrans.jl")
# include("SurrealAlt2.jl")

end # module 




 
