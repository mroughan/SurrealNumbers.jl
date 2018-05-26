module SurrealNumbers

import Base: convert, promote, promote_rule, hash, show, read, expand, 
    size, sign, round, floor, ceil, mod, trunc,
    isinteger, isinf, isnan, isfinite, isodd, iseven,
    <, <=, zero, one, ==, <, -, +, *, /

export Surreal, SurrealFinite, Dyadic
export ExistingSurreals
export <=, zero, one, ==, ≅, ≇, <, -, +, *, ϕ, ≀, ↑, ↓, 
    convert, promote, hash,
    pf, pff, spf, show, surreal2dot, surreal2dag, surreal2node, surreal2tex, read, expand,
    generation, canonicalise, iscanonical,
    unique2!, size, size_u, n_zeros, depth_max, depth_av, list_n, count_n, depth, 
    sign, round, floor, ceil, trunc, mod, 
    isinteger, isinf, isnan, isfinite, isodd, iseven, isdivisible

abstract type Surreal <: Real end 


# include("Dyadic.jl")
include("SurrealFinite.jl")
# include("SurrealTrans.jl")
 
end # module 




 
