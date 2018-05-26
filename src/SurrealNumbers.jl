module SurrealNumbers

importall Base 

export Surreal, SurrealFinite, Dyadic
export ExistingSurreals
export <=, zero, one, ==, ≅, ≇, <, -, +, *, ϕ, ≀, ↑, ↓, 
    convert, promote, hash,
    pf, pff, spf, show, surreal2dot, surreal2dag, surreal2node, surreal2tex, read, expand,
    generation, canonicalise, iscanonical,
    unique2!, size, n_zeros, depth_max, depth_av, count_n, depth, 
    sign, round, floor, ceil, mod, 
    isinteger, isinf, isnan, isfinite, isodd, iseven, isdivisible


abstract type Surreal <: Real end 


# include("Dyadic.jl")
include("SurrealFinite.jl")
# include("SurrealTrans.jl")
 
end # module 




 
