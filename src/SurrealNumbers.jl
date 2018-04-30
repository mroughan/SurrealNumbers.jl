module SurrealNumbers

importall Base 

export Surreal, SurrealFinite, Dyadic
export <=, zero, one, ==, ≅, ≇, <, -, +, *, ϕ, ≀, ↑, ↓, 
    convert, promote, hash, pf, pff, spf, show, generation, canonicalise, unique2!, 
    sign, round, floor, ceil,
    isinteger, isinf, isnan, isfinite


abstract type Surreal <: Real end 


# include("Dyadic.jl")
include("SurrealFinite.jl")
# include("SurrealTrans.jl")
 
end # module 




 
