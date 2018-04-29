module SurrealNumbers

importall Base 

export Surreal, SurrealFinite, Dyadic
export <=, zero, one, ==, ≅, ≇, <, -, +, *, ϕ, ≀, ↑, ↓, 
    convert, promote, hash, pf, pff, show, generation, simplify, sign,
    round, floor, ceil, isinteger, unique2!,
    isinf, isnan, isfinite


abstract type Surreal <: Real end 


# include("Dyadic.jl")
include("SurrealFinite.jl")
# include("SurrealTrans.jl")
 
end # module 




 
