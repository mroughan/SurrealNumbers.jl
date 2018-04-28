importall Base

struct SurrealSequence
    finite :: Boolean
    sup :: Surreal
    inf :: Surreal
    limsup :: Surreal
    liminf :: Surreal
end
<=(x::SurrealSequence, y::SurrealSequence) = x.sup <= y.inf
function +(x::SurrealSequence, y::SurrealSequence)
    # can't add supremums
    # need more information about a sequence
    #   could have it in sorted order
    
end

struct SurrealTrans <: Real
    name::String
    X_L :: SurrealSequence
    X_R :: SurrealSequence
    # constructor should check that X_L <= X_R
    function SurrealTrans(name::String, X_L::SurrealSequence, X_R::SurrealSequence)
        if X_L <= X_R
            return new(name, X_L, X_R) 
        else
            error("Must have X_L <= X_R")
        end 
    end
end 


#
# though to build these usig lazy evaluation, but there are things
# I am not sure how to do then
#
# struct SurrealTrans <: Real
#    name::String
#    X_L :: Function # provide a map from an index to the surreal number 
#    X_R :: Function # provide a map from an index to the surreal number
# end
# 
# ω = SurrealTrans("ω", x->x,   x->ϕ)
# ε = SurrealTrans("ε", x->1/x, x->ϕ)

