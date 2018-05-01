# example forms 
#   epsilon = < {1,1/2,1/4,...} | {} >
#   omega = < {0,1,2,3,...} | {} >
# 
#   could do non-dyadic rationals, e.g., 1/3 by truncating binary approximations
#

#
# the idea here is to build the infinite bits of these using lazy evaluation
#       sim. to lazy indexing into a sequence
#           https://github.com/MikeInnes/Lazy.jl
#           http://mikeinnes.github.io/2017/06/05/lazy.html
#
# but we'd still have to be able to compare sup and inf such sequences to make it all work
# so maybe there is a better mathematical way to specify the sequences such that this is easy
#  

struct SurrealTrans <: Surreal
    shorthand::String
    X_L::SurrealSequence
    X_R::SurrealSequence
    # constructor should check that X_L <= X_R
    function SurrealTrans(name::String, X_L::SurrealSequence, X_R::SurrealSequence)
        if X_L < X_R
            return new(name, X_L, X_R) 
        else
            error("Must have X_L <= X_R")
        end 
    end
end 

struct SurrealSequence <: Function
#  provide a map from an index to the surreal number 
#     but its a bit restricted
#        i.e., for inf sequence make it a map from Naturals -> SurrealFinite
#              but for finite or empty sequences do something else?
#  seem similar to indexing
end

# 
# ω = SurrealTrans("ω", x->x,   x->ϕ)
# ε = SurrealTrans("ε", x->1/x, x->ϕ)
# 

# do + using
#  x+y = < X^L + Y^L | X^R + Y^R >
#
# where
#    X + Y = {  x_i + y_j | for all x_i in X, and y_j in Y }
# so
#
#    [X+Y]_I = x_i + y_j
# where we take indexes (i,j) from table
#              j
#      1 3 6 10
#      2 5 9  
#    i 4 8 
#      7 
#      
# and so in (the entry in the table is I)
#
