mutable struct SurrealFinite <: Surreal
    # mutable means I can set the has value when I need it, instead of
    #   (i) always having to calculate it recursively
    #   (ii) calculating it for even temporary surreals
    shorthand::String 
    L::Array{SurrealFinite,1}  
    R::Array{SurrealFinite,1}
    h::UInt64 # this is only set the first time the hash function is called
    # constructor should check that L < R
    function SurrealFinite(shorthand::String, L::Array{SurrealFinite}, R::Array{SurrealFinite}, h::UInt64)
        global Count
        Count['c'] += 1 
        if length(L) > 1
            L = sort( unique(L), by=x->(x,hash(x)) ) # use hash as tie break in sort so that order is deterministic
        end 
        if length(R) > 1
            R = sort( unique(R), by=x->(x,hash(x)) ) # use hash as tie break in sort so that order is deterministic
        end 
        # println("L = $L, R = $R") 
        # use the fact they are sorted to not do a complete comparison
        # also means that == is much easier than if they were unsorted
        if isempty(L) || isempty(R) || L[end] < R[1]
            return new(shorthand, L, R, h)
        else 
            error("Surreal number must have L < R (currently they are $(L) and $(R) )") 
        end  
    end 
end 
SurrealFinite(shorthand::String, L::Array, R::Array ) =
    SurrealFinite( shorthand, convert(Array{SurrealFinite},L), convert(Array{SurrealFinite},R), zero(UInt64))
SurrealFinite(L::Array, R::Array ) =
    SurrealFinite( "", convert(Array{SurrealFinite},L), convert(Array{SurrealFinite},R), zero(UInt64))
≀(L::Array, R::Array) = SurrealFinite(L::Array, R::Array )
 
const SurrealShort  = SurrealFinite 
const SurrealDyadic = SurrealFinite 

function hash(x::SurrealFinite, h::UInt)
    if x.h == 0
        x.h = hash(x.L, convert(UInt64, 0) ) * hash(x.R, convert(UInt64, 1))
    end
    return hash(x.h, h)   
end
function hash(X::Array{SurrealFinite}, h::UInt)
    if isempty(X)
        hash(0, h) 
    elseif length(X) == 1
        hash(X[1], h) 
    else
        hash(X[1], h) * hash( X[2:end], h) # note, order is technically not important
    end
end
# really slow version: hash(X::Array{SurrealFinite}, h::UInt) = isempty(X) ? hash(0,h) : hash( prod(hash.(X, h )), h)
hash(x::SurrealFinite) = hash(x, zero(UInt64) )
hash(x::SurrealFinite, h::Integer) = hash(x, convert(UInt64, h) )
hash(X::Array{SurrealFinite}) = hash(X, zero(UInt64) )
 
# global dictionaries to avoid repeated calculations of the same things
#   in particular to speed up calculations, but also to ensure that resulting
#   structures are actually DAGs, i.e., pointers to the same "number" point at the same bit of memory
# i guess eventually these should be replaced by let-closures: https://docs.julialang.org/en/v0.6.2/manual/variables-and-scoping/
const ExistingSurreals = Dict{UInt64, SurrealFinite}()
const ExistingConversions = Dict{UInt64,Rational}()
const ExistingCanonicals = Dict{Rational,UInt64}()
# ExistingProducts = Dict{SurrealFinite, Dict{SurrealFinite,SurrealFinite}}()
const ExistingProducts   = Dict{UInt64, Dict{UInt64,UInt64}}()
const ExistingSums       = Dict{UInt64, Dict{UInt64,UInt64}}() 
const ExistingNegations  = Dict{UInt64, UInt64}() 
const Count = Dict{Char, Integer}('+'=>0, '*'=>0, '-'=>0, 'c'=>0, '='=>0, '≦'=>0)

# function size(d::Dict)
#     [length(d[k]) for k in sort(collect(keys(d)))]        
# end 
function reset!(d::Dict)
    for k in keys(d)
        d[k] = 0
    end
end 

function clearcache() 
    global ExistingSurreals 
    global ExistingConversions 
    global ExistingCanonicals
    global ExistingProducts
    global ExistingSums
    global ExistingNegations 
    global Count
    empty!(ExistingSurreals)
    empty!(ExistingConversions)
    empty!(ExistingCanonicals)
    empty!(ExistingProducts)
    empty!(ExistingSums)
    empty!(ExistingNegations)
    reset!(Count)
    return 1
end

function convert(::Type{SurrealFinite}, n::Int )  
    global ExistingSurreals 
    global ExistingCanonicals 
    global SurrealZero
    global SurrealOne  
    if n==0 
        result = SurrealZero 
    elseif n==1 
        result = SurrealOne
    elseif haskey(ExistingCanonicals, Rational(n)) 
        return ExistingSurreals[ ExistingCanonicals[Rational(n)] ]
    elseif n>1
        result = SurrealFinite(string(n), [convert(SurrealFinite, n-1)], ϕ )
    else
        result = SurrealFinite(string(n), ϕ, [convert(SurrealFinite, n+1)] )
    end
    hr = hash(result)
    if haskey(ExistingSurreals, hr)
       result = ExistingSurreals[hr] # don't double up on memory
    else
       ExistingSurreals[hr] = result
    end
    ExistingCanonicals[Rational(n)] = hr 
    return result 
    # should check to make sure abs(n) is not too big, i.e., causes too much recursion
end 
function convert(::Type{SurrealFinite}, r::Rational )  
    global ExistingSurreals 
    global ExistingCanonicals 
    if haskey(ExistingCanonicals, r)  
        return ExistingSurreals[ ExistingCanonicals[r] ]
    elseif isinteger(r)  
        result = convert(SurrealFinite, convert(Int64,r) )
    elseif ispow2(r.den)
        # Non-integer dyadic numbers
        n = convert(Int64, round( log2(r.den) ))
        if abs(r) + n < 60 # 60 is a bit arbitrary -- could experiment more on real limits
            result = SurrealFinite("$r", [convert(SurrealFinite, r - 1//2^n)], [convert(SurrealFinite, r + 1//2^n)] )
        else 
            error("Generation too large")
        end 
    else 
        error("we can't do these yet as they require infinite sets")
        # return convert(SurrealFinite, r.num) // convert(SurrealFinite, r.den)
    end  
    hr = hash(result)
    if haskey(ExistingSurreals, hr)
       result = ExistingSurreals[hr] # don't double up on memory
    else
       ExistingSurreals[hr] = result
    end
    ExistingCanonicals[r] = hr 
    return result
end
function convert(::Type{SurrealFinite}, f::Float64 ) 
    if isinteger(f) 
        return convert(SurrealFinite, convert(Int64,f) )
    elseif isfinite(f)
        # WARNING this could be very slow and memory hungry as exact representations could
        #    be up to 2^52 generations out, so can even cause seg fault, e.g., for 1122342342.23422522
        s = sign(f)
        c = significand(f) 
        q = exponent(f)
        @static if VERSION < v"0.7.0"
            t = bits(f)[13:end]
        else
            t = bitstring(f)[13:end]
        end
        p1 = parse(Int, "0b1" * t)
        p2 = p1 * 2.0^(-52 + q)
        r = convert(Rational,s) * p1 // 2^(52 - q)  # rational representation of the Float
        if abs(r.num/r.den) + log2(abs(r.den)) < 60 
            # fractional case 
            return convert(SurrealFinite, r)
        else 
            error("generation too large")
        end
    else 
        @static if VERSION < v"0.7.0"
            throw(DomainError())
        else
            throw(DomainError(f))
        end
    end 
end 
dali(x) = convert(SurrealFinite, x)
SurrealFinite(x) = convert(SurrealFinite, x)

function convert(::Type{Rational}, s::SurrealFinite ) 
    global ExistingConversions
    h = hash(s)
    if haskey(ExistingConversions, h)
        return ExistingConversions[h]
    else 
        if s ≅ zero(s)
            result = 0 // 1
        elseif s ≅ one(s)
            result = 1 // 1
        elseif s < zero(s)
            result = -convert(Rational, -s)
        elseif (sf = floor(Integer, s)) ≅ s
            result = Rational( sf )
        else # 0 < x < 1
            # do a binary search from top down, first valid is the simplest
            xl = isempty(s.L) ? convert(SurrealFinite, sf)   : maximum(s.L)
            xr = isempty(s.R) ? convert(SurrealFinite, sf+1) : minimum(s.R) 
            not_end = true
            k = 0
            a = sf 
            b = sf+1
            while not_end && k < 24 # 24 is arbitrary, but it would be painful to go lower
                k += 1
                d = (b+a) // 2
                # print("k=$k, a=$a, b=$b, d=$d \n")
                
                c = convert(SurrealFinite,  d)
                if xl < c < xr
                    result = Rational( d )
                    not_end = false 
                elseif c <= xl
                    a = d 
                elseif c >= xr
                    b = d
                else
                    error("this case should not happen")
                end 
            end
        end 
    end 
    ExistingConversions[h] = result
    return result
end


# some catch alls
convert(::Type{T}, s::SurrealFinite ) where {T <: AbstractFloat} = convert(T, convert(Rational, s) )
convert(::Type{T}, s::SurrealFinite ) where {T <: Integer} = convert(T, convert(Rational, s) )
convert(::Type{Rational{T}}, s::SurrealFinite ) where {T <: Integer} = convert(Rational{T}, convert(Rational, s) )
AbstractFloat(x::SurrealFinite) = convert(AbstractFloat, x)
# Int64(x::SurrealFinite) = convert(Int64, x)

function convert(::Type{String}, s::SurrealFinite ) 
    # try to work out a nice way to print it
    if isinteger(s)
        return string(convert(Int, s))
    else
        r = convert(Rational, s)
        return string(r.num) * "/" * string(r.den)
    end
end

# promote all numbers to surreals for calculations
promote_rule(::Type{T}, ::Type{SurrealFinite}) where {T<:Real} = SurrealFinite
 

@static if VERSION < v"0.7.0"
    const ϕ = Array{SurrealFinite,1}(0) # empty array of SurrealFinites
    const ∅ = ϕ 
else
    const ϕ = Array{SurrealFinite,1}(undef,0) # empty array of SurrealFinites
    const ∅ = ϕ
end 
const SurrealZero = SurrealFinite("0", ϕ, ϕ ) 
const SurrealOne  = SurrealFinite("1", [ zero(SurrealFinite) ], ϕ ) 
const SurrealMinusOne  = SurrealFinite("-1", ϕ, [ zero(SurrealFinite) ] ) 
const SurrealTwo  = SurrealFinite("2", [SurrealOne], ϕ ) 
const SurrealThree  = SurrealFinite("3", [SurrealTwo], ϕ )  
const SurrealMinusTwo  = SurrealFinite("-2", ϕ, [SurrealOne] ) 
zero(::SurrealFinite) = SurrealZero # always use the same zero
one(::SurrealFinite)  = SurrealOne  # always use the same one
# ↑ = one(SurrealFinite)  # this causes an error???
# ↓ = -one(SurrealFinite)  

# relations
#   these are written in terms of the definition, but could
#   rewrite in terms of max/min to make marginally faster
#   or in terms of set operations to make more succinct
function leq(x::SurrealFinite, y::SurrealFinite, processed_list::Dict{UInt64,Bool})
#    for t in x.L
#        if y <= t
#            return false
#        end
#    end
#    for t in y.R
#        if t <= x
#            return false 
#        end 
#    end
    # global Count
    # Count['≦'] += 1

    k = hash(x,0) * hash(y,1)
    if !haskey(processed_list,k)
        # if !isempty(x.L) && y <= x.L[end] 
        if !isempty(x.L) && leq(y, x.L[end],  processed_list)
            processed_list[k] = false
        # elseif !isempty(y.R) && x >= y.R[1] 
        elseif !isempty(y.R) && leq(y.R[1], x,  processed_list)
            processed_list[k] = false 
        else
            processed_list[k] = true
        end 
    end 
    return processed_list[k]
end 
<=(x::SurrealFinite, y::SurrealFinite) = leq(x, y, Dict{UInt64,Bool}())
<(x::SurrealFinite, y::SurrealFinite) = x<=y && !(y<=x)
# ===(x::SurrealFinite, y::SurrealFinite) = x<=y && y<x # causes an error
#   === is 'egal', and hardcoded for mutables to test they are same object in memory
≅(x::SurrealFinite, y::SurrealFinite) = x<=y && y<=x
≅(x::Real, y::Real) = ≅(promote(x,y)...)
≇(x::SurrealFinite, y::SurrealFinite) = !( x ≅ y ) 
≇(x::Real, y::Real) = ≇(promote(x,y)...)
# ==(x::SurrealFinite, y::SurrealFinite) = size(x.L) == size(y.L) &&
#                                         size(x.R) == size(y.R) &&
#                                         all(x.L .== y.L) &&
#                                         all(x.R .== y.R)

# much faster version because doesn't evaluate the whole array every time
#   potential bug because sort order of equivalent values is not defined
function equals(x::SurrealFinite, y::SurrealFinite, processed_list::Dict{UInt64,Bool})
    # global Count
    # Count['='] += 1
    k = hash(x,0) * hash(y,1)
    if !haskey(processed_list,k)
        if length(x.L) != length(y.L) || length(x.R) != length(y.R)
            processed_list[k] = false
            return processed_list[k]
        end
        for i=1:length(x.L)
#             if x.L[i] != y.L[i]
            if !equals(x.L[i], y.L[i], processed_list)
                processed_list[k] = false
                return processed_list[k]
            end  
        end
        for i=1:length(x.R)
#            if x.R[i] != y.R[i]
            if !equals(x.R[i], y.R[i], processed_list)
                processed_list[k] = false
                return processed_list[k]
            end 
        end 
        processed_list[k] = true
    end
    return processed_list[k]
end
==(x::SurrealFinite, y::SurrealFinite) = x===y || equals(x, y, Dict{UInt64,Bool}())
iszero(x::SurrealFinite) = isempty(x.L) && isempty(x.R)
 
# comparisons between sets (i.e., arrays) are all-to-all, so
#   (1) don't have to be the same size
#   (2) the < is not exactly the same as the < defined above
function <=(X::Array{SurrealFinite}, Y::Array{SurrealFinite} )
    if isempty(X) || isempty(Y)
        return true
    else
        for x in X
            for y in Y
                if !(x <= y) 
                    return false
                end 
            end 
        end
        return true
    end
end
function <(X::Array{SurrealFinite}, Y::Array{SurrealFinite} )
    if isempty(X) || isempty(Y)
        return true
    else
        for x in X
            for y in Y 
                if !( x < y )
                    return false
                end 
            end
        end
        return true
    end
end 

# unary operators
function -(x::SurrealFinite)
    global ExistingSurreals 
    global ExistingNegations
    global Count
    hx = hash(x) # build the dictionary in terms of hashs, because it is used quite a bit
    if iszero(x) 
        result = x
    elseif haskey(ExistingNegations, hx)
        return ExistingSurreals[ ExistingNegations[hx] ]
    elseif isempty(x.shorthand)
        result = SurrealFinite("", -x.R, -x.L )
    elseif x.shorthand[1] == '-'
        result = SurrealFinite(x.shorthand[2:end], -x.R, -x.L )
    else
        result = SurrealFinite("-"*x.shorthand, -x.R, -x.L )
    end  
    hr = hash(result)
    if haskey(ExistingSurreals, hr)
        result = ExistingSurreals[hr] # don't double up on memory
    else
        ExistingSurreals[hr] = result
    end
    if !haskey(ExistingSurreals, hx)
        ExistingSurreals[hx] = x
    end 
    ExistingNegations[hx] = hr
    ExistingNegations[hr] = hx
    Count['-'] += 1
    return result
end

# somewhat incomplete (for reasons described), slightly a cheat, and potentially very slow
function /(x::SurrealFinite, y::SurrealFinite)
    xr = convert(Rational, x)
    yr = convert(Rational, y) 
    if y ≅ zero(y)
        error(InexactError)
    elseif y ≅ 2
        return x * convert(SurrealFinite,  1 // 2)
    elseif isinteger(y) && ispow2(yr.num)
        return x * convert(SurrealFinite,  1 // yr) 
    elseif isinteger(xr.num/yr)
        return convert(SurrealFinite,  (xr.num/yr) // xr.den  ) 
    else
        error(InexactError)        
    end
end

# binary operators
function +(x::SurrealFinite, y::SurrealFinite)
    global ExistingSurreals 
    global ExistingSums
    global Count
    hx = hash(x) # build the dictionary in terms of hashs, because it is used quite a bit
    hy = hash(y) #   and these amortise the cost of initial calculation of hashs 
    if haskey(ExistingSums, hx) && haskey(ExistingSums[hx], hy)
        return ExistingSurreals[ ExistingSums[hx][hy] ]
    elseif iszero(x)
        result = y   
    elseif iszero(y)
        result = x 
    else 
#       result = SurrealFinite( [x.L .+ y; x .+ y.L],
        #                               [x.R .+ y; x .+ y.R] )
        shorthand = "($(x.shorthand) + $(y.shorthand))"
        result = SurrealFinite(shorthand,
                               [[x + y for x in x.L]; [x + y for y in y.L]],
                               [[x + y for x in x.R]; [x + y for y in y.R]] )
    end    
    hr = hash(result)
    # println(" $x+$y:   $hr   $result")   
    if haskey(ExistingSurreals, hr)
       result = ExistingSurreals[hr] # don't double up on memory
    else
       ExistingSurreals[hr] = result
    end 
    # ExistingSurreals[hr] = result  # swap to just using this if we want to see all the ways an addition can given the same result

    if !haskey(ExistingSums, hx)
        ExistingSums[hx] = Dict{UInt64,UInt64}()
    end
    if !haskey(ExistingSums, hy)
        ExistingSums[hy] = Dict{UInt64,UInt64}() 
    end
    ExistingSums[hx][hy] = hr
    ExistingSums[hy][hx] = hr
    Count['+'] += 1
    return result
end
 
# +(x::SurrealFinite, y::SurrealFinite) = SurrealFinite([x.L .+ y; x .+ y.L],
#                                                       [x.R .+ y; x .+ y.R] )

# can't do like this because of empty arrays I think
#+(x::SurrealFinite, y::SurrealFinite) = SurrealFinite([x.L + y.L],
#                                                      [x.R + y.R] )
+(X::Array{SurrealFinite}, Y::Array{SurrealFinite}) = vec([s+t for s in X, t in Y])
 
-(x::SurrealFinite, y::SurrealFinite) = x + -y
-(X::Array{SurrealFinite}, Y::Array{SurrealFinite}) = X + -Y

function *(x::SurrealFinite, y::SurrealFinite)
    global ExistingSurreals 
    global ExistingProducts
    global Count
    hx = hash(x)
    hy = hash(y)
    if haskey(ExistingProducts, hx) && haskey(ExistingProducts[hx], hy)
        return ExistingSurreals[ ExistingProducts[hx][hy] ]
    elseif iszero(x) || iszero(y) 
        result = zero(x)
#     elseif x == 1 
#        result = y 
#    elseif y == 1 
#        result = x
    else  
        # tmp1 = vec([s*y + x*t - s*t for s in x.L, t in y.L])
        # tmp2 = vec([s*y + x*t - s*t for s in x.R, t in y.R])
        # tmp3 = vec([s*y + x*t - s*t for s in x.L, t in y.R])
        # tmp4 = vec([s*y + x*t - s*t for s in x.R, t in y.L])
        # L = [ tmp1; tmp2 ] 
        # R = [ tmp3; tmp4 ] 
        shorthand = "($(x.shorthand) * $(y.shorthand))"
        L = [vec([s*y + x*t - s*t for s in x.L, t in y.L]);
             vec([s*y + x*t - s*t for s in x.R, t in y.R])]
        R = [vec([s*y + x*t - s*t for s in x.L, t in y.R]);
             vec([s*y + x*t - s*t for s in x.R, t in y.L])]
        # L = [subtimes( x, y, x.L, y.L ); subtimes( x, y, x.R, y.R )]
        # R = [subtimes( x, y, x.L, y.R ); subtimes( x, y, x.R, y.L )]
        result = SurrealFinite(shorthand, L, R)
    end
    if !haskey(ExistingProducts, hx)
        ExistingProducts[hx] = Dict{UInt64,UInt64}()
    end
    if !haskey(ExistingProducts, hy)
        ExistingProducts[hy] = Dict{UInt64,UInt64}()
    end 
    hr = hash(result)
    if haskey(ExistingSurreals, hr)
        result = ExistingSurreals[hr] # don't double up on memory
    else 
        ExistingSurreals[hr] = result
    end
    ExistingProducts[hx][hy] = hr
    ExistingProducts[hy][hx] = hr
    Count['*'] += 1 
    return result
end
#function subtimes( x, y, X::Array{SurrealFinite}, Y::Array{SurrealFinite} )
#    n = length(X)
#    m = length(Y)
#    result = Array{SurrealFinite}(m*n)
#    for i=1:n
#        for j=1:m
#            result[(i-1)*m + j] = X[i]*y + x*Y[j] - X[i]*Y[j]
#        end
#    end 
#    return result
#end 
*(x::SurrealFinite, Y::Array{SurrealFinite}) = [ x*s for s in Y ]
*(X::Array{SurrealFinite}, y::SurrealFinite) = y*X


#####################################################3
# read/write routines for surreals

# print the first level in full (ignoring the top level shorthand if present)
#   these should be replaced using "expand"
pf(io::IO, x::SurrealFinite) = println(io, "{ ", x.L, " | ", x.R, " }")
pf(x::SurrealFinite) = pf(stdout, x) 

"""
    expand(x::SurrealFinite; level=0) 

 Writes a surreal as a string with varying levels of expansion.

## Arguments
* `x::SurrealFinite`: the number of elements to expand
* `level=0`: the amount of expansion
    
    + 0 : write shorthand if it exists, or ``\\{ X_L \\| X_R \\}`` if not
    + 1 : ``\\{ X_L \\| X_R \\}``
    + 2 : expand out ``X_L`` and ``X_R`` recursively

## Examples
```jldoctest
julia> expand( convert(SurrealFinite, 2))
"2"
julia> expand( convert(SurrealFinite, 2); level=1)
"{ 1 | ∅ }"
julia> expand( convert(SurrealFinite, 2); level=2)
"{ { { ∅ | ∅ } | ∅ } | ∅ }"
```
"""
function expand(x::SurrealFinite; level=0)
    if level==0
        s = x.shorthand != "" ? x.shorthand : expand(x; level=1)
    elseif level==1
        tmpL = isempty(x.L) ? "∅" : join(convert.(String, x.L), ',')
        tmpR = isempty(x.R) ? "∅" : join(convert.(String, x.R), ',')
        s = "{ " * tmpL * " | " * tmpR  * " }"
        return s 
    elseif level>=1
        return "{ " * expand(x.L;level=level) * " | " * expand(x.R;level=level) * " }"
    end
end
expand(X::Array{SurrealFinite}; level=0) = isempty(X) ? "∅" : join(expand.(X; level=level), ',') 
 
# isnumeric(s::AbstractString) = ismatch(r"^-?\d*.\d*$", s) || ismatch(r"^-?\d*//-?\d*$", s)

# this has to parse a surreal written into a string 
function convert(::Type{SurrealFinite}, s::AbstractString ) 
    # interpret, (i) numbers as canonical, (ii) phi, \phi, ϕ correctly, (iii) structure
    @static if VERSION < v"0.7.0"
        s = replace(s, r"\s+", "") # remove white space
        s = replace(s, r"phi|\\phi|ϕ|∅|\{\}", "") # replace phi or \phi with empty set
    else
        s = replace(s, r"\s+" => "") # remove white space
        s = replace(s, r"phi|\\phi|ϕ|∅|\{\}" => "") # replace phi or \phi with empty set
    end
    # println("     s1 = $s")
    not_end = true
    basic_number = r"\{([^{}|]*)\|([^{|}]*)\}"
    while not_end
        @static if VERSION < v"0.7.0"
            if ismatch(basic_number, s)
                s = replace(s, basic_number, s"SurrealFinite([\1],[\2])") # remove comments
            else 
                not_end = false
            end
        else
            if occursin(basic_number, s)
                s = replace(s, basic_number => s"SurrealFinite([\1],[\2])") # remove comments
            else 
                not_end = false
            end
        end
    end
    # println("     s2 = $s")
    return eval(Meta.parse(s))
end

# read in the full format
function read(io::IO, ::Type{SurrealFinite}, n::Int=1)
    @static if VERSION < v"0.7.0"
        X = Array{SurrealFinite,1}(n)
    else 
        X = Array{SurrealFinite,1}(undef,n) 
    end  
    k = 1
    while !eof(io) && k<=n
        # println(" $k:   $line")
        @static if VERSION < v"0.7.0"
            line = replace(readline(io), r"#.*", s"") # remove comments
            if ismatch(r"\S", line)
                X[k] = convert(SurrealFinite, line)
                k += 1
            end 
        else
            line = replace(readline(io), r"#.*" => s"") # remove comments
            if occursin(r"\S", line)
                X[k] = convert(SurrealFinite, line)
                k += 1
            end
        end
    end 
    return X
end
# read(filename::AbstractString, args...) = open(io->read(io, args...), filename)

# write out a string suitable for inclusion into latex docs
function surreal2tex(io::IO, x::SurrealFinite; level=0)
    s = expand(x; level=level)
    @static if VERSION < v"0.7.0"
        s = replace(s, r"\{", " \\{ ") 
        s = replace(s, r"ϕ", " \\emptyset ") 
        s = replace(s, r"∅", " \\emptyset ") 
        s = replace(s, r"\|", " \\mid ") 
        s = replace(s, r"\}", " \\} ")
        s = replace(s, r"(\d+)//(\d+)", s"\\frac{\1}{\2}")
    else
        s = replace(s, r"\{" => " \\{ ") 
        s = replace(s, r"ϕ" => " \\emptyset ") 
        s = replace(s, r"∅" => " \\emptyset ") 
        s = replace(s, r"\|" => " \\mid ") 
        s = replace(s, r"\}" => " \\} ")
        s = replace(s, r"(\d+)//(\d+)" => s"\\frac{\1}{\2}")
    end
    println(io,s)
end
surreal2tex(x::SurrealFinite; level=0) = surreal2tex(stdout, x; level=level)
 
# standard show will use shorthand when available
function show(io::IO, x::SurrealFinite)
    if io==stdout && x.shorthand != ""       
        @static if VERSION < v"0.7.0"
            print_with_color(:bold, io, x.shorthand ) # could be :red
        else
            printstyled(io, x.shorthand; bold=true)
        end
    elseif x.shorthand != ""
        print(io, x.shorthand )
    else
        # print( io, "<", x.L, ":", x.R, ">")
        print( io, "{ ", x.L, " | ", x.R, " }")
    end
end
# show(io::IO, X::Array{SurrealFinite}) = print(io, "{", join(X, ", "), "}")
function show(io::IO, X::Array{SurrealFinite})
    if isempty(X)
        print(io, "∅")
    else
        print(io, join(X, ", "))
    end
end
# special "canonicalised" output
spf(x::SurrealFinite) = print("{ ", canonicalise.(x.L), " | ", canonicalise.(x.R), " }")

""" 
    surreal2dag(x::SurrealFinite)
    surreal2dag(io::IO, x::SurrealFinite)

 Writes a surreal representation as a DAG out in DOT format for drawing using GraphVis,
 and returns the number of nodes in the graph. Returns the number of nodes in the graph. 

## Arguments
* `io::IO`: output stream, default is stdout
* `x::SurrealFinite`: the number to write out
    
## Examples
```jldoctest
julia> surreal2dag(convert(SurrealFinite, 0))
digraph "0.0" {
   node_1 [shape=none,margin=0,label=
         <<TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0" CELLPADDING="4">
         <TR><TD COLSPAN="2">0</TD></TR>
         <TR><TD PORT="L"> ∅ </TD><TD PORT="R"> ∅ </TD></TR>
         </TABLE>>,
         ]; 
}
1
```
"""
function surreal2dag(io::IO, x::SurrealFinite)
    println(io, "digraph \"", float(x), "\" {")
    k = 1
    SurrealsinDAG = Dict{SurrealFinite,Int}()
    m = surreal2dag_f(io, x, k, SurrealsinDAG)
    println(io, "}")
    return m
end
function surreal2dag_f(io::IO, x::SurrealFinite, k::Integer, SurrealsinDAG)
    m = k 
    if !haskey(SurrealsinDAG, x)
        SurrealsinDAG[x] = m 
        surreal2node(io, x, k) 
        c = 1
        for s in x.L
            if !haskey(SurrealsinDAG, s)
                m += 1
                # println(io, "   node_$k:L -> node_$m;")
                println(io, "   node_$k:\"" *  convert(String, s) * "," * string(c) * "\" -> node_$m [color=\"red3\"];")
                m = surreal2dag_f(io, s, m, SurrealsinDAG)
            else
                println(io, "   node_$k:\"" *  convert(String, s) * "," * string(c)  * "\" -> node_$(SurrealsinDAG[s]) [color=\"red3\"];") 
            end
            c += 1 
        end
        c = 1
        for s in x.R
            if !haskey(SurrealsinDAG, s)
                m += 1
                # println(io, "   node_$k:R -> node_$m;")
                println(io, "   node_$k:\"" *  convert(String, s) * "," * string(c) * "\" -> node_$m [color=\"blue3\"];")
                m = surreal2dag_f(io, s, m, SurrealsinDAG)
            else
                println(io, "   node_$k:\"" *  convert(String, s) * "," * string(c) * "\" -> node_$(SurrealsinDAG[s]) [color=\"blue3\"];")
            end
            c += 1 
        end
    end 
    return m 
end
surreal2dag(x::SurrealFinite) = surreal2dag(stdout, x)


""" 
    surreal2dot(x::SurrealFinite)
    surreal2dot(io::IO, x::SurrealFinite)

 Writes a surreal representation as a tree out in DOT format for drawing using GraphVis,
 and returns the number of nodes in the graph.

## Arguments
* `io::IO`: output stream, default is stdout
* `x::SurrealFinite`: the number to write out
    
## Examples
```jldoctest
julia> surreal2dot(convert(SurrealFinite, 1))
digraph "1.0" {
   node_1 [shape=none,margin=0,label=
         <<TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0" CELLPADDING="4">
         <TR><TD COLSPAN="2">1</TD></TR>
         <TR><TD PORT="L"> <TABLE BORDER="0" CELLBORDER="0" CELLPADDING="0"><TR><TD PORT="0,1"> 0 </TD> &nbsp; </TR></TABLE> </TD><TD PORT="R"> ϕ </TD></TR>
         </TABLE>>,
         ];
   node_1:"0,1" -> node_2;
   node_2 [shape=none,margin=0,label=<<B>0</B>>]
}
2
```
"""
function surreal2dot(io::IO, x::SurrealFinite)
    println(io, "digraph \"", float(x), "\" {")
    k = 1
    m = surreal2dot_f(io, x, k)
    println(io, "}")
    return m
end
function surreal2dot_f(io::IO, x::SurrealFinite, k::Integer)
    m = k
    if x == zero(x)
        println(io, "   node_$k [shape=none,margin=0,label=<<B>0</B>>]")
    else
        #if x.shorthand==""
        #    S = convert(String, x)
        #else
        #    S = x.shorthand
        #end
        surreal2node(io, x, k)
        c = 1
        for s in x.L
            m += 1
            # println(io, "   node_$k:L -> node_$m;")
            println(io, "   node_$k:\"" *  convert(String, s) * "," * string(c) * "\" -> node_$m;")
            m = surreal2dot_f(io, s, m)
            c += 1 
        end
        c = 1
        for s in x.R
            m += 1
            # println(io, "   node_$k:R -> node_$m;")
            println(io, "   node_$k:\"" *  convert(String, s) * "," * string(c)  * "\" -> node_$m;")
            m = surreal2dot_f(io, s, m)
           c += 1 
         end
    end 
    return m 
end
function surreal2node(io::IO, x::SurrealFinite, k::Integer; extra_args::String="")
    S = convert(String, x)
    if isempty(x.L)
        L = "∅" 
    else
        L = "<TABLE BORDER=\"0\" CELLBORDER=\"0\" CELLPADDING=\"0\"><TR>"
        c = 1
        for s in x.L
            tmp = convert(String, s)
            L *= "<TD PORT=\"$tmp," * string(c) * "\"> " * tmp  * " </TD> &nbsp; "
            c += 1
        end
        L *= "</TR></TABLE>" 
    end 
    if isempty(x.R)
        R = "∅" 
    else
        R = "<TABLE BORDER=\"0\" CELLBORDER=\"0\" CELLPADDING=\"0\"><TR>"
        c = 1
        for s in x.R
            tmp = convert(String, s)
            R *= "<TD PORT=\"$tmp," * string(c) * "\"> " * tmp * " </TD> &nbsp; "
            c += 1 
        end
        R *= "</TR></TABLE>" 
    end
    print(io, "   ")
    if k>=0
        label = "$k"
    else
        label = "m$k"
    end
    println(io, """
                node_$label [shape=none,margin=0,label=
                         <<TABLE BORDER=\"0\" CELLBORDER=\"1\" CELLSPACING=\"0\" CELLPADDING=\"4\">
                         <TR><TD COLSPAN=\"2\">$S</TD></TR>
                         <TR><TD PORT=\"L\"> $L </TD><TD PORT=\"R\"> $R </TD></TR>
                         </TABLE>>,$extra_args
                         ];""")
end
surreal2dot(x::SurrealFinite) = surreal2dot(stdout, x)

#######################################################

""" 
    generation(x::SurrealFinite)

 Finds the birthday of a surreal number, which is 1 + the max of any of its components.

## Arguments
* `x::SurrealFinite`: the number to operate on
     
## Examples
```jldoctest
julia> generation( convert(SurrealFinite, 1) )
1
``` 
"""
function generation(x::SurrealFinite)
    if x==zero(x)
        return 0
    else
        return max( maximum( generation.( [x.L; 0]) ),
                    maximum( generation.( [x.R; 0]) )) + 1
    end
end 

""" 
    canonicalise(s::SurrealFinite)

 Convert a surreal number form into its equivalent canonical form. It performs the 
 conversion by converting to a Rational and then back to a Surreal.  
     
## Examples 
```jldoctest
julia> convert(SurrealFinite, 1) - convert(SurrealFinite, 1)
{ { ϕ | { ϕ | ϕ } } | { { ϕ | ϕ } | ϕ } }
julia> pf( canonicalise( convert(SurrealFinite, 1) - convert(SurrealFinite, 1) ) )
{ ϕ | ϕ }
``` 
"""
canonicalise(s::SurrealFinite) = convert(SurrealFinite, convert(Rational, s))
iscanonical(s::SurrealFinite) = canonicalise(s) == s

function unique2!( X::Array{SurrealFinite} )
    # our own unique that is based on ≅
    sort!(X)
    i = 1
    while i < length(X) 
        if X[i] ≅ X[i+1]
            splice!( X, i )
        end
        i += 1
    end
end

###### standard math routines ##############################

sign(x::SurrealFinite) = x<zero(x) ? -one(x) : x>zero(x) ? one(x) : zero(x)
# abs(x::SurrealFinite) = x<zero(x) ? -x : x

# subtraction is much slower than comparison, so got rid of old versions
# these aren't purely surreal arithmetic, but everything could be, just would be slower
function floor(T::Type, s::SurrealFinite)
    if s < zero(s)
        if s >= -one(s)
            return -one(T)
        elseif s >= -2
            return convert(T, -2)
        else
            k = 1
            while s < -2^(k+1) && k < 12 # 12 is arbitrary, but it would be painful to go higher
                k += 1
            end 
            if k==12
                error("s is too large for the current floor function")
            end
            a = -2^(k+1)
            b = -2^k
            for i=1:k
                d = (a+b) / 2
                if s < d
                    b = d
                else
                    a = d
                end
            end
            return convert(T, a) # N.B. returns canonical form of floor 
        end
    elseif s < one(s) 
        return zero(T)
    elseif s < 2
        return one(T)
    else
        # start with geometric search to bound the number
        k = 1
        while s >= 2^(k+1) && k < 12 # 12 is arbitrary, but it would be painful to go higher
            k += 1
        end  
        if k==12
            error("s is too large for the current floor function")
        end
        # now a binary search to narrow it down
        a = 2^k
        b = 2^(k+1)
        for i=1:k
            d = (a+b) / 2
            if s < d
                b = d
            else
                a = d
            end
         end
        return convert(T, a) # N.B. returns canonical form of floor 
    end
end
function ceil(T::Type, s::SurrealFinite)
    if isinteger(s)
        return floor(T,s)
    else
        return floor(T,s) + 1
    end
end

# implicit rounding mode is 'RoundNearestTiesUp'
#   to bbe consistent, should do the other rounding modes, and a precision, but the latter is hard
function round(T::Type, s::SurrealFinite)
    return floor(T, s + 1//2)
end

trunc(T::Type, s::SurrealFinite) = s>=0 ? floor(T,s) : -floor(T,-s)

# simple versions 
floor(s::SurrealFinite) = floor(SurrealFinite, s)
ceil(s::SurrealFinite) = ceil(SurrealFinite, s)
round(s::SurrealFinite) = round(SurrealFinite, s)
trunc(s::SurrealFinite) = trunc(SurrealFinite, s)

# this should still be rewritten in terms of searches
function mod(s::SurrealFinite, n::SurrealFinite)
    if !isinteger(n)
        error("n should be an integer")
    end
    if s ≅ zero(s) 
        return s 
    elseif s < zero(s) 
        return mod(s + n, n)
    elseif s >= n
        return mod(s - n, n)
    else 
        return s
    end
end

function isinteger(s::SurrealFinite)
    if s ≅ floor(s)
        return true
    else
        return false
    end
end

isdivisible(s::SurrealFinite, n::SurrealFinite) = isinteger(s) ? mod(s,n) ≅ zero(s) : false 
isodd(s::SurrealFinite)  = isinteger(s) ? !isdivisible(s, convert(SurrealFinite,2) ) : false 
iseven(s::SurrealFinite) = isinteger(s) ?  isdivisible(s, convert(SurrealFinite,2) ) : false 
# ispow2

isinf(s::SurrealFinite) = false 
isnan(s::SurrealFinite) = false
isfinite(s::SurrealFinite) = true

parents(s::SurrealFinite) = [s.L; s.R]
function ⪯(s::SurrealFinite, t::SurrealFinite)# is s an ancestor of t, or equal to it
    if s == t
        return true 
    elseif iszero(s) 
        return true 
    else
        for p in parents(t)
            if s⪯p
                return true
            end
        end
    end
    return false
end
isancestor(s::SurrealFinite, t::SurrealFinite) = s ⪯ t && s!=t
≺(s::SurrealFinite, t::SurrealFinite) = isancestor(s, t)
≻(s::SurrealFinite, t::SurrealFinite) = t ≺ s
⪰(s::SurrealFinite, t::SurrealFinite) = t ⪯ s


############################################

# extra analysis functions

# structure to hold stats on a surreal form
struct SurrealDAGstats
    nodes::Integer # nodes in its DAG
    tree_nodes::Integer # nodes in a tree-based representation of the graph
    edges::Integer # edges in its DAG
    generation::Integer # generation/birthday of the surreal
    longest_path::Array{SurrealFinite,1} 
    paths::Int128 # number of paths from source to sink
    value::Rational
    minval::Rational 
    maxval::Rational
end
 
#    LP=true means keep track of a longest path
#    V =true means keep track of values
# have these as switches, because they are expensive to calculate on large DAGs 
function dag_stats(s::SurrealFinite, processed_list::Dict{SurrealFinite,SurrealDAGstats}; LP=false, V=true) 
    if s == zero(s)    
        nodes = 1  
        tree_nodes = 1
        edges = 0
        generation = 0
        longest_path = LP ? [zero(s)] : ϕ # important this be \phi and not an arbitrary empty array
        paths = 1
        value = 0 
        minval = 0  
        maxval = 0 
    else
        P = parents(s) 
        nodes = 1
        tree_nodes = 1
        edges = length(P)
        generation = 0
        longest_path = LP ? [zero(s)] : ϕ
        paths = 0 
        value = V ? convert(Rational, s) : 0
        minval = value 
        maxval = value
        for p in P 
            if !haskey(processed_list, p)
                (stats,l) = dag_stats(p, processed_list; LP=LP, V=V) 
                nodes += stats.nodes 
                edges += stats.edges 
                minval = min(minval, stats.minval)
                maxval = max(maxval, stats.maxval)
            else
                stats = processed_list[p]
            end 
            if stats.generation > generation
                generation = stats.generation 
                longest_path = LP ? copy(stats.longest_path) : ϕ
            end  
            tree_nodes += stats.tree_nodes
            paths += stats.paths
        end
        generation += 1 
        longest_path = LP ? [longest_path; s] : ϕ
    end  
    stats = SurrealDAGstats(nodes, tree_nodes, edges, generation, longest_path, paths, value, minval, maxval) 
    processed_list[s] = stats # could make this a reverse list, ...
    return ( stats, processed_list )  
end 
dag_stats(s::SurrealFinite; LP=false, V=true) = dag_stats(s, Dict{SurrealFinite,SurrealDAGstats}(); LP=LP, V=V)[1] 
nodes(s::SurrealFinite) = dag_stats(s).nodes 
edges(s::SurrealFinite) = dag_stats(s).edges
paths(s::SurrealFinite) = dag_stats(s).paths
tree_nodes(s::SurrealFinite) = dag_stats(s).tree_nodes 
function breadth(s::SurrealFinite) 
    d = dag_stats(s)
    return d.maxval - d.minval
end
function width(s::SurrealFinite) 
    # not yet defined exactly what I mean here --
    #    maybe maximum number of nodes for a fixed generation?
    0
end

# for debugging uniqueness of nodes in memory
#    check whether an entry in processed list is unique according to hash tables
#    the values in U should all be 1, but if we do say 2*4, without the caching in +,
#       we'll get a '2' because identical forms are being allocated multiple times in memory
function increment( U::Dict{UInt64, Int64}, h::UInt64 )
    if haskey(U, h)
        U[h] += 1
    else
        U[h] = 1
    end
end
function uniqueness(s::SurrealFinite,
                    P::Dict{UInt64, UInt64},
                    Q::Dict{UInt64, SurrealFinite},
                    U::Dict{UInt64, Int64},
                    V::Dict{UInt64, SurrealFinite} )
    h = hash(s)
    increment(U, h)
    @static if VERSION < v"0.7.0"
        os = object_id(s)
    else
        os = objectid(s)
    end 
    P[ os ] = h
    Q[ os ] = s
    V[h] = s
    for p in parents(s)
        @static if VERSION < v"0.7.0"
            op = object_id(p)
        else
            op = objectid(p)
        end  
        if !haskey( P, op )
            uniqueness(p, P, Q, U, V )
        end  
    end 
    return (P, Q, U, V ) 
end 
uniqueness(s::SurrealFinite) = uniqueness(s,
                                          Dict{UInt64, UInt64}(),
                                          Dict{UInt64, SurrealFinite}(),
                                          Dict{UInt64, Int64}(),
                                          Dict{UInt64, SurrealFinite}() )
uniqueness_max( U::Dict{UInt64, Int64} ) = maximum( values(U) )
function uniqueness_failure(  U::Dict{UInt64, Int64}, V::Dict{UInt64, SurrealFinite} )
    U2 = filter( (k,v) -> v>1 , U) # not Julia 1.0 compliant
    V2 = [V[i] for i in keys(U2) ]
    return V2
end

###########################################################
### work towards a formal DAG representation and tools

function surrealDAG(s::SurrealFinite)
    dag_dict = dag_stats(s, Dict{SurrealFinite,SurrealDAGstats}())[2] 
    # dag = Array{SurrealDAGstats}(length(dag_dict)) # sort from simplest to least simple
    dag = Array{SurrealFinite}(length(dag_dict)) # sort from simplest to least simple
    i = 1 
    for x in sort( collect(keys(dag_dict)); by = x->(dag_dict[x].generation,dag_dict[x].value,dag_dict[x].nodes) )
        dag[i] = x
        i += 1  
    end
    return dag # an array containing sorted list of all nodes of the surreal form
end
function dag_add(x::SurrealFinite, y::SurrealFinite)
    # bit of a cheat -- not reworking links, just relying on sub-sums
    # just intended to check the sums and cartesian products
    dag_x = surrealDAG(x)
    n_x = length(dag_x)
    dag_y = surrealDAG(y) 
    n_y = length(dag_y)
    dag_z = Array{SurrealFinite}( n_x * n_y )
    for i=1:n_x
        for j=1:n_y
            dag_z[(i-1)*n_y + j] = dag_x[i] + dag_y[j]
        end
    end
    return unique(dag_z)
end


