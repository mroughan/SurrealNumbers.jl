struct SurrealFinite <: Surreal
    shorthand::String
    X_L :: Array{SurrealFinite,1} 
    X_R :: Array{SurrealFinite,1} 
    # constructor should check that X_L < X_R
    function SurrealFinite(shorthand::String, X_L::Array{SurrealFinite}, X_R::Array{SurrealFinite})
        unique2!( X_L ) # make elements unique and sort them in increasing order
        unique2!( X_R ) # make elements unique and sort them in increasing order
        if X_L < X_R
            return new(shorthand, X_L, X_R)
        else 
            error("Surreal number must have X_L < X_R.")
        end
    end
end
SurrealFinite(shorthand::String, X_L::Array, X_R::Array ) =
    SurrealFinite( shorthand, convert(Array{SurrealFinite},X_L), convert(Array{SurrealFinite},X_R) )
SurrealFinite(X_L::Array, X_R::Array ) =
    SurrealFinite( "", convert(Array{SurrealFinite},X_L), convert(Array{SurrealFinite},X_R) )
≀(X_L::Array, X_R::Array ) = SurrealFinite(X_L::Array, X_R::Array )

hash(x::SurrealFinite, h::UInt) = hash( hash(x.X_L, h) * hash(x.X_R, h), h )
hash(X::Array{SurrealFinite}, h::UInt) = isempty(X) ? hash(0) : hash( hash(X[1]) * hash( X[2:end]), h )
hash(x::SurrealFinite) = hash(x, convert(UInt64,0) )
hash(X::Array{SurrealFinite}) = hash(X, convert(UInt64,0) )

function convert(::Type{SurrealFinite}, n::Int ) 
    if n==0 
        return SurrealFinite("0", ϕ, ϕ )
    elseif n>0
        return SurrealFinite(string(n), [convert(SurrealFinite, n-1)], ϕ )
    else
        return SurrealFinite(string(n), ϕ, [convert(SurrealFinite, n+1)] )
    end
    # should check to make sure abs(n) is not too big, i.e., causes too much recursion
end 
function convert(::Type{SurrealFinite}, r::Rational )  
    if isinteger(r) 
        return convert(SurrealFinite, convert(Int64,r) )
    elseif isinteger( log2(r.den) )
        # Non-integer dyadic numbers
        n = convert(Int64, round( log2(r.den) ))
        if abs(r) + n < 60 # 60 is a bit arbitrary -- could experiment more on real limits
            return SurrealFinite("$r", [convert(SurrealFinite, r - 1//2^n)], [convert(SurrealFinite, r + 1//2^n)] )
        else 
            error("Generation too large")
        end 
    else 
        error("we can't do these yet as they require infinite sets")
        # return convert(SurrealFinite, r.num) // convert(SurrealFinite, r.den)
    end  
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
        t = bits(f)[13:end]
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
        error(DomainError)
    end 
end 
function convert(::Type{Rational}, s::SurrealFinite )  
    if s ≅ zero(s)
        return 0 // 1
    elseif s ≅ one(s)
        return 1 // 1
    elseif s < zero(s)
        return -convert(Rational, -s)
    elseif s > one(s)
        convert(Rational, s - one(s) ) + 1 
    else # 0 < x < 1
        # do a binary search from top down, first valid is the simplest
        xl = maximum(s.X_L)
        xr = minimum(s.X_R)
        not_end = true
        k = 0
        a = 0
        b = 1
        while not_end && k < 24 # 24 is arbitrary, but it would be painful to go lower
            k += 1
            d = (b+a) // 2
            # print("k=$k, a=$a, b=$b, d=$d \n")

            c = convert(SurrealFinite,  d)
            if xl < c < xr
                return d
            elseif c <= xl
                a = d 
            elseif c >= xr
                b = d
            else
                error("this case should not happen")
            end 
        end
        error("should not get to here -- down to n/2^24")
    end
end

# some catch alls
convert(::Type{AbstractFloat}, s::SurrealFinite ) = float( convert(Rational, s) )
convert(::Type{Integer}, s::SurrealFinite ) = Int( convert(Rational, s) )
## convert{T<:Real}(::Type{T}, s::SurrealFinite ) =  convert(T, convert(Rational, s) )

# promote all numbers to surreals for calculations
promote_rule{T<:Real}(::Type{T}, ::Type{SurrealFinite}) = SurrealFinite

ϕ = Array{SurrealFinite,1}(0) # empty array of SurrealFinites
zero(::SurrealFinite) = SurrealFinite("0", ϕ, ϕ )
one(::SurrealFinite) = SurrealFinite("1", [ zero(SurrealFinite) ], ϕ )
# ↑ = one(SurrealFinite)  # this causes an error???
# ↓ = -one(SurrealFinite) 

# relations
#   these are written in terms of the definition, but could
#   rewrite in terms of max/min to make marginally faster
function <=(x::SurrealFinite, y::SurrealFinite)
    for t in x.X_L
        if y <= t
            return false
        end
    end
    for t in y.X_R
        if t <= x
            return false
        end
    end
    return true 
end
<(x::SurrealFinite, y::SurrealFinite) = x<=y && !(y≅x)
# ===(x::SurrealFinite, y::SurrealFinite) = x<=y && y<x # causes an error
≅(x::SurrealFinite, y::SurrealFinite) = x<=y && y<=x
≇(x::SurrealFinite, y::SurrealFinite) = !( x ≅ y ) 
==(x::SurrealFinite, y::SurrealFinite) = size(x.X_L) == size(y.X_L) &&
                                         size(x.X_R) == size(y.X_R) &&
                                         all(x.X_L .== y.X_L) &&
                                         all(x.X_R .== y.X_R)  

# comparisons between arrays are all-to-all, so don't have to be the same size
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
    if isempty(x.shorthand)
        SurrealFinite("", -x.X_R, -x.X_L )
    elseif x.shorthand == "0"
        zero(x)
    elseif x.shorthand[1] == '-'
        SurrealFinite(x.shorthand[2:end], -x.X_R, -x.X_L )
    else
        SurrealFinite("-"*x.shorthand, -x.X_R, -x.X_L )
    end 
end  
function /(x::SurrealFinite, y::SurrealFinite)
    xr = convert(Rational, x)
    yr = convert(Rational, y) 
    if y ≅ zero(y)
        error(InexactError)
    elseif isinteger(log2(abs(yr)))
        return x * convert(SurrealFinite,  1 // yr) 
    elseif isinteger(xr.num / yr)
        return convert(SurrealFinite,  (xr.num/yr) // xr.den) 
    else
        error(InexactError)        
    end
end

# binary operators
+(x::SurrealFinite, y::SurrealFinite) = SurrealFinite([x.X_L .+ y; x .+ y.X_L],
                                                      [x.X_R .+ y; x .+ y.X_R] )
# can't do like this because of empty arrays I think
#+(x::SurrealFinite, y::SurrealFinite) = SurrealFinite([x.X_L + y.X_L],
#                                                      [x.X_R + y.X_R] )
+(X::Array{SurrealFinite}, Y::Array{SurrealFinite}) = vec([s+t for s in X, t in Y])
 
-(x::SurrealFinite, y::SurrealFinite) = x + -y
-(X::Array{SurrealFinite}, Y::Array{SurrealFinite}) = X + -Y

function *(x::SurrealFinite, y::SurrealFinite)
    XL = x.X_L
    XR = x.X_R
    YL = y.X_L
    YR = y.X_R
    tmp1 = vec([s*y + x*t - s*t for s in XL, t in YL])
    tmp2 = vec([s*y + x*t - s*t for s in XR, t in YR])
    tmp3 = vec([s*y + x*t - s*t for s in XL, t in YR])
    tmp4 = vec([s*y + x*t - s*t for s in XR, t in YL])

    left  = [ tmp1; tmp2]
    right = [ tmp3; tmp4]
    return SurrealFinite("", left, right)
end
*(x::SurrealFinite, Y::Array{SurrealFinite}) = return [ x*s for s in Y ]
*(X::Array{SurrealFinite}, y::SurrealFinite) = y*X

# print commands
pf(x::SurrealFinite) = print("<", x.X_L, ":", x.X_R, ">") 
pff(x::SurrealFinite) = print("<", pff.(x.X_L), ":", pff.(x.X_R), ">") 
function show(io::IO, x::SurrealFinite)
    if x.shorthand != ""
        print_with_color(:bold, io, x.shorthand ) # could be :red
    else
        print( io, "<", x.X_L, ":", x.X_R, ">")
    end
end
function show(io::IO, X::Array{SurrealFinite})
    print(io, "{", join(X, ", "), "}")
end

# generation or birth day calculation
function generation(x::SurrealFinite)
    if x≅zero(x)
        return 0
    else
        return max( maximum( generation.( [x.X_L; 0]) ),
                    maximum( generation.( [x.X_R; 0]) )) + 1
    end
end 

# this is a bit of a cheat, but I'm not smart enough to work out how to do it otherwise
canonicalise(s::SurrealFinite) = convert(SurrealFinite, convert(Rational, s))

sign(x::SurrealFinite) = x<zero(x) ? -one(x) : x>zero(x) ? one(x) : zero(x)
# abs(x::SurrealFinite) = x<zero(x) ? -x : x

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

function isinteger(s::SurrealFinite)
    if s ≅ zero(s) 
        return true
    elseif s <= -one(s)
        return isinteger(s + one(s))
    elseif s >= one(s)
        return isinteger(s - one(s))
    else
        return false
    end
end

function floor(s::SurrealFinite)
    if zero(s) <= s < one(s) 
        return zero(s)
    elseif s < zero(s)
        return floor(s + one(s)) - one(s)
    elseif s >= one(s)
        return floor(s - one(s)) + one(s)
    end
end

function ceil(s::SurrealFinite)
    if zero(s) <= s < one(s) 
        return zero(s)
    elseif s < zero(s)
        return ceil(s + one(s)) - one(s)
    elseif s >= one(s)
        return ceil(s - one(s)) + one(s)
    end
end


# not the more general form of rounding defined in Julia -- should fix
function round(s::SurrealFinite)
    if s ≅ zero(s) 
        return s
    elseif s <= -one(s)
        return round(s + one(s)) - one(s)
    elseif s >= one(s)
        return round(s - one(s)) + one(s)
    elseif s >= convert(SurrealFinite, 1//2)
        return one(s)
    elseif s >= convert(SurrealFinite, -1//2)
        return zero(s)
    else
        return -ones(s)
    end
end

isinf(s::SurrealFinite) = false
isnan(s::SurrealFinite) = false
isfinite(s::SurrealFinite) = true

