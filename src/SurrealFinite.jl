struct SurrealFinite <: Surreal
    shorthand::String
    X_L :: Array{SurrealFinite,1} 
    X_R :: Array{SurrealFinite,1} 
    # constructor should check that X_L <= X_R 
    function SurrealFinite(shorthand::String, X_L::Array{SurrealFinite}, X_R::Array{SurrealFinite})
        unique2!( X_L ) 
        unique2!( X_R )
        # X_L = sort( unique( X_L ) )
        # X_R = sort( unique( X_R ) )
        # X_L = isempty(X_L) ? ϕ : [maximum(X_L)]
        # X_R = isempty(X_R) ? ϕ : [minimum(X_R)]
        if X_L <= X_R
            return new(shorthand, X_L, X_R) 
        else 
            error("Must have X_L <= X_R")
        end 
    end
end 
SurrealFinite(shorthand::String, X_L::Array, X_R::Array ) =
    SurrealFinite( shorthand, convert(Array{SurrealFinite},X_L), convert(Array{SurrealFinite},X_R) )
SurrealFinite(X_L::Array, X_R::Array ) =
    SurrealFinite( "", convert(Array{SurrealFinite},X_L), convert(Array{SurrealFinite},X_R) )
≀(X_L::Array, X_R::Array ) = SurrealFinite(X_L::Array, X_R::Array )

# convert(::Type{Array{SurrealFinite}}, x::Array{Any}) = [ convert(SurrealFinite,x[i]) for i=1:length(x) ]
# convert(::Type{Array{Real}},    x::Array{Any}) = [ convert(Real,x[i]) for i=1:length(x) ]

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
        # integers
        return convert(SurrealFinite, convert(Int64,r) )
    elseif isinteger( log2(r.den) )
        # powers of 2 in denominator
        n = convert(Int64, round( log2(r.den) ))
        if abs(r) + n < 60 # 60 is a bit arbitrary -- could experiment more on real limits
            return SurrealFinite("$r", [convert(SurrealFinite, r - 1//2^n)], [convert(SurrealFinite, r + 1//2^n)] )
        else 
            error("generation too large")
        end 
    else 
        # error("we can't do these yet as they require infinite sets")
        return convert(SurrealFinite, r.num) // convert(SurrealFinite, r.den)
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
    else
        # do a binary search from top down, first valid is the one
        # xl = convert(Rational, maximum(s.X_L))
        # yl = convert(Rational, minimum(s.X_R))
        xl = maximum(s.X_L)
        xr = minimum(s.X_R)
        not_end = true
        k = 0
        a = 0
        b = 1
        while not_end && k < 24
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

function convert(::Type{AbstractFloat}, s::SurrealFinite )  
    return float( convert(Rational, s) )
end
 
promote_rule(::Type{Float64}, ::Type{SurrealFinite}) = SurrealFinite
promote_rule(::Type{Float32}, ::Type{SurrealFinite}) = SurrealFinite
# prolly need some more of these

ϕ = Array{SurrealFinite,1}(0) # empty array of SurrealFinites
zero(::SurrealFinite) = SurrealFinite("0", ϕ, ϕ )
one(::SurrealFinite) = SurrealFinite("1", [ zero(SurrealFinite) ], ϕ )
# ↑ = one(SurrealFinite)  # cause error???
# ↓ = -one(SurrealFinite) 
    
# relations
#    rewrite in terms of max/min, or really sup/inf
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
function <=(X::Array{SurrealFinite}, Y::Array{SurrealFinite} )
    if isempty(X) || isempty(Y)
        return true
    else
        for x in X
            for y in Y
                if ! (x <= y) 
                    return false
                end 
            end
        end
        return true
    end
end 
# ===(x::SurrealFinite, y::SurrealFinite) = x<=y && y<x # causes an error
≅(x::SurrealFinite, y::SurrealFinite) = x<=y && y<=x

≇(x::SurrealFinite, y::SurrealFinite) = !( x ≅ y ) 

==(x::SurrealFinite, y::SurrealFinite) = size(x.X_L) == size(y.X_L) &&
                             size(x.X_R) == size(y.X_R) &&
                             all(x.X_L .== y.X_L) &&
                             all(x.X_R .== y.X_R)  
<(x::SurrealFinite, y::SurrealFinite) = x<=y && !(y≅x)

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
function /(x::SurrealFinite)
    # not sure if I can do reciprocals with finite sequences
    
end

# binary operators
+(x::SurrealFinite, y::SurrealFinite) = SurrealFinite([x.X_L .+ y; x .+ y.X_L],
                                                      [x.X_R .+ y; x .+ y.X_R] )
-(x::SurrealFinite, y::SurrealFinite) = x + -y
function *(x::SurrealFinite, y::SurrealFinite)
    XL = x.X_L
    XR = x.X_R
    YL = y.X_L
    YR = y.X_R
    left = [ XL*y .+ x*YL .- XL*YL;
             XR*y .+ x*YR .- XR*YR
            ]
    right = [ XL*y .+ x*YR .- XL*YR;
              XR*y .+ x*YL .- XR*YL 
             ]
    return SurrealFinite("", left, right)
end
function *(x::SurrealFinite, Y::Array{SurrealFinite})
    return [ x*s for s in Y ]
end
function *(X::Array{SurrealFinite}, y::SurrealFinite)
    return [ s*y for s in X ]
end
function *(X::Array{SurrealFinite}, Y::Array{SurrealFinite})
    # return [ X*s for s in Y ]
    return vec([s*t for s in X, t in Y])
end

# /(x::SurrealFinite, y::SurrealFinite) = 1
# not sure if I can do division with finite sequences

function pf(x::SurrealFinite) 
    print("<", x.X_L, ":", x.X_R, ">") 
end
function pff(x::SurrealFinite)
    # tmp = isempty(x.X_L) ? ϕ : pff(x.X_L)
    print("<", pff.(x.X_L), ":", pff.(x.X_R), ">") 
end 
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


function generation(x::SurrealFinite)
    if x≅zero(x)
        return 0
    else
        gen = 0
        for z in x.X_L
            tmp = generation(z)
            if tmp > gen
                gen = tmp
            end
        end
        for z in x.X_R
            tmp = generation(z)
            if tmp > gen
                gen = tmp
            end
        end
        return gen+1
    end
end 

function simplify(s::SurrealFinite)
    # this is a bit of a cheap, but also simplify isn't unique otherwise
    return convert(SurrealFinite, convert(Rational, s))
end

sign(x::SurrealFinite) = x<zero(x) ? -one(x) : x>zero(x) ? one(x) : zero(x)
abs(x::SurrealFinite) = x<zero(x) ? -x : x

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

# old versions that required transform back to real
# isinteger(s::SurrealFinite) = isinteger( convert(AbstractFloat, s) )
# round(s::SurrealFinite) = convert(SurrealFinite, round( convert(AbstractFloat, s) ))

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

# not the more general form of rounding defined in Julia
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

