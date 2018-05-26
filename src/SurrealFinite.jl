struct SurrealFinite <: Surreal
    shorthand::String 
    L::Array{SurrealFinite,1} 
    R::Array{SurrealFinite,1} 
    # constructor should check that L < R
    function SurrealFinite(shorthand::String, L::Array{SurrealFinite}, R::Array{SurrealFinite})
        # unique2!( L ) # make elements unique and sort them in increasing order
        # unique2!( R ) # make elements unique and sort them in increasing order
        L = sort(unique( L ))
        R = sort(unique( R ))
        # println("L = $L, R = $R")
        # use the fact they are sorted to not doa complete comparison
        if isempty(L) || isempty(R) || L[end] < R[1]
            return new(shorthand, L, R)
        else 
            error("Surreal number must have L < R.")
        end
    end 
end
SurrealFinite(shorthand::String, L::Array, R::Array ) =
    SurrealFinite( shorthand, convert(Array{SurrealFinite},L), convert(Array{SurrealFinite},R) )
SurrealFinite(L::Array, R::Array ) =
    SurrealFinite( "", convert(Array{SurrealFinite},L), convert(Array{SurrealFinite},R) )
≀(L::Array, R::Array) = SurrealFinite(L::Array, R::Array )

# put hash(x,0) into construct, based on recursive construction
# hash does a call to this, but hash(x.hash, h)
# keep the array versions
# note whether we get the expected performance improvements

hash(x::SurrealFinite, h::UInt) = hash( hash(x.L, h) * hash(x.R, h), h )
hash(X::Array{SurrealFinite}, h::UInt) = isempty(X) ? hash(0,h) : hash( hash(X[1],h) * hash( X[2:end],h), h )
#hash(X::Array{SurrealFinite}, h::UInt) = isempty(X) ? hash(0,h) : hash( prod(hash.(X, h )), h)
hash(x::SurrealFinite) = hash(x, convert(UInt64,0) )
hash(X::Array{SurrealFinite}) = hash(X, convert(UInt64,0) )

# global dictionary for use in reducing cost of repeated conversions back to real numbers
ExistingSurreals = Dict{SurrealFinite,Rational}()

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
    elseif ispow2(r.den)
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
    global ExistingSurreals
    if !haskey(ExistingSurreals, s)
        if s ≅ zero(s)
            ExistingSurreals[s] = 0 // 1
        elseif s ≅ one(s)
            ExistingSurreals[s] = 1 // 1
        elseif s < zero(s)
            ExistingSurreals[s] = -convert(Rational, -s)
        elseif (sf = floor(Integer, s)) ≅ s
            ExistingSurreals[s] = sf 
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
                    ExistingSurreals[s] = d
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
    return ExistingSurreals[s]
end

# some catch alls
convert(::Type{AbstractFloat}, s::SurrealFinite ) = float( convert(Rational, s) )
convert{T <: AbstractFloat}(::Type{T}, s::SurrealFinite ) = convert(T, convert(Rational, s) )
convert(::Type{Integer}, s::SurrealFinite ) = Int( convert(Rational, s) )
convert{T <: Integer}(::Type{T}, s::SurrealFinite ) = convert(T, convert(Rational, s) )
## convert{T<:Real}(::Type{T}, s::SurrealFinite ) =  convert(T, convert(Rational, s) )
function convert(::Type{String}, s::SurrealFinite ) 
    # try to work out a nice way to print it
    if isinteger(s)
        return string(convert(Integer, s))
    else
        r = convert(Rational, s)
        return string(r.num) * "/" * string(r.den)
    end
end

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
#   or in terms of set operations to make more succinct
function <=(x::SurrealFinite, y::SurrealFinite)
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
    if !isempty(x.L) && y <= x.L[end] 
        return false
    end
    if !isempty(y.R) && x >= y.R[1] 
        return false
    end
    return true 
end
<(x::SurrealFinite, y::SurrealFinite) = x<=y && !(y<=x)
# ===(x::SurrealFinite, y::SurrealFinite) = x<=y && y<x # causes an error
≅(x::SurrealFinite, y::SurrealFinite) = x<=y && y<=x
≅(x::Real, y::Real) = ≅(promote(x,y)...)
≇(x::SurrealFinite, y::SurrealFinite) = !( x ≅ y ) 
≇(x::Real, y::Real) = ≇(promote(x,y)...)
==(x::SurrealFinite, y::SurrealFinite) = size(x.L) == size(y.L) &&
                                         size(x.R) == size(y.R) &&
                                         all(x.L .== y.L) &&
                                         all(x.R .== y.R)  

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
#        if !(X[end] <= Y[1])
#            return false
#        end
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
#        # assume sort
#        if !(X[end] < Y[1])
#            return false
#        end
        return true
    end
end 

# unary operators
function -(x::SurrealFinite)
    if isempty(x.shorthand)
        SurrealFinite("", -x.R, -x.L )
    elseif x.shorthand == "0"
        zero(x)
    elseif x.shorthand[1] == '-'
        SurrealFinite(x.shorthand[2:end], -x.R, -x.L )
    else
        SurrealFinite("-"*x.shorthand, -x.R, -x.L )
    end 
end
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
        return convert(SurrealFinite,  (xr.num/yr) // xr.den) 
    else
        error(InexactError)        
    end
end

# binary operators
+(x::SurrealFinite, y::SurrealFinite) = SurrealFinite([x.L .+ y; x .+ y.L],
                                                      [x.R .+ y; x .+ y.R] )
# can't do like this because of empty arrays I think
#+(x::SurrealFinite, y::SurrealFinite) = SurrealFinite([x.L + y.L],
#                                                      [x.R + y.R] )
+(X::Array{SurrealFinite}, Y::Array{SurrealFinite}) = vec([s+t for s in X, t in Y])
 
-(x::SurrealFinite, y::SurrealFinite) = x + -y
-(X::Array{SurrealFinite}, Y::Array{SurrealFinite}) = X + -Y

function *(x::SurrealFinite, y::SurrealFinite)
    if x ≅ 0 || y ≅ 0
        return zero(x)
#    elseif x ≅ 1
#         return y
#    elseif y ≅ 1 
#        return x
    else
        # println("x = $x = ", float(x), ", y = $y = ", float(y))
        # print("   x = ")
        # pf(x)
        # println()
        # print("   y = ")
        # pf(y)
        # println()
        
        tmp1 = vec([s*y + x*t - s*t for s in x.L, t in y.L])
        # println("  tmp1 = $tmp1")
        tmp2 = vec([s*y + x*t - s*t for s in x.R, t in y.R])
        # println("  tmp2 = $tmp2")
        tmp3 = vec([s*y + x*t - s*t for s in x.L, t in y.R])
        # println("  tmp3 = $tmp3")
        tmp4 = vec([s*y + x*t - s*t for s in x.R, t in y.L])
        # println("  tmp4 = $tmp4")
        L = [ tmp1; tmp2]
        # println("    L = $L")
        # spf.(L)
        # println()
        R = [ tmp3; tmp4]
        # println("    R = $R")
        # spf.(R)
        # println()
        # println( " L < R ", L<R)
        return SurrealFinite("", L, R)
    end
end
*(x::SurrealFinite, Y::Array{SurrealFinite}) = return [ x*s for s in Y ]
*(X::Array{SurrealFinite}, y::SurrealFinite) = y*X


#####################################################3
# read/write routines for surreals

# print the first level in full (ignoring the top level shorthand if present)
#   these should be replaced using "expand"
pf(io::IO, x::SurrealFinite) = println(io, "{ ", x.L, " | ", x.R, " }")
pf(x::SurrealFinite) = pf(STDOUT, x) 

# expand at some level, 0 steps (if possible), 1 step, or completely
function expand(x::SurrealFinite; level=0)
    if level==0
        s = x.shorthand != "" ? x.shorthand : expand(x; level=1)
    elseif level==1
        tmpL = isempty(x.L) ? "ϕ" : join(convert.(String, x.L), ',')
        tmpR = isempty(x.R) ? "ϕ" : join(convert.(String, x.R), ',')
        s = "{ " * tmpL * " | " * tmpR  * " }"
        return s 
    elseif level>=1
        return "{ " * expand(x.L;level=level) * " | " * expand(x.R;level=level) * " }"
    end
end
expand(X::Array{SurrealFinite}; level=0) = isempty(X) ? "ϕ" : join(expand.(X; level=level), ',') 
 
# isnumeric(s::AbstractString) = ismatch(r"^-?\d*.\d*$", s) || ismatch(r"^-?\d*//-?\d*$", s)

# this has to parse a surreal written into a string 
function convert(::Type{SurrealFinite}, s::AbstractString ) 
    # interpret, (i) numbers as canonical, (ii) phi, \phi, ϕ correctly, (iii) structure
    s = replace(s, r"\s+", "") # remove white space
    s = replace(s, r"phi|\\phi|ϕ|\{\}", "") # replace phi or \phi with empty set
    # println("     s1 = $s")
    not_end = true
    basic_number = r"\{([^{}|]*)\|([^{|}]*)\}"
    while not_end
        if ismatch(basic_number, s)
            s = replace(s, basic_number, s"SurrealFinite([\1],[\2])") # remove comments
        else 
            not_end = false
        end
    end
    # println("     s2 = $s")
    return eval(parse(s))
end

# read in the full format
function read(io::IO, ::Type{SurrealFinite}, n::Int=1)
    X = Array{SurrealFinite}(n,)
    k = 1
    while !eof(io) && k<=n
        line = replace(readline(io), r"#.*", s"") # remove comments
        # println(" $k:   $line")
        if ismatch(r"\S", line)
            X[k] = convert(SurrealFinite, line)
            k += 1
        end
    end 
    return X
end
# read(filename::AbstractString, args...) = open(io->read(io, args...), filename)

# write out a string suitable for inclusion into latex docs
function surreal2tex(io::IO, x::SurrealFinite; level=0)
    s = expand(x; level=level)
    s = replace(s, r"\{", "\\\{") 
    s = replace(s, r"\|", "\\mid") 
    s = replace(s, r"\}", "\\\}")
    s = replace(s, r"(\d+)//(\d+)", s"\\frac{\1}{\2}")
    println(io,s)
end
surreal2tex(x::SurrealFinite; level=0) = surreal2tex(STDOUT, x; level=level)

# standard show will use shorthand when available
function show(io::IO, x::SurrealFinite)
    if io==STDOUT && x.shorthand != ""
        print_with_color(:bold, io, x.shorthand ) # could be :red
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
        print(io, "ϕ")
    else
        print(io, join(X, ", "))
    end
end
# special "canonicalised" output
spf(x::SurrealFinite) = print("{ ", canonicalise.(x.L), " | ", canonicalise.(x.R), " }")

function surreal2dag(io::IO, x::SurrealFinite)
    println(io, "digraph \"", float(x), "\" {")
    k = 1
    SurrealsinDAG = Dict{SurrealFinite,Int}()
    # for a in keys(ExistingSurreals) 
    #    delete!(ExistingSurreals, a)
    # end
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
                println(io, "   node_$k:\"" *  convert(String, s) * "," * string(c) * "\" -> node_$m;")
                m = surreal2dag_f(io, s, m, SurrealsinDAG)
            else
                println(io, "   node_$k:\"" *  convert(String, s) * "," * string(c)  * "\" -> node_$(SurrealsinDAG[s]);") 
            end
            c += 1 
        end
        c = 1
        for s in x.R
            if !haskey(SurrealsinDAG, s)
                m += 1
                # println(io, "   node_$k:R -> node_$m;")
                println(io, "   node_$k:\"" *  convert(String, s) * "," * string(c) * "\" -> node_$m;")
                m = surreal2dag_f(io, s, m, SurrealsinDAG)
            else
                println(io, "   node_$k:\"" *  convert(String, s) * "," * string(c) * "\" -> node_$(SurrealsinDAG[s]);")
            end
            c += 1 
        end
    end 
    return m 
end
surreal2dag(x::SurrealFinite) = surreal2dag(STDOUT, x)


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
    # L = isempty(x.L) ? "ϕ" : "" * join( convert.(String, x.L), ",</TD><TD> ") *"</TD>
    # R = isempty(x.R) ? "ϕ" : join( convert.(String, x.R), ", ")
    if isempty(x.L)
        L = "ϕ" 
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
        R = "ϕ" 
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
surreal2dot(x::SurrealFinite) = surreal2dot(STDOUT, x)

#######################################################

# generation or birth day calculation
function generation(x::SurrealFinite)
    if x==zero(x)
        return 0
    else
        return max( maximum( generation.( [x.L; 0]) ),
                    maximum( generation.( [x.R; 0]) )) + 1
    end
end 

# this is a bit of a cheat, but I'm not smart enough to work out how to do it otherwise
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

 
############################################

# extra analysis functions

# could do a better implementation of this
size_u(s::SurrealFinite) = length(unique(count_n(s)))
 
function size(x::SurrealFinite)
    if x==zero(x)
        return 1
    else
        return 1 + sum(size.(x.L)) + sum(size.(x.R))
    end
end
function n_zeros(x::SurrealFinite)
    if x==zero(x)
        return 1
    else
        return sum(n_zeros.(x.L)) + sum(n_zeros.(x.R))
    end
end 
function count_n(x::Surreal)
    # return [c; [convert(Rational, x)]; count_n.(x.L); count_n.(x.R) ] 
    tmp = [ convert(Rational, x) ]
    for s in x.L
        tmp = [tmp; count_n(s) ]
    end
    for s in x.R 
        tmp = [tmp; count_n(s) ]
    end
    # println( "x = ", convert(Rational, x))
    # println( join(tmp, ", ") )
    return tmp 
end
 
function depth(x::Surreal)
    if x==zero(x)
        return [0]
    else
        tmp = [ ]
        for s in x.L
            tmp = [tmp; depth(s) ]
        end
        for s in x.R
            tmp = [tmp; depth(s) ]
        end
        return 1 + tmp
    end
end
depth_max(x::SurrealFinite) = maximum( depth(x) )
depth_av(x::Surreal) = mean( depth(x) )





