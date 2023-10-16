const HashType = UInt64
mutable struct SurrealFinite <: Surreal
    # mutable means I can set the hash value when I need it, instead of
    #   (i) always having to calculate it recursively
    #   (ii) calculating it for even temporary surreals
    # but now version doesn't need it, its just faster for reasons???, and avoids some errors??? 
    shorthand::String 
    L::Vector{SurrealFinite}  
    R::Vector{SurrealFinite}
    maxL::Union{SurrealFinite,Missing}
    minR::Union{SurrealFinite,Missing}
    h::HashType # hash value, this is only set the first time the hash function is called
    # constructor should check that L < R
    function SurrealFinite(shorthand::String, L::Vector{SurrealFinite}, R::Vector{SurrealFinite}, 
                            maxL::Union{SurrealFinite,Missing}, minR::Union{SurrealFinite,Missing}, h::HashType)
        global Count
        global CountUncached
        global ExistingSurreals 
        Count['n'] += 1 
        inc_depth()
        if length(L) > 1
            L = sort( unique(L) )
            maxL = L[end]
            # L = unique(L)
            # maxL = maximum(L)
        elseif length(L) == 1
            maxL = L[1]
        else 
            maxL = missing
        end 
        if length(R) > 1
            R = sort( unique(R) )
            minR = R[1]
            # R = unique(R)
            # minR = minimum(R)
        elseif length(R) == 1
            minR = R[1]
        else
            minR = missing
        end 
        # println("L = $L, R = $R")
        # use the fact they are sorted to not do a complete comparison for <=
        # also means that == is much easier than if they were unsorted
        if isempty(L) || isempty(R) || maxL < minR
            h = hash(L, R)
            # tmp = new(shorthand, L, R, 0)
            # h_tmp = hash(tmp) # can't have hash value as an input, has to be calculated
            while haskey(ExistingSurreals, h) &&  
                (!quick_compare(ExistingSurreals[h].L, L) || 
                 !quick_compare(ExistingSurreals[h].R, R))
                # linear probing for a free slot in the Cache
                # ideally would compare the actual surreals, but the surreal hasn't been constructed yet
                # and a comparison of L (or R) won't work because it is order dependent, and we aren'T
                # deterministically sorting equal valued surreals
                Count['c'] += 1
                h += 1
                # hard to test this because collisions happen so rarely -- need a means to cause one
            end
            if !haskey(ExistingSurreals, h)
                # make sure every surreal we have calculated a hash for is registered
                CountUncached['n'] += 1
                ExistingSurreals[h] = new(shorthand, L, R, maxL, minR, h)
            else
                # if this version got created separately, then reuse the existing version so this is a no op           
            end
            dec_depth()
            return ExistingSurreals[h]
        else 
            error("Surreal number must have L < R (currently the extreme values are $(maxL) and $(minR) )") 
        end  
    end 
end
SurrealFinite(shorthand::String, L::AbstractArray, R::AbstractArray ) =
    SurrealFinite( shorthand, convert(Vector{SurrealFinite},L), convert(Vector{SurrealFinite},R), missing, missing, zero(HashType))
SurrealFinite( L::AbstractArray, R::AbstractArray ) = SurrealFinite( "", L, R)
≀(L::AbstractArray, R::AbstractArray) = SurrealFinite( L, R )

function quick_compare(A1::Vector{SurrealFinite}, A2::Vector{SurrealFinite})
    # try to do a quick comparison of two sets of surreals to see if they are equal
    # could create a == for vectors of surreals, but this is trying to work without a nest
    # so it can be called in a constructor
    if length(A1) != length(A2)
        return false
    elseif hash(A1) != hash(A2)
        return false
    else
        S1 = sort(A1, by=x->hash(x)) # sort just by hash values for unique ordering
        S2 = sort(A2, by=x->hash(x)) # sort just by hash values for unique ordering
        for (i,s) in enumerate(S1)
            if hash(s) != hash(S2[i]) # hash values of components will be unique
                return false
            end
        end
        return true
    end 
end 

const SurrealShort  = SurrealFinite 
const SurrealDyadic = SurrealFinite 

function hash(x::SurrealFinite)
    if x.h == 0
        x.h = hash(x.L, zero(UInt64) ) * hash(x.R, one(UInt64) )
        # this branch is only used on construction, in fact, not even there now
        # note we don't include shorthand or hash in calculation of hash
        #    -- they former doesn't matter, the latter would be circular
    end
    return x.h
end

function hash( L::Vector{SurrealFinite}, R::Vector{SurrealFinite})
    return hash(L, zero(HashType) ) * hash(R, one(HashType) )
end
hash(x::SurrealFinite, h::Integer) = hash(x, convert(HashType, h) )
function hash(X::Vector{SurrealFinite}) # effectively depth first???
    H = hash(zero(HashType))
    for x in X
        H ⊻= hash(x) # order should not be important, so don't do hash(x,H)
    end
    return H
end
hash(X::Vector{SurrealFinite}, h::UInt64) = hash( hash(X), h )
# function hash_old(X::Vector{SurrealFinite}, h::UInt64) # effectively breadth first???
#     if isempty(X)
#         hash(zero(UInt64), h) 
#     elseif length(X) == 1
#         hash(X[1], h) 
#     else
#         hash(X[1], h) * hash( X[2:end], h) # note, order is technically not important
#     end
# end
# really slow version: hash(X::Array{SurrealFinite}, h::UInt) = isempty(X) ? hash(0,h) : hash( prod(hash.(X, h )), h)

# global dictionaries to avoid repeated calculations of the same things
#   in particular to speed up calculations, but also to ensure that resulting
#   structures are actually DAGs, i.e., pointers to the same "number" point at the same bit of memory
# i guess eventually these should be replaced by let-closures: https://docs.julialang.org/en/v0.6.2/manual/variables-and-scoping/
# originally had things like 
#     ExistingProducts = Dict{SurrealFinite, Dict{SurrealFinite,SurrealFinite}}()
# but seem simpler to have an extra level of redirect to ensure fixed size, but maybe that is misguided?
const ExistingSurreals  = Dict{HashType, SurrealFinite}()
const ExistingConversions = Dict{HashType,Rational}()
const ExistingCanonicals = Dict{Rational,HashType}()
const ExistingCanonicalsC = Dict{HashType,Bool}()
# const ExistingLEQ       = Dict{HashType, Dict{HashType,Bool}}() 
const ExistingLEQ       = SwissDict{HashType, SwissDict{HashType,Bool}}() # small improvement from SwissDict, particullary on mem use
const ExistingLEQ1      = SwissDict{HashType, SwissDict{HashType,Bool}}() # small improvement from SwissDict, particullary on mem use
const ExistingLEQ2      = SwissDict{HashType, SwissDict{HashType,Int64}}() # small improvement from SwissDict, particullary on mem use
const ExistingLEQ3      = SwissDict{Tuple{HashType, HashType},Bool}() # flat version is nearly twice as big???
const ExistingLEQ4      = SwissDict{HashType, Set{HashType}}() # store LEQ/GT in least significant bit
const ExistingLEQ5      = SwissDict{HashType, SwissDict{HashType,Nothing}}() # store LEQ/GT in least significant bit
const ExistingEQ        = Dict{HashType, Dict{HashType,Bool}}() 
const ExistingProducts  = Dict{HashType, Dict{HashType,HashType}}()
# const ExistingSums      = Dict{HashType, Dict{HashType,HashType}}() 
const ExistingSums      = SwissDict{HashType, SwissDict{HashType,HashType}}() # small improvement obtained from SwissDict
const ExistingSums2     = SwissDict{HashType, HashType}() # flat version, not used at the moment
const ExistingSums3     = SwissDict{SurrealFinite, SwissDict{SurrealFinite,SurrealFinite}}() # for testing redirection
const ExistingSums4     = SwissDict{SurrealFinite, SwissDict{SurrealFinite,HashType}}() # for testing redirection
const ExistingNegations = Dict{HashType, HashType}() 
const ExistingFloors    = Dict{HashType, HashType}() 
const ExistingIntegers    = Dict{HashType, Bool}() 
const ExistingCanonicalIntegers  = Dict{HashType, Bool}() 

const Count         = Dict{Char, Integer}('+'=>0, '*'=>0, '-'=>0, 'n'=>0, 'c'=>0, '='=>0, '≤'=>0)
const CountUncached = Dict{Char, Integer}('+'=>0, '*'=>0, '-'=>0, 'n'=>0, 'c'=>0, '='=>0, '≤'=>0)
            # note that "uncached also excludes trivial cases like 0*x
function cache_hit_percent()
    global Count
    global CountUncached
    for k in keys(Count)
        hit_rate = (Count[k] - CountUncached[k]) / Count[k]
        println( " Hit % for operation $(k) is ", 100*hit_rate)
    end
end   

const RecursionDepth = Dict{String, Integer}("current_depth"=>0, "max_depth"=>0)
    # depth should never be <0, but use signed int for debugging
function inc_depth()
    global RecursionDepth
    RecursionDepth["current_depth"] += 1
    if RecursionDepth["current_depth"] > RecursionDepth["max_depth"]
        RecursionDepth["max_depth"] = RecursionDepth["current_depth"]
    end 
end
function dec_depth()
    global RecursionDepth
    RecursionDepth["current_depth"] -= 1
end 

# const ϕ = SVector{0,SurrealFinite}()
const ϕ = Array{SurrealFinite,1}(undef,0) # empty array of SurrealFinites
const ∅ = ϕ

const SurrealZero = SurrealFinite("0", ϕ, ϕ )
const SurrealOne  = SurrealFinite("1", [SurrealZero], ϕ ) 
const SurrealTwo  = SurrealFinite("2", [SurrealOne], ϕ ) 
const SurrealThree  = SurrealFinite("3", [SurrealTwo], ϕ )  
const SurrealMinusOne  = SurrealFinite("-1", ϕ, [SurrealZero] ) 
const SurrealMinusTwo  = SurrealFinite("-2", ϕ, [SurrealMinusOne] ) 
zero(::SurrealFinite) = SurrealZero # always use the same zero
one(::SurrealFinite)  = SurrealOne  # always use the same one
zero(::Type{T}) where T<:Surreal = SurrealZero # use the fast version, not Julia's default with conversion
one(::Type{T}) where T<:Surreal  = SurrealOne  # use the fast version, not Julia's default with conversion
# ↑ = one(SurrealFinite)  # this causes an error???
# ↓ = -one(SurrealFinite)  


const check_collision_flag = false # set this to be true to do a (slow) diagnostic check of hash collisions
# need to do a better sort here when this is set to true

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
    global ExistingCanonicalsC
    global ExistingLEQ
    global ExistingLEQ1
    global ExistingLEQ2
    global ExistingLEQ3
    global ExistingLEQ4
    global ExistingLEQ5
    # global ExistingGT
    global ExistingEQ
    global ExistingProducts
    global ExistingSums
    global ExistingSums2
    global ExistingSums3
    global ExistingSums4
    global ExistingNegations
    global ExistingFloors    
    global ExistingIntegers  
    global ExistingCanonicalIntegers  
    global Count
    global CountUncached
    global RecursionDepth
    global ExistingSurrealDAGstats

    # can't really empty this cache, or it breaks things, and doesn't reduce costs much anyway
    # empty!(ExistingSurreals)
    
    empty!(ExistingConversions)
    empty!(ExistingCanonicals)
    empty!(ExistingCanonicalsC)
    empty!(ExistingLEQ)
    empty!(ExistingLEQ1)
    empty!(ExistingLEQ2)
    empty!(ExistingLEQ3)
    empty!(ExistingLEQ4)
    empty!(ExistingLEQ5)
    # empty!(ExistingGT)
    empty!(ExistingEQ)
    empty!(ExistingProducts)
    empty!(ExistingSums)
    empty!(ExistingSums2)
    empty!(ExistingSums3)
    empty!(ExistingSums4)
    empty!(ExistingNegations)
    empty!(ExistingFloors)
    empty!(ExistingIntegers)
    empty!(ExistingCanonicalIntegers)
    empty!(ExistingSurrealDAGstats)
    
    reset!(Count)
    reset!(CountUncached)
    reset!(RecursionDepth)
    return 1
end

function cache_stats_summary( )
    # this function is really just a reminder of how to do this
    df = Main.varinfo( SurrealNumbers, r"Existing*" )
    return df
end 

function cache_stats( C::AbstractDict{S, T} ) where {T <: Any, S <: Number}
# function cache_stats( C::Dict{S, T} ) where T <: Any where S <: Number # also works
    total_bytes1 = summarysize(C)
    total_entries1 = length(C)
    total_bytes2 = 0
    total_entries2 = 0
    for k in keys(C)
        total_bytes2 += summarysize( C[k] )
        total_entries2 += length( C[k] )     
    end 
    return total_bytes1, total_bytes2, total_entries1, total_entries2
end 

function cache_stats_LEQcount( )
    C = ExistingLEQ2
    m = 1000
    x = collect( -m: 1 : m )
    y = zeros(Int, size(x) )
    for k1 in keys(C)
        for k2 in keys( C[k1] )
            y[ C[k1][k2] + m + 1 ] += 1
        end
    end
    k = findall(y .> 0)
    return x[k], y[k]
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
    ExistingCanonicals[Rational(n)] = hr 
    # these would save time later, but seem like cheating
    #   ExistingIntegers[hr] = true
    #   ExistingFloors[hr] = n
    #
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
            r_l = r - 1//2^n
            r_r = r + 1//2^n
            if r_l >= r_r
                error("   weirdness: $r -> ($r_l, $r_r)")
            end
            result = SurrealFinite("$r", [convert(SurrealFinite, r_l)], [convert(SurrealFinite, r_r)] )
        else 
            error("Generation too large")
        end 
    else 
        error("we can't do these yet as they require infinite sets")
        # return convert(SurrealFinite, r.num) // convert(SurrealFinite, r.den)
    end  
    hr = hash(result)
    ExistingCanonicals[r] = hr 
    # these would save time later, but seem like cheating
    #   ExistingIntegers[hr] = isinteger(r)
    #   ExistingFloors[hr] = floor(r)
    # 
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
            throw(DomainError(f,"the value must be finite"))
        end
    end 
end 
dali(x) = convert(SurrealFinite, x)
SurrealFinite(x) = convert(SurrealFinite, x)

function convert(::Type{T}, s::SurrealFinite ) where T <: Integer
    global ExistingConversions
    h = hash(s)
    if haskey(ExistingConversions, h)
        return convert(T, ExistingConversions[h])
    elseif iszero(s)
        result = 0
    elseif equivtozero(s)
        result = 0
    elseif s < 0 
        result = -convert(T, -s)
    elseif isinteger(s) # implicitly > 0
        f = floor(s) # this returns a canonical form of the integer
        result = convert(T, f.maxL) + 1
        # recursive dive through all left sets might not fast for big values, except these will each be quick
        # and all prior isinteger and floors will already be calculated by first recursive calls to these fns
    else    
        throw( InexactError( :convert, T, s) )
    end
    ExistingConversions[h] = convert(Rational, result)
    # ExistingConversions[hash(-s)] = -result 
    return result
end 

function convert(::Type{Rational}, s::SurrealFinite ) 
    global ExistingConversions
    h = hash(s)
    if haskey(ExistingConversions, h)
        return ExistingConversions[h]
    else 
        if equivtozero(s)
            result = 0 // 1
        elseif s ≅ one(s)
            result = 1 // 1
        elseif s < zero(s)
            result = -convert(Rational, -s)
        elseif isinteger(s)
            result = Rational( convert(Integer, s )  )
        else
            # do a binary search from top down, first valid is the simplest
            sf = floor(s)
            xl = isempty(s.L) ? sf   : s.maxL
            xr = isempty(s.R) ? sf+1 : s.minR 
            not_end = true
            k = 0
            a = convert( Integer, sf )
            b = a + 1
            while not_end && k < 24 # 24 is arbitrary, but it would be painful to go lower
                k += 1
                d = (b+a) / 2 # a and b are real numbers, so we can do division easily
                              # but also d will be the value we return, so must be a Rational
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
                    error("this case should never happen")
                end 
            end
            if k >= 24
                error("number's exponent too large to convert (ie exponent <= -24)")
            end
        end 
    end 
    ExistingConversions[h] = result
    # ExistingConversions[hash(-s)] = -result 
    return result
end


# some catch alls
convert(::Type{T}, s::SurrealFinite ) where {T <: AbstractFloat} = convert(T, convert(Rational, s) )
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
function overbarString( s::SurrealFinite )
    # use this sparingly
    if iscanonical(s)
        #return "<SPAN STYLE=\"text-decoration:overline\">" * convert(String, s) * "</SPAN>"
        # return convert(String, s) * "\u0305"
        return "<o>" * convert(String, s) * "</o>"
    else
        return convert(String, s)
    end
end

# promote all numbers to surreals for calculations
promote_rule(::Type{T}, ::Type{SurrealFinite}) where {T<:Real} = SurrealFinite

# relations: the latest version uses the fact that letf and right sets are sorted
function leq_without_cache(x::SurrealFinite, y::SurrealFinite)
    # slightly misnamed -- it can use the cache for recursed calculations, just not at the top level
    global CountUncached
    CountUncached['≤'] += 1
    if !isempty(x.L) && leq(y, x.maxL)
        return false
    elseif !isempty(y.R) && leq(y.minR, x)
        return false
    else
        return true
    end
end
function leq_0(x::SurrealFinite, y::SurrealFinite)
    # simple direct cache
    global Count
    Count['≤'] += 1
    global ExistingLEQ
    hx = hash(x) 
    hy = hash(y)
    if !haskey(ExistingLEQ, hx)
        ExistingLEQ[hx] = Dict{HashType,Bool}()
    end         
    if !haskey(ExistingLEQ[hx], hy)
        ExistingLEQ[hx][hy] = leq_without_cache(x, y)
    end 
    return ExistingLEQ[hx][hy]
end
function leq_02(x::SurrealFinite, y::SurrealFinite)
    # simple direct cache
    global Count
    Count['≤'] += 1
    global ExistingLEQ3
    hx = hash(x) 
    hy = hash(y)        
    if !haskey(ExistingLEQ3, (hx,hy))
        ExistingLEQ3[(hx,hy)] = leq_without_cache(x, y)
    end 
    return ExistingLEQ3[(hx,hy)]
end
function leq_01(x::SurrealFinite, y::SurrealFinite)
    # simple direct cache
    global Count
    Count['≤'] += 1
    global ExistingLEQ
    hx = hash(x) 
    hy = hash(y)
    if haskey(ExistingLEQ, hx) && haskey(ExistingLEQ[hx], hy)
        return ExistingLEQ[hx][hy]
    end

    # start by checking integer approximations
    if haskey(ExistingFloors, hx) && haskey(ExistingFloors, hy)
        if ExistingFloors[hx] < ExistingFloors[hy]
            return true
        elseif ExistingFloors[hx] > ExistingFloors[hy]
            return false
        end
    end

    # if floors aren't already calculated, or fx==fy, then do a more accurate comparison, and only cache these cases
    if !haskey(ExistingLEQ, hx)
        ExistingLEQ[hx] = Dict{HashType,Bool}()
    end         
    if !haskey(ExistingLEQ[hx], hy)
        ExistingLEQ[hx][hy] = leq_without_cache(x, y)
    end 
    return ExistingLEQ[hx][hy]
end
function leq_1(x::SurrealFinite, y::SurrealFinite)
    # this version uses a second, hopefully smaller cache that might speed things up (but doesn't)
    #   but it really needs that the big cache be pruned at some point, because memory efficiency is the big deal here
    global Count
    Count['≤'] += 1
    global ExistingLEQ
    global ExistingLEQ1
    hx = hash(x)
    hy = hash(y)
    if haskey(ExistingLEQ1, hx) && haskey(ExistingLEQ1[hx], hy)
        return ExistingLEQ1[hx][hy]
    end
    if !haskey(ExistingLEQ, hx)
        ExistingLEQ[hx] = Dict{HashType,Bool}()
    end
    if !haskey(ExistingLEQ1, hx)
        ExistingLEQ1[hx] = Dict{HashType,Bool}()
    end
    if !haskey(ExistingLEQ[hx], hy)
        ExistingLEQ[hx][hy] = leq_without_cache(x, y)
    else
        ExistingLEQ1[hx][hy] = ExistingLEQ[hx][hy]
    end 
    return ExistingLEQ[hx][hy]
end
function leq_2(x::SurrealFinite, y::SurrealFinite)
    # this version counts repeated hits on the cache
    global Count
    Count['≤'] += 1
    global ExistingLEQ
    global ExistingLEQ2
    hx = hash(x)
    hy = hash(y)
    if haskey(ExistingLEQ2, hx) && haskey(ExistingLEQ2[hx], hy)
        ExistingLEQ2[hx][hy] += sign( ExistingLEQ2[hx][hy] )
        if ExistingLEQ2[hx][hy] > 0
            return true
        elseif ExistingLEQ2[hx][hy] < 0
            return false
        else
            error("The second level cache should never be zero")
        end
    end
    if !haskey(ExistingLEQ, hx)
        ExistingLEQ[hx] = Dict{HashType,Bool}()
    end         
    if !haskey(ExistingLEQ2, hx)
        ExistingLEQ1[hx] = Dict{HashType,Int}()
    end
    if !haskey(ExistingLEQ[hx], hy)
        ExistingLEQ[hx][hy] = leq_without_cache(x, y)
    else
        if ExistingLEQ[hx][hy]
            ExistingLEQ2[hx][hy] = 1
        else
            ExistingLEQ2[hx][hy] = -1
        end
    end 
    return ExistingLEQ[hx][hy]
end
function leq_4(x::SurrealFinite, y::SurrealFinite)
    # use 2-layer nest Dict/Set with least sig bit of object in set providing T/F
    global Count
    Count['≤'] += 1
    global ExistingLEQ4
    hx = hash(x) 
    hy = hash(y)
    if !haskey(ExistingLEQ4, hx)
        ExistingLEQ4[hx] = Set{HashType}()
    end
    if hy in ExistingLEQ4[hx]
        result = true
    elseif hy+1 in ExistingLEQ4[hx]
        result = false
    else
        result = leq_without_cache(x, y)
        if result
            push!( ExistingLEQ4[hx], hy )   
        else
            push!( ExistingLEQ4[hx], hy+1 ) # could introduce yet another sort of hash collision :(, but with low probability
        end
    end
    return result
end
function leq_5(x::SurrealFinite, y::SurrealFinite)
    # use 2-layer nest Dict/Set with least sig bit of object in set providing T/F
    global Count
    Count['≤'] += 1
    global ExistingLEQ5
    hx = hash(x) 
    hy = hash(y)
    if !haskey(ExistingLEQ5, hx)
        ExistingLEQ5[hx] = Dict{HashType, Nothing}()
    end
    if haskey( ExistingLEQ5[hx], hy)
        result = true
    elseif haskey(ExistingLEQ5[hx], hy+1)
        result = false
    else
        result = leq_without_cache(x, y)
        if result
            ExistingLEQ5[hx][hy] = nothing  
        else
            ExistingLEQ5[hx][hy+1] = nothing # could introduce yet another sort of hash collision :(, but with low probability
        end
    end
    return result
end

leq(x::SurrealFinite, y::SurrealFinite) = (inc_depth(); b=leq_5(x,y); dec_depth(); return b)

# using SparseArrays doesn't seem to work as there is a massive hit to create a big enough sparse array even though empty
# const ExistingLEQ1 = spzeros(Bool, 2^16, 2^16)
# const ExistingGT1  = spzeros(Bool, 2^16, 2^16)

# using sets is slower, and uses more memory ????
# const ExistingLEQ2 = Dict{HashType, Set{HashType}}()
# const ExistingGT = Dict{HashType, Set{HashType}}()
# function leq2(x::SurrealFinite, y::SurrealFinite)
#     global Count
#     Count['≤'] += 1
#     global ExistingLEQ2
#     global ExistingGT
#     hx = hash(x)
#     hy = hash(y)
#     if !haskey(ExistingLEQ2, hx)
#         ExistingLEQ2[hx] = Set{HashType}()
#     end         
#     if !haskey(ExistingGT, hx)
#         ExistingGT[hx] = Set{HashType}()
#     end
#     if in(hy, ExistingLEQ2[hx] )
#         return true
#     elseif in(hy, ExistingGT[hx] )
#         return false
#     else
#         result = leq_without_cache(x, y)
#     end 
#     if result
#         push!( ExistingLEQ2[hx], hy)
#     else    
#         push!( ExistingGT[hx], hy)
#     end 
#     return result
# end

<=(x::SurrealFinite, y::SurrealFinite) = leq(x, y)
<(x::SurrealFinite, y::SurrealFinite) = !(y<=x)
#### can't us ≡ because it is used for ===
# ===(x::SurrealFinite, y::SurrealFinite) = x<=y && y<x # causes an error
#     === is 'egal', and hardcoded for mutables to test they are same object in memory
≅(x::SurrealFinite, y::SurrealFinite) = x<=y && y<=x
≅(x::Real, y::Real) = ≅(promote(x,y)...)
≇(x::SurrealFinite, y::SurrealFinite) = !( x ≅ y ) 
≇(x::Real, y::Real) = ≇(promote(x,y)...)


# sort( X::Array{SurrealFinite} ) = sort( X, lt = (x,y) -> x == y ? hash(x)<hash(y) : x<y )
# sort( X::Array{SurrealFinite} ) = sort( X, lt = (x,y) -> x ≅ y ? hash(x)<hash(y) : x<y )
# NB do it this way because using "by =" implies use of equals signs in sort, and we use equals for identity, not equal value
#
# old reasoning was that I needed a deterministic sorted list to be able to compare "equal" things
# but sorting takes reassignments, which take time (lots it seems)
# so rather, when we create a new surreal in (most, but check all) ways we can creaete it, we check its hash
# and if we find it already has a form, then we use that, so the sorted order is based on the first created version
# this makes the results of code non-isolated -- the output of a result can depend on what has gone before
# but within one session everything will be consistent.
# just be aware, data written by one sequence of calculations can then be different the next time
# but to be fair, hashes are not stable between versions of Julia, so relying on the hash for a deterministic sort
# was always a little troublesome.
#
# But we do need the correct sorted versions when including a check for collisions, or we see artificial collisions
#
# anyway, upshot is that
#   1. we use the default sort
#   2. elements of left and right sets are sorted by value, and equal value items might appear in any order, but should be consistent in a session
#       unless deliberately broken (CHECK)
#   3. hash of a vector has to be order independent (because it is really a set)
# 

# old, slow, direct version
# ==(x::SurrealFinite, y::SurrealFinite) = size(x.L) == size(y.L) &&
#                                         size(x.R) == size(y.R) &&
#                                         all(x.L .== y.L) &&
#                                         all(x.R .== y.R)
# 
function equals(x::SurrealFinite, y::SurrealFinite)
    global ExistingEQ
    global Count
    global CountUncached
    Count['='] += 1
    inc_depth()
    hx = hash(x)
    hy = hash(y)
    # if hx == hy
    #     return true # trust that hashes don't collide
    # else    
    #    return false
    # end 
    if !haskey(ExistingEQ, hx)
        ExistingEQ[hx] = Dict{HashType,Bool}()
    end         
    # if !haskey(ExistingEQ, hy)
    #     ExistingEQ[hy] = Dict{HashType,Bool}() # it may be commutative, but when =, it doesn't matter
    # end
    if !haskey(ExistingEQ[hx], hy)
        CountUncached['='] += 1
        if hx != hy
            ExistingEQ[hx][hy] = false
            # ExistingEQ[hy][hx] = false
        elseif length(x.L) != length(y.L) || length(x.R) != length(y.R)
            ExistingEQ[hx][hy] = false
            # ExistingEQ[hy][hx] = false
        else
            for i=1:length(x.L)
                if !equals(x.L[i], y.L[i])
                    ExistingEQ[hx][hy] = false
                    # ExistingEQ[hy][hx] = false
                    dec_depth()
                    return ExistingEQ[hx][hy]
                end
            end
            for i=1:length(x.R)
                if !equals(x.R[i], y.R[i])
                    ExistingEQ[hx][hy] = false 
                    # ExistingEQ[hy][hx] = false
                    dec_depth()
                    return ExistingEQ[hx][hy]                    
                end 
            end 
            ExistingEQ[hx][hy] = true
            if !haskey(ExistingLEQ, hx)
                ExistingLEQ[hx] = Dict{HashType,Bool}()
            end         
            ExistingLEQ[hx][hy] = true # we get this for free, seems like lookup in LEQ would save memory, but doesn't save time :(
            # ExistingEQ[hy][hx] = true
        end
    end 
    dec_depth()
    return ExistingEQ[hx][hy]
end
==(x::SurrealFinite, y::SurrealFinite) = x===y || equals(x, y)
iszero(x::SurrealFinite) = isempty(x.L) && isempty(x.R) # tests for canonical zero 
equivtozero(x::SurrealFinite) = x ≅ zero(x)             # tests equiv to zero 

# comparisons between sets (i.e., arrays) are all-to-all, so
#   (1) don't have to be the same size
#   (2) the < is not exactly the same as the < defined above
# but these shouldn't be used for left- and right set comparisons because these sets are sorted
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
    global CountUncached
    Count['-'] += 1
    inc_depth()
    hx = hash(x) # build the dictionary in terms of hashs, because it is used quite a bit
    if iszero(x) 
        result = x
    elseif haskey(ExistingNegations, hx)
        dec_depth()
        return ExistingSurreals[ ExistingNegations[hx] ]
    elseif isempty(x.shorthand)
        result = SurrealFinite("", -x.R, -x.L )
    elseif x.shorthand[1] == '-'
        result = SurrealFinite(x.shorthand[2:end], -x.R, -x.L )
    else
        result = SurrealFinite("-"*x.shorthand, -x.R, -x.L )
    end  
    hr = hash(result)
    # if haskey(ExistingSurreals, hr)
    #     result = ExistingSurreals[hr] # don't double up on memory
    # else
    #     ExistingSurreals[hr] = result
    # end
    # if !haskey(ExistingSurreals, hx)
    #     ExistingSurreals[hx] = x
    # end 
    ExistingNegations[hx] = hr
    ExistingNegations[hr] = hx
    CountUncached['-'] += 1
    dec_depth()
    return result
end

# # binary operators

# function non_hierarchical_cache_sum(x::SurrealFinite, y::SurrealFinite)
#     # this is used only to test whether the hierarchical cache actually improved performance
#     global ExistingSurreals 
#     global ExistingSums2
#     global Count
#     global CountUncached
#     Count['+'] += 1
#     hx = hash(x) # build the dictionary in terms of hashs, because it is used quite a bit
#     hy = hash(y) #   and these amortise the cost of initial calculation of hashs
#     h = hx * hy
#     if haskey(ExistingSums2, h) 
#         return ExistingSurreals[ ExistingSums2[h] ]
#     elseif iszero(x)
#         result = y   
#     elseif iszero(y)
#         result = x 
#     else 
# #       result = SurrealFinite( [x.L .+ y; x .+ y.L],
#         #                               [x.R .+ y; x .+ y.R] )
#         shorthand = "($(x.shorthand) + $(y.shorthand))"
#         result = SurrealFinite(shorthand,
#                                [[x + y for x in x.L]; [x + y for y in y.L]],
#                                [[x + y for x in x.R]; [x + y for y in y.R]] )
#     end
#     hr = hash(result)
#     if haskey(ExistingSurreals, hr)
#        result = ExistingSurreals[hr] # don't double up on memory, make sure we always point at the existing surreal form
#     else
#        ExistingSurreals[hr] = result
#     end
#     ExistingSums2[h] = hr
#     CountUncached['+'] += 1
#     return result
# end


function +(x::SurrealFinite, y::SurrealFinite)
    global ExistingSurreals 
    global ExistingSums
    global Count
    global CountUncached
    global check_collision_flag
    commutative = false # set to true to store y+x, everytime we calculate x+y
    Count['+'] += 1
    inc_depth()
    hx = hash(x) # build the dictionary in terms of hashs, we retain some control over them
    hy = hash(y) # which has been useful, but at some point, should probably drop the indirection

    if haskey(ExistingSums, hx) && haskey(ExistingSums[hx], hy)
        dec_depth()
        return ExistingSurreals[ ExistingSums[hx][hy] ]
    elseif !commutative && haskey(ExistingSums, hy) && haskey(ExistingSums[hy], hx)
        dec_depth()
        return ExistingSurreals[ ExistingSums[hy][hx] ]
    elseif iszero(x)
        dec_depth()
        return y
    elseif iszero(y)
        dec_depth()
        return x
    else 
        shorthand = "($(x.shorthand) + $(y.shorthand))"
        # trying to create the sum in close to sorted order didn't help - just wasted an extra op.
        #if x <= y
            result = SurrealFinite(shorthand,
                    [[x_l + y for x_l in x.L]; [x + y_l for y_l in y.L]],
                    [[x_r + y for x_r in x.R]; [x + y_r for y_r in y.R]] )
        #else
        #    result = SurrealFinite(shorthand,
        #            [[x + y_l for y_l in y.L]; [x_l + y for x_l in x.L]],
        #            [[x + y_r for y_r in y.R]; [x_r + y for x_r in x.R]] )
        #end
        #       result = SurrealFinite( [x.L .+ y; x .+ y.L],   # dot syntax was a fair but slower than comprehensions
        #                               [x.R .+ y; x .+ y.R] )
    end
    CountUncached['+'] += 1
    hr = hash(result)
    # if haskey(ExistingSurreals, hr)
    #     if check_collision_flag
    #         if ExistingSurreals[hr] == result
    #             # only do this check when debugging because we also need to change the sort to a more expensive version 
    #             result = ExistingSurreals[hr]
    #             # don't double up on memory, reuse existing surreal by changing reference here
    #             # note that different additions could result in the same surreal form, and we don't want to have a difference
    #             # version for each way of creating this
    #         else
    #             error("HASH collision :( -- ($hx,$hy,$hr,$hash(result)) $( (pf(ExistingSurreals[hx]), pf(ExistingSurreals[hy]), pf(ExistingSurreals[hr]), pf(result)) )")
    #         end 
    #     else
    #         result = ExistingSurreals[hr]
    #     end
    # else
    #     ExistingSurreals[hr] = result
    # end

    if !haskey(ExistingSums, hx)
        ExistingSums[hx] = Dict{HashType,HashType}()
        # ExistingSums3[x] = Dict{SurrealFinite,SurrealFinite}()
        # ExistingSums4[x] = Dict{SurrealFinite,HashType}()
    end
    ExistingSums[hx][hy] = hr
    # ExistingSums3[x][y] = result
    # ExistingSums4[x][y] = hr
     
    if commutative
        if !haskey(ExistingSums, hy)
            ExistingSums[hy] = Dict{HashType,HashType}() 
        end
        ExistingSums[hy][hx] = hr
    end
    dec_depth()
    return result
end

# can't do like this because of empty arrays I think
#     +(x::SurrealFinite, y::SurrealFinite) = SurrealFinite([x.L + y.L],[x.R + y.R] )
# but also would fail to do caching

+(X::Array{SurrealFinite}, Y::Array{SurrealFinite}) = vec([s+t for s in X, t in Y])

-(x::SurrealFinite, y::SurrealFinite) = x + -y
-(X::Array{SurrealFinite}, Y::Array{SurrealFinite}) = X + -Y

function *(x::SurrealFinite, y::SurrealFinite)
    global ExistingSurreals 
    global ExistingProducts
    global Count
    global CountUncached
    Count['*'] += 1 
    inc_depth()
    hx = hash(x)
    hy = hash(y)

    if haskey(ExistingProducts, hx) && haskey(ExistingProducts[hx], hy)
        dec_depth()
        return ExistingSurreals[ ExistingProducts[hx][hy] ]
    elseif iszero(x) || iszero(y) 
        dec_depth()
        return SurrealZero
    elseif x == SurrealOne 
        dec_depth()
        return y 
    elseif y == SurrealOne
        dec_depth()
        return x
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
        ExistingProducts[hx] = Dict{HashType,HashType}()
    end
    if !haskey(ExistingProducts, hy)
        ExistingProducts[hy] = Dict{HashType,HashType}()
    end 
    hr = hash(result)
    # if haskey(ExistingSurreals, hr)
    #     result = ExistingSurreals[hr] # don't double up on memory
    # else 
    #     ExistingSurreals[hr] = result
    # end
    ExistingProducts[hx][hy] = hr
    ExistingProducts[hy][hx] = hr
    CountUncached['*'] += 1 
    dec_depth()
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


# somewhat incomplete (for reasons described), its a cheat, and potentially very, very slow
# but it can handle some obvious cases like dividing by 2
function /(x::SurrealFinite, y::SurrealFinite)
    xr = convert(Rational, x)
    yr = convert(Rational, y) 
    if y ≅ zero(y)
        error(InexactError(:/, "x/y", "inexact division") )
    elseif y ≅ 2
        return x * convert(SurrealFinite,  1 // 2)
    elseif isinteger(y) && ispow2(yr.num)
        return x * convert(SurrealFinite,  1 // yr) 
    elseif isinteger(xr.num/yr)
        return convert(SurrealFinite,  (xr.num/yr) // xr.den  ) 
    else
        error(InexactError(:/, "x/y", "inexact division"))        
    end
end


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
    surreal2dag(x::SurrealFinite; direction::String="forward")
    surreal2dag(io::IO, x::SurrealFinite; direction::String="forward")

 Writes a surreal representation as a DAG out in DOT format for drawing using GraphVis,
 and returns the number of nodes in the graph. Returns the number of nodes in the graph. 

## Arguments
* `io::IO`: output stream, default is stdout
* `x::SurrealFinite`: the number to write out
* `direction::String="forward"`: by default arrows point towards parents ("forward"), to go other way  use "back"

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
function surreal2dag(io::IO, x::SurrealFinite; direction::String="forward")
    if !(direction=="forward" || direction=="back")
        error("direction should be \"forward\" or \"back\".")
    end
    println(io, "digraph \"", float(x), "\" {")
    k = 1
    SurrealsinDAG = Dict{SurrealFinite,Int}()
    m = surreal2dag_f(io, x, k, SurrealsinDAG; direction=direction)
    println(io, "}")
    return m
end
function surreal2dag_f(io::IO, x::SurrealFinite, k::Integer, SurrealsinDAG; direction::String="forward")
    m = k 
    if !haskey(SurrealsinDAG, x)
        SurrealsinDAG[x] = m 
        surreal2node(io, x, k) 
        c = 1
        for s in x.L
            if !haskey(SurrealsinDAG, s)
                m += 1
                # println(io, "   node_$k:L -> node_$m;")
                println(io, "   node_$k:\"" *  convert(String, s) * "," * string(c) * "\" -> node_$m [color=\"red3\", dir=$direction];")
                m = surreal2dag_f(io, s, m, SurrealsinDAG; direction=direction)
            else
                println(io, "   node_$k:\"" *  convert(String, s) * "," * string(c)  * "\" -> node_$(SurrealsinDAG[s]) [color=\"red3\", dir=$direction];") 
            end
            c += 1 
        end
        c = 1
        for s in x.R
            if !haskey(SurrealsinDAG, s)
                m += 1
                # println(io, "   node_$k:R -> node_$m;")
                println(io, "   node_$k:\"" *  convert(String, s) * "," * string(c) * "\" -> node_$m [color=\"blue3\", dir=$direction];")
                m = surreal2dag_f(io, s, m, SurrealsinDAG; direction=direction)
            else
                println(io, "   node_$k:\"" *  convert(String, s) * "," * string(c) * "\" -> node_$(SurrealsinDAG[s]) [color=\"blue3\", dir=$direction];")
            end
            c += 1 
        end
    end 
    return m 
end
surreal2dag(x::SurrealFinite; direction::String="forward") = surreal2dag(stdout, x; direction=direction)


""" 
    surreal2dot(x::SurrealFinite)
    surreal2dot(io::IO, x::SurrealFinite)

 Writes a surreal representation as a tree out in DOT format for drawing using GraphVis,
 and returns the number of nodes in the graph.

## Arguments
* `io::IO`: output stream, default is stdout
* `x::SurrealFinite`: the number to write out
* `direction::String="forward"`: by default arrows point towards parents ("forward"), to go other way  use "back"

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
function surreal2dot(io::IO, x::SurrealFinite; direction::String="forward")
    if !(direction=="forward" || direction=="back")
        error("direction should be \"forward\" or \"back\".")
    end
    println(io, "digraph \"", float(x), "\" {")
    k = 1
    m = surreal2dot_f(io, x, k; direction=direction)
    println(io, "}")
    return m
end
function surreal2dot_f(io::IO, x::SurrealFinite, k::Integer; direction::String="forward")
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
            println(io, "   node_$k:\"" *  convert(String, s) * "," * string(c) * "\" -> node_$m [color=\"red3\", dir=$direction];")
            m = surreal2dot_f(io, s, m; direction=direction)
            c += 1 
        end
        c = 1
        for s in x.R
            m += 1
            # println(io, "   node_$k:R -> node_$m;")
            println(io, "   node_$k:\"" *  convert(String, s) * "," * string(c)  * "\" -> node_$m [color=\"blue3\", dir=$direction];")
            m = surreal2dot_f(io, s, m; direction=direction)
           c += 1 
         end
    end 
    return m 
end
function surreal2node(io::IO, x::SurrealFinite, k::Integer; extra_args::String="")
    S = convert(String, x)
    S2 = overbarString(x)
    if isempty(x.L)
        # L = "<FONT POINT-SIZE=\"20\">∅</FONT>" 
        L = "Ø" 
    else
        L = "<TABLE BORDER=\"0\" CELLBORDER=\"0\" CELLPADDING=\"0\"><TR>"
        c = 1
        for s in x.L
            tmp = convert(String, s)
            tmp2 = overbarString(s)
            L *= "<TD CELLPADDING=\"1pt\" PORT=\"$tmp," * string(c) * "\"> " * tmp2  * " </TD> &nbsp; "
            c += 1
        end
        L *= "</TR></TABLE>" 
    end 
    if isempty(x.R)
        R = "Ø" 
        # R = "∅" 
        # R = "&#8709;" 
        # R = "&#x2205;" 
    else
        R = "<TABLE BORDER=\"0\" CELLBORDER=\"0\" CELLPADDING=\"0\"><TR>"
        c = 1
        for s in x.R
            tmp = convert(String, s)
            tmp2 = overbarString(s)
            R *= "<TD CELLPADDING=\"1pt\"  PORT=\"$tmp," * string(c) * "\"> " * tmp2 * " </TD> &nbsp; "
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
                         <TR><TD  CELLPADDING=\"5pt\" COLSPAN=\"2\">$S2</TD></TR>
                         <TR><TD PORT=\"L\"> $L </TD><TD PORT=\"R\"> $R </TD></TR>
                         </TABLE>>,$extra_args
                         ];""")
end
surreal2dot(x::SurrealFinite; direction::String="forward") = surreal2dot(stdout, x; direction=direction)

#######################################################

""" 
    generation(x::SurrealFinite)

 Finds the birthday of a surreal number, which is 1 + the max of any of its components.

## Arguments
* `x::SurrealFinite`: the number to operate on
     
##
```jldoctest
julia> generation( convert(SurrealFinite, 1) )
1
``` 
"""
generation(s::SurrealFinite) = dag_stats(s).generation
# function generation(x::SurrealFinite)
    #
    # old version, based on definition could be very slow as doesn't cache
    # if iszero(x)
    #     return 0
    # else
    #     return max( maximum( generation.( [x.L; 0]) ),
    #                 maximum( generation.( [x.R; 0]) )) + 1
    # end
# end 

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
function canonicalise(s::SurrealFinite) 
    # should have a better, more direct approach
    if isinteger(s)
        convert(SurrealFinite, convert(Int64, s))
        # following should work, but don't want to intefere with floor, which should use this
        # a = SurrealFinite("floor($(s.shorthand))", [ floor(s.maxL) ],  ∅)
    else
        convert(SurrealFinite, convert(Rational, s))
    end
end
function iscanonical_old(s::SurrealFinite) 
    global ExistingCanonicalsC
    hs = hash(s)
    if !haskey(ExistingCanonicalsC, hs)
        ExistingCanonicalsC[hs] = canonicalise(s) == s
    end
    return ExistingCanonicalsC[hs]
end
function iscanonical(s::SurrealFinite)
    global ExistingCanonicalsC
    global SurrealZero
    hs = hash(s)
    if !haskey(ExistingCanonicalsC, hs)
        if iscanonicalinteger(s)
            ExistingCanonicalsC[hs] = true
        elseif isinteger(s)
            ExistingCanonicalsC[hs] = false
        elseif length(s.L)>1 || length(s.R)>1
            ExistingCanonicalsC[hs] = false
        elseif length(s.L)==0 || length(s.R)==0
            ExistingCanonicalsC[hs] = false
        elseif s.L[1] < SurrealZero && s.R[1] > SurrealZero
            ExistingCanonicalsC[hs] = false
        elseif iscanonical(s.L[1]) && iscanonical(s.R[1]) 
            # s.L[1] -- s -- s.R[1] must occur in a sequence with same difference between them
            # m1 = s - s.L[1]
            # m2 = s.R[1] - s
            # or we could do s+s = s.L[1] + s.R[1]
            if s + s ≅ s.L[1] + s.R[1]
                ExistingCanonicalsC[hs] = true
            else
                ExistingCanonicalsC[hs] = false
            end
        else
            ExistingCanonicalsC[hs] = false
        end
    end
    return ExistingCanonicalsC[hs]
end

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
# N.B. returns canonical form of floor
function floor_0(T::Type, s::SurrealFinite)
    global ExistingFloors
    global ExistingSurreals
    global SurrealZero
    global SurrealOne
    global SurrealTwo
    global SurrealMinusOne
    global SurrealMinusTwo
    if haskey(ExistingFloors, hash(s))
        return convert(T, ExistingSurreals[ ExistingFloors[ hash(s) ] ])
    elseif s < SurrealZero
        if s >= SurrealMinusOne
            ExistingFloors[ hash(s) ] = hash(SurrealMinusOne)
            return convert(T, SurrealMinusOne)
        elseif s >= SurrealMinusTwo
            ExistingFloors[ hash(s) ] = hash(SurrealMinusTwo)
            return convert(T, SurrealMinusTwo)
        else
            k = 1
            while s < -2^(k+1) && k < 12 # 12 is arbitrary, but it would be painful to go higher
                k += 1
            end 
            if k==12
                error("s is too large (and negative) for the current floor function")
            end
            a = -2^(k+1)
            b = -2^k
            for i=1:k
                d = (a+b) / 2  # this is real number arithmetic to avoid dividing by 2 in surreals
                if s < d       # implicit conversions to surreals to keep everything canonical
                    b = d
                else
                    a = d
                end
            end
            a = dali(a)
            ExistingFloors[ hash(s) ] = hash(a)
            return convert(T, a) 
        end
    elseif s < SurrealOne
        ExistingFloors[ hash(s) ] = hash(SurrealZero)
        return zero(T)
    elseif s < SurrealTwo
        ExistingFloors[ hash(s) ] = hash(SurrealOne)
        return one(T)
    else
        # start with geometric search to bound the number
        k = 1
        while s >= 2^(k+1) && k < 12 # 12 is arbitrary, but it would be super painful to go higher
            k += 1
        end  
        if k==12
            error("s is too large for the current floor function")
        end
        # now a binary search to narrow it down
        a = 2^k
        b = 2^(k+1)
        for i=1:k
            d = (a+b) / 2  # this is real number arithmetic to avoid dividing by 2 in surreals
            if s < d       # implicit conversions to surreals to keep everything canonical
                b = d
            else
                a = d
            end
        end
        a = dali(a) # a wasn't implicitly converted yet
        ExistingFloors[ hash(s) ] = hash(a)
        return convert(T, a)
    end
end
function floor_1(T::Type, s::SurrealFinite)
    global ExistingFloors
    global ExistingSurreals
    global SurrealZero
    global SurrealOne
    global SurrealTwo
    global SurrealMinusOne
    global SurrealMinusTwo
    if haskey(ExistingFloors, hash(s))
        return convert(T, ExistingSurreals[ ExistingFloors[ hash(s) ] ])
    elseif s >= SurrealZero
        if  s < SurrealOne
            ExistingFloors[ hash(s) ] = hash(SurrealZero)
            return zero(T)
        elseif s < SurrealTwo
            ExistingFloors[ hash(s) ] = hash(SurrealOne)
            return one(T)
        else # could do a special case for canonical integers maybe? 
            if isinteger(s) # can include this in the new version only -- old version of isinteger used floor
                #   println("   floor of an integer $s")
                # a = canonicalise(s), can't use canonicalise because it converts, which uses floor
                a = SurrealFinite("floor($(s.shorthand))", [ floor(s.maxL) ],  ∅)
            else
                # only the max of the left integers matters, and we take the floor of that
                #   println("   floor of a non integer $s")
                a = floor( s.maxL )
            end
            ExistingFloors[ hash(s) ] = hash(a)
            return convert(T, a)
        end
    elseif s >= SurrealMinusOne # negative s is implicit
        ExistingFloors[ hash(s) ] = hash(SurrealMinusOne)
        return convert(T, SurrealMinusOne)
    elseif s >= SurrealMinusTwo # negative s is implicit
        ExistingFloors[ hash(s) ] = hash(SurrealMinusTwo)
        return convert(T, SurrealMinusTwo)
    else # negative s is implicit
        a = - ceil( -s )
        ExistingFloors[ hash(s) ] = hash(a)
        return convert(T, a)
    end
end

floor(T::Type, s::SurrealFinite) = floor_1(T, s)
floor(s::SurrealFinite) = floor(SurrealFinite, s)
function ceil(T::Type, s::SurrealFinite)
    if isinteger(s)
        return floor(T,s)
    else
        return floor(T,s+1) 
    end
end
ceil(s::SurrealFinite) = ceil(SurrealFinite, s)

trunc(T::Type, s::SurrealFinite) = s>=0 ? floor(T,s) : -floor(T,-s)
trunc(s::SurrealFinite) = trunc(SurrealFinite, s)

# implicit rounding mode is 'RoundNearestTiesUp'
#   to be consistent, should do a precision, but is hard, and probably should just throw an error
round(s::SurrealFinite, r::RoundingMode{:ToZero})  = trunc(s)
round(s::SurrealFinite, r::RoundingMode{:Down})    = floor(s)
round(s::SurrealFinite, r::RoundingMode{:Up})      = ceil(s)
round(s::SurrealFinite, r::RoundingMode{:Nearest}) = floor(s + 1//2)

# fallbacks
# floor(::Type{T}, x::AbstractFloat) where {T<:Integer} = trunc(T,round(x, RoundDown))
# ceil(::Type{T}, x::AbstractFloat) where {T<:Integer} = trunc(T,round(x, RoundUp))
# round(::Type{T}, x::AbstractFloat) where {T<:Integer} = trunc(T,round(x, RoundNearest))


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

function iscanonicalinteger(s::SurrealFinite)
    global ExistingCanonicalInteger
    global ExistingInteger
    global SurrealZero
    
    hs = hash( s )
    if haskey(ExistingCanonicalIntegers, hs)
        return ExistingCanonicalIntegers[hs]
    elseif iszero(s)
        ExistingCanonicalIntegers[ hs ] = true
        ExistingIntegers[ hs ] = true
    elseif isempty(s.L) && length(s.R)==1 && iscanonicalinteger(s.R[1]) && s.R[1] <= SurrealZero
        ExistingCanonicalIntegers[ hs ] = true
        ExistingIntegers[ hs ] = true
    elseif isempty(s.R) && length(s.L)==1 && iscanonicalinteger(s.L[1]) && s.L[1] >= SurrealZero
        ExistingCanonicalIntegers[ hs ] = true
        ExistingIntegers[ hs ] = true
    else
        ExistingCanonicalIntegers[ hs ] = false
    end
    return ExistingCanonicalIntegers[ hs ]
end
function isinteger_old(s::SurrealFinite)
    global ExistingInteger

    hs = hash( s )
    if haskey(ExistingIntegers, hs)
        return ExistingIntegers[hs]
    elseif iscanonicalinteger(s)
        return true
    elseif s ≅ floor(s)
        ExistingIntegers[ hs ] = true
        return true
    else
        ExistingIntegers[ hs ] = false
        return false
    end
end
function isinteger(s::SurrealFinite)
    global ExistingInteger
    global SurrealOne

    hs = hash( s )
    if haskey(ExistingIntegers, hs)
        return ExistingIntegers[hs]
    elseif iscanonicalinteger(s)
        result = true
    elseif length(s.L)==0 || length(s.R)==0
        result = true
    elseif s.maxL + SurrealOne < s.minR # if the gap is more than 1, then an integer will always fit
        result = true
    elseif isinteger(s.maxL) || isinteger(s.minR) 
        # gap is <= 1, so they can't be either side of an integer in this case
        result = false
    elseif floor(s.maxL) < floor(s.minR)
        # the left and right parts are on either side of an integer, and we can't use floor on s but on its components is OK?
        result = true  
    else
        result = false
    end
    ExistingIntegers[ hs ] = result
    return result
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
abstract type SurrealStats end

# structure to hold stats on a surreal form
struct SurrealDAGstats <: SurrealStats
    nodes::Int64 # nodes in its DAG
    tree_nodes::Int64 # nodes in a tree-based representation of the graph
    edges::Int64 # edges in its DAG
    generation::Int32 # generation/birthday of the surreal, i.e., shortest path to root at zero
    longest_path::Array{SurrealFinite,1} 
    paths::BigInt # number of paths from source to sink
    value::Rational
    minval::Rational 
    maxval::Rational
    n_zeros::Int64 # number of nodes with value zero: only calculate if V switch is true
    max_width::Int64
    gap::Float64 # diff between MaxL and minR as a value, noting that it could be Inf, so use Floats
end
const ExistingSurrealDAGstats = Dict{HashType, SurrealDAGstats}()

mutable struct SurrealDegreeDist <: SurrealStats
    s::SurrealFinite
    nodes::Int64
    in_degree::Dict{UInt64, Int64}
    out_degree::Dict{UInt64, Int64}
    left_in_degree::Dict{UInt64, Int64}
    left_out_degree::Dict{UInt64, Int64}
    right_in_degree::Dict{UInt64, Int64}
    right_out_degree::Dict{UInt64, Int64}
end
SurrealDegreeDist(s::SurrealFinite) = SurrealDegreeDist( s, 0, Dict{UInt64, Int64}(), Dict{UInt64, Int64}(), Dict{UInt64, Int64}(), Dict{UInt64, Int64}(), Dict{UInt64, Int64}(), Dict{UInt64, Int64}() )

function ==(a::T, b::T) where T <: SurrealStats
    for name in fieldnames(T)
        if getfield(a,name) != getfield(b,name)
            return false
        end
    end 
    return true
end 

#   keep local cache for this calculation because some stats of each form in this DAG are different depending on context
function degree_stats( s::SurrealFinite, stats::SurrealDegreeDist, processed_list::Dict{SurrealFinite,Bool}; ) 
    h = hash(s)
    P = parents(s)
    stats.nodes += 1
    stats.in_degree[h] = length(P) # in degree is the number of parents of a node
    stats.left_in_degree[h] = length(s.L) # in degree is the number of parents of a node
    stats.right_in_degree[h] = length(s.R) # in degree is the number of parents of a node
    if !haskey(stats.out_degree, h)
        stats.out_degree[h] = 0 # its not currently a parent of anything
        stats.left_out_degree[h] = 0 # its not currently a parent of anything
        stats.right_out_degree[h] = 0 # its not currently a parent of anything
    end
    for p in P
        if !haskey(stats.out_degree, hash(p))
            stats.out_degree[hash(p)] = 1
            if p in s.L
                stats.left_out_degree[hash(p)] = 1
                stats.right_out_degree[hash(p)] = 0
            elseif p in s.R
                stats.left_out_degree[hash(p)] = 0
                stats.right_out_degree[hash(p)] = 1
            else
                error("  parent is neighter left or right")
            end
        else
            stats.out_degree[hash(p)] += 1
            if p in s.L
                stats.left_out_degree[hash(p)] += 1
            elseif p in s.R
                stats.right_out_degree[hash(p)] += 1
            else
                error("  parent is neighter left or right")
            end
        end
        if !haskey(processed_list, p)
            (stats, processed_list) = degree_stats( p, stats, processed_list )
        end  
    end
    processed_list[s] = true
    return ( stats, processed_list )  
end
function print_degree_stats( stats )
    println("degree stats for $(stats.s) with $(stats.nodes) nodes")
    for x in keys(stats.in_degree)
        println("   \"$(ExistingSurreals[x])\", $(stats.in_degree[x]), $(stats.out_degree[x]), $(stats.left_in_degree[x]), $(stats.left_out_degree[x]), $(stats.right_in_degree[x]), $(stats.right_out_degree[x])") 
    end
#     println("  in degrees")
#    for x in keys(stats.in_degree)
#        println("   $(ExistingSurreals[x]): $(stats.in_degree[x])") 
#   end
#    println("  out degrees")
#    for x in keys(stats.out_degree)
#        println("   $(ExistingSurreals[x]): $(stats.out_degree[x])")         
#    end
end
function degree_dist( stats )
    in_degrees = [ stats.in_degree[x] for x in keys(stats.in_degree)]
    out_degrees = [ stats.out_degree[x] for x in keys(stats.in_degree)]
    left_in_degrees = [ stats.left_in_degree[x] for x in keys(stats.in_degree)]
    left_out_degrees = [ stats.left_out_degree[x] for x in keys(stats.in_degree)]
    right_in_degrees = [ stats.right_in_degree[x] for x in keys(stats.in_degree)]
    right_out_degrees = [ stats.right_out_degree[x] for x in keys(stats.in_degree)]
    
    max_degree = max( maximum(in_degrees), maximum(out_degrees) )
    in_deg_dist = zeros(Int64, max_degree+1)
    out_deg_dist = zeros(Int64, max_degree+1)
    left_in_deg_dist = zeros(Int64, max_degree+1)
    left_out_deg_dist = zeros(Int64, max_degree+1)
    right_in_deg_dist = zeros(Int64, max_degree+1)
    right_out_deg_dist = zeros(Int64, max_degree+1)
    
    for n in in_degrees
        in_deg_dist[n+1] += 1
    end
    for n in out_degrees
        out_deg_dist[n+1] += 1
    end 
    for n in left_in_degrees
        left_in_deg_dist[n+1] += 1
    end
    for n in left_out_degrees
        left_out_deg_dist[n+1] += 1
    end 
    for n in right_in_degrees
        right_in_deg_dist[n+1] += 1
    end
    for n in right_out_degrees
        right_out_deg_dist[n+1] += 1
    end 
    return in_deg_dist, out_deg_dist, left_in_deg_dist, left_out_deg_dist, right_in_deg_dist, right_out_deg_dist
end

#    LP=true means keep track of a longest path
#    V =true means keep track of values
# have these as switches, because they are expensive to calculate on large DAGs
function dag_stats( s::SurrealFinite, processed_list::Dict{SurrealFinite,SurrealDAGstats}; LP=false, V=true) 
    if s == zero(s)    
        nodes = 1  
        tree_nodes = 1
        edges = 0
        generation = 0
        longest_path = LP ? [zero(s)] : ϕ # important this be \phi and not an arbitrary empty array
        paths = 1
        value = 0 
        gap = 0.0
        minval = 0  
        maxval = 0 
        n_zeros = V ? 1 : 0
        max_width = 0
    else
        P = parents(s) 
        nodes = 1
        tree_nodes = 1
        edges = length(P)
        generation = 0
        longest_path = LP ? [zero(s)] : ϕ
        paths = 0 
        value = V ? convert(Rational, s) : 0
        if V
            if ismissing(s.minR)
                gap = Inf
            elseif ismissing(s.maxL)
                gap = -Inf
            else
                gap = convert(Rational, s.minR) - convert(Rational, s.maxL)
            end
        else    
            gap = NaN
        end
        minval = value 
        maxval = value
        n_zeros = V ? Int64(value == 0) : 0
        max_width = length(P)
        for p in P 
             if !haskey(processed_list, p)
                (stats,l) = dag_stats(p, processed_list; LP=LP, V=V) 
                nodes += stats.nodes 
                edges += stats.edges 
                minval = min(minval, stats.minval)
                maxval = max(maxval, stats.maxval)
                n_zeros += stats.n_zeros
                max_width = max(max_width, stats.max_width)
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
    stats = SurrealDAGstats(nodes, tree_nodes, edges, generation, longest_path, paths, value, minval, maxval, n_zeros, max_width, gap) 
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
function increment!( d::Dict{S, T}, k::S, i::T) where {T<:Real, S<:Any}
    if haskey(d, k)
        d[k] += i
    else
        d[k] = i
    end
end
increment!(d::Dict{S, T}, k::S ) where {T<:Real, S<:Any} = increment!( d, k, one(T))

function decrement!( d::Dict{S, T}, k::S, i::T) where {T<:Real, S<:Any}
    if haskey(d, k)
        d[k] -= i
    else
        d[k] = -i
    end
end
decrement!(d::Dict{S, T}, k::S ) where {T<:Real, S<:Any} = increment!( d, k, one(T))

function uniqueness(s::SurrealFinite,
                    P::Dict{HashType, HashType},
                    Q::Dict{HashType, SurrealFinite},
                    U::Dict{HashType, Int64},
                    V::Dict{HashType, SurrealFinite} )
    h = hash(s)
    increment!(U, h)
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
                                          Dict{HashType, HashType}(),
                                          Dict{HashType, SurrealFinite}(),
                                          Dict{HashType, Int64}(),
                                          Dict{HashType, SurrealFinite}() )
uniqueness_max( U::Dict{HashType, Int64} ) = maximum( values(U) )
function uniqueness_failure(  U::Dict{HashType, Int64}, V::Dict{HashType, SurrealFinite} )
    
    @static if VERSION < v"1.0.0"
        U2 = filter( (k,v) -> v>1 , U) # not Julia 1.0 compliant
    else
        U2 = filter( ((k,v),) -> v>1 , U) 
    end 
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


