# compare the time of two different implementations of floor
using IndexedTables
using SurrealNumbers

# old code for reference, e.g., in case you want to rerun some timing tests to understand
# why current code is written the way it is
function floor2(s::SurrealFinite)
    if zero(s) <= s < one(s) 
        return zero(s)
    elseif s < zero(s)
        return floor(s + one(s)) - one(s)
    elseif s >= one(s)
        return floor(s - one(s)) + one(s)
    end
end

function isinteger2(s::SurrealFinite)
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


# not the more general form of rounding defined in Julia -- should fix
function round2(s::SurrealFinite)
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
        return -one(s)
    end
end

# run a simple bit of code to ensure compilation
convert(SurrealFinite, 2) * convert(SurrealFinite, 1)

# setup the tests
ns = -4 : 0.125 : 4
n = length(ns)
m = 1000 * ones(size(ns))
x1 = Array{SurrealFinite,1}(undef, n) 
x2 = Array{SurrealFinite,1}(undef, n) 
t1 = Array{Float64,1}(undef, n) 
t2 = Array{Float64,1}(undef, n) 
bytes2 = Array{Float64,1}(undef, n) 
bytes1 = Array{Float64,1}(undef, n) 
gctime1 = Array{Float64,1}(undef, n) 
gctime2 = Array{Float64,1}(undef, n)
s = convert(SurrealFinite, 1)
floor(s)
floor2(s)
for i=1:n
    println("  working on ", ns[i])
    t1[i] = Inf
    t2[i] = Inf
    b1 = 0 
    g1 = 0 
    mem1 = 0
    val1 = 0
    b2 = 0 
    g2 = 0
    mem2 = 0
    val2 = 0
    for j=1:m[i]
        s = convert(SurrealFinite, ns[i] + 0.5)
        # val, time, b, g, mem = @timed convert(Rational, s) 
        clearcache()
        val1, time1, b1, g1, mem1 = @timed floor( s )
        clearcache()
        val2, time2, b2, g2, mem2 = @timed floor2( s )
        t1[i] = min(t1[i], time1)
        t2[i] = min(t2[i], time2)
    end
    x1[i] = val1
    bytes1[i] = b1
    gctime1[i] = g1
    x2[i] = val2
    bytes2[i] = b2
    gctime2[i] = g2
end
# t1 = t1./m[1:n]
bytes1 = bytes1./m[1:n]
# t2 = t2./m[1:n]
bytes2 = bytes2./m[1:n]

t1_p_ms = 1.0e6 * t1
t2_p_ms = 1.0e6 * t2
 
tab1 = table((i=ns,
                time1=t1_p_ms,
                time2=t2_p_ms,
                bytes1=bytes1,
                bytes2=bytes2);
           pkey = [:i])




