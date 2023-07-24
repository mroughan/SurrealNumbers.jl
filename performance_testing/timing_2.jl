# compare the time taken to do simple subtractions compared to <
# with a view to making some routines that use subtraction a lot faster
using IndexedTables
using SurrealNumbers
convert(SurrealFinite, 2) * convert(SurrealFinite, 1)

#ns = 6:-1:0
ns = -10 : 0.5 : 10
n = length(ns)
m = 1000 * ones(n)
x1 = Array{SurrealFinite,1}(undef, n) 
x2 = Array{Bool,1}(undef, n) 
t1 = Array{Float64,1}(undef, n) 
t2 = Array{Float64,1}(undef, n) 
tt1 = Array{Float64,1}(undef, n) 
tt2 = Array{Float64,1}(undef, n) 
bytes2 = Array{Float64,1}(undef, n) 
bytes1 = Array{Float64,1}(undef, n) 
gctime1 = Array{Float64,1}(undef, n) 
gctime2 = Array{Float64,1}(undef, n)
s = convert(SurrealFinite, 1)
s > 0
s - 1
for i=1:n
    println("  working on ", ns[i])
    t1[i] = 0
    t2[i] = 0
    tt1[i] = Inf
    tt2[i] = Inf
    b1 = 0 
    g1 = 0 
    mem1 = 0
    val1 = 0
    b2 = 0 
    g2 = 0
    mem2 = 0
    val2 = false
    for j=1:m[i]
        s = convert(SurrealFinite, ns[i])
        # val, time, b, g, mem = @timed convert(Rational, s) 
        clearcache()
        val1, time1, b1, g1, mem1 = @timed s - 1
        clearcache()
        val2, time2, b2, g2, mem2 = @timed s > 1
        t1[i] += time1
        t2[i] += time2
        tt1[i] = min( tt1[i], time1 )
        tt2[i] = min( tt2[i], time2 )
    end
    x1[i] = val1
    bytes1[i] = b1
    gctime1[i] = g1
    x2[i] = val2
    bytes2[i] = b2
    gctime2[i] = g2
end
t1 = t1./m[1:n]
bytes1 = bytes1./m[1:n]
t2 = t2./m[1:n]
bytes2 = bytes2./m[1:n]

t1_p_ms = 1.0e6 * t1
t2_p_ms = 1.0e6 * t2
tt1_p_ms = 1.0e6 * tt1
tt2_p_ms = 1.0e6 * tt2
 
tab1 = table( (i=ns,
                time_subtraction = tt1_p_ms,
                time_comparison  = tt2_p_ms,
                bytes_subtraction = bytes1,
                bytes_comparison  = bytes2);
           pkey = [:i])

           



