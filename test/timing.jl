using IndexedTables
using SurrealNumbers
convert(SurrealFinite, 2) * convert(SurrealFinite, 1)

ns = [0; 1; 2; 3; 4; 5]
n = length(ns)
m = [100; 100; 100; 30; 1; 1; 1]
x = Array{SurrealFinite,1}(n) 
t = Array{Float64,1}(n) 
bytes = Array{Float64,1}(n) 
gctime = Array{Float64,1}(n) 
for i=1:n
    println(" 2 x ", ns[i])
    t[i] = 0
    b = 0
    g = 0
    mem = 0
    val = 0
    for j=1:m[i]
        val, time, b, g, mem = @timed convert(SurrealFinite, 2) * convert(SurrealFinite, ns[i])
        t[i] += time
    end
    x[i] = val
    bytes[i] = b
    gctime[i] = g
end
s = size.(x)
d = depth_av.(x)
f = n_zeros.(x)
t = t./m[1:n]
bytes = bytes./m[1:n]

t_p_ms = 1.0e6 * t ./ s

t1 = table(@NT( i=0:n-1,
                size=s,
                depth=d,
                leaves=f,
                time=t,
                time_per_n=t_p_ms,
                bytes=bytes);
           pkey = [:i])




