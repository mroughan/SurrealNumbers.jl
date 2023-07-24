# estimate times for 2x table
using IndexedTables
using SurrealNumbers
convert(SurrealFinite, 2) * convert(SurrealFinite, 1)

ns = [0; 1; 2; 3; 4; 5; 6; 7; 8]
n = length(ns)
m = [100; 100; 100; 100; 100; 100; 30; 10; 3]
x = Array{SurrealFinite,1}(undef, n) 
t = Array{Float64,1}(undef, n) 
bytes = Array{Float64,1}(undef, n) 
gctime = Array{Float64,1}(undef, n) 
label = Array{String,1}(undef,n)
for i=1:n
    label[i] = "2 x $(ns[i])"
    println(label[i])
    t[i] = 0
    b = 0
    g = 0
    mem = 0
    val = 0
    for j=1:m[i]
        clearcache() 
        val, time, b, g, mem = @timed convert(SurrealFinite, 2) * convert(SurrealFinite, ns[i])
        t[i] += time
    end
    x[i] = val
    bytes[i] = b
    gctime[i] = g
end
stats = dag_stats.(x)
no = [stats[i].nodes for i=1:length(stats)]
ed = [stats[i].edges for i=1:length(stats)]
p = [stats[i].paths for i=1:length(stats)]
g = [stats[i].generation for i=1:length(stats)]
f = [stats[i].n_zeros for i=1:length(stats)]

t = t./m[1:n]
bytes = bytes./m[1:n]
t_p_ms = 1.0e6 * t ./ no

t1 = table(( i=0:n-1,
                label=label,
                nodes=no,
                edges=ed,
                depth=g,
                leaves=f,
                time_s=t,
                time_per_n_ms=t_p_ms,
                bytes=bytes);
           pkey = [:i])




