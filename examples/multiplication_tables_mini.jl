# using IndexedTables
using DataFrames
@static if VERSION < v"0.7.0"
else
    using Printf # @sprintf moved here from 0.6 -> 1.0
end
using FileIO
using LatexPrint
using SurrealNumbers
# using CSVFiles # earlier Julia version seem to need this explicitly
using CSV
out_dir = "Data"
fig_dir = "Figs"
include("write_machine_details.jl")

convert.(SurrealFinite, 1:12)
convert(SurrealFinite, 1) * convert(SurrealFinite, 2)
a = convert(SurrealFinite, 2) * convert(SurrealFinite, 2)
b = convert(SurrealFinite, 1) - convert(SurrealFinite, 1)

n2 = 12
n3 = 5
n4 = 2
n_12 = 12
n_14 = 3

# # short time version
# n2 = 10
# n3 = 3
# n4 = 2
# n_12 = 8
# n_14 = 3

clearcache()
ns1 = convert.(SurrealFinite, [0; 1:n2;        2;      1:n3;       1:n4;          1:n_12;          1:n_14; 1//2; 1//2; 1//4; a; b])
ns2 = convert.(SurrealFinite, [2; 2*ones(n2); -2; 3*ones(n3); 4*ones(n4); 1//2*ones(n_12); 1//4*ones(n_14); 1//2; 1//4; 1//4; 2; 2])
# ns1 = convert.(SurrealFinite, [0:n2-1;           1:n3;       1:n4;         1:n_12;           1:n_14; 1//2; 1//2])
# ns2 = convert.(SurrealFinite, [2*ones(n2); 3*ones(n3); 4*ones(n4); 1//2*ones(n_12); 1//4*ones(n_14); 1//2; 1//4])
# ns1 = convert.(SurrealFinite, [0; 1; 2; 3; 4; 5; 6; 7; 8; 3; a])
# ns2 = convert.(SurrealFinite, [2; 2; 2; 2; 2; 2; 2; 2; 2; 3; 2])
# ns1 = [2; 4]
# ns2 = [1/2; 1/4]
n = length(ns1)

# m = [10; 10; 10; 10; 10; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1]
m = ones(Integer, n)

ns1 = convert.(SurrealFinite, [3])
ns2 = convert.(SurrealFinite, [4])
# ns1 = convert.(SurrealFinite, [12])
# ns2 = convert.(SurrealFinite, [1//2])
ns1 = convert.(SurrealFinite, [12])
ns2 = convert.(SurrealFinite, [2])
n = length(ns1)

x = Array{SurrealFinite,1}(undef,n)
t = Array{Float64,1}(undef,n)
bytes = Array{Float64,1}(undef,n)
gctime = Array{Float64,1}(undef,n)
c_plus = Array{Integer,1}(undef,n)
c_times = Array{Integer,1}(undef,n)
c_created = Array{Integer,1}(undef,n)
c_equals = Array{Integer,1}(undef,n)
c_leq = Array{Integer,1}(undef,n)

i = 0
println("calculating products")
for i=1:n
    println("   ", ns1[i], " x ", ns2[i])
    t[i] = 0
    local b = 0
    local g = 0
    mem = 0
    val = 0
    for j=1:m[i]
        clearcache()
        val, time, b, g, mem = @timed ns1[i] * ns2[i];
        t[i] += time
    end
    x[i] = val
    bytes[i] = b
    gctime[i] = g

    c_plus[i] = Count['+']
    c_times[i] = Count['*']
    c_created[i] = Count['c']
    c_equals[i] = Count['=']
    c_leq[i] = Count['≤']
end
println("calculating stats")
# s = dag_stats.( x )
@static if VERSION < v"0.7.0"
    s = Array{SurrealDAGstats}(n)
else
    s = Array{SurrealDAGstats}(undef, n)
end
for i=1:n
    println("   ", ns1[i], " x ", ns2[i])
    s[i] = dag_stats( x[i]; V=true )
end
no = [s[i].nodes for i=1:length(s)]
ed = [s[i].edges for i=1:length(s)]
# mw = [s[i].max_width for i=1:length(s)]
nz = [s[i].n_zeros for i=1:length(s)]
p = [s[i].paths for i=1:length(s)]
d = (ed.+1) ./ no
g = [s[i].generation for i=1:length(s)]
b = [s[i].maxval - s[i].minval for i=1:length(s)]
t = t./m[1:n]
maxv = [s[i].maxval for i=1:length(s)]
minv = [s[i].minval for i=1:length(s)]
version = [VERSION for i=1:length(s)]

t1 = DataFrame(version=version,
               n1= convert.(Rational, ns1),
               n2= convert.(Rational, ns2),
               generation=g,
               nodes=no,
               edges=ed,
               minval=minv,    # not calculating values
               maxval=maxv,    # not calculating values
               n_zeros=nz,     # not calculating values
               paths=float.(p),
               density=d,
               time=t,
               bytes=bytes,
               adds=c_plus,
               prods=c_times,
               created=c_created,
               equals=c_equals,
               leq=c_leq
               )

output_file = joinpath(@__DIR__, out_dir, "multiplication_table_mini_$(VERSION).csv")
CSV.write(output_file, t1)

machine_file = joinpath(@__DIR__, out_dir, "multiplication_tables__mini_$(VERSION).machine")
write_machine_details(machine_file)
