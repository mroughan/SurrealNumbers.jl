using IndexedTables
@static if VERSION < v"0.7.0"
else
    using Printf # @sprintf moved here from 0.6 -> 1.0
end
using FileIO
using LatexPrint
using SurrealNumbers
using CSVFiles # earlier Julia version seem to need this explicitly
out_dir = "Data"
fig_dir = "Figs"

convert.(SurrealFinite, 1:12)
convert(SurrealFinite, 1) * convert(SurrealFinite, 2)
a = convert(SurrealFinite, 2) * convert(SurrealFinite, 2)
b = convert(SurrealFinite, 1) - convert(SurrealFinite, 1)


clearcache()
# ns1 = convert.(SurrealFinite, [2;  3; 12;   4; 4    ])
# ns2 = convert.(SurrealFinite, [12; 5; 1//2; 4; 1//4 ])
# also want 2*2*2 and x-x 

# ns1 = convert.(SurrealFinite, [2;   12;   3  ])
# ns2 = convert.(SurrealFinite, [12;  1//2; 5  ])

examples = [ :(dali(2)*dali(12)); :(dali(12)*dali(1//2)); :(dali(2)*dali(2)*dali(2)); :(dali(200)-dali(200)) ]
n = length(examples)

# m = [10; 10; 10; 10; 10; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1]
m = ones(Integer, n)


x = Array{SurrealFinite,1}(undef,n)
t = Array{Float64,1}(undef,n)
bytes = Array{Float64,1}(undef,n)
gctime = Array{Float64,1}(undef,n)
c_plus = Array{Integer,1}(undef,n)
c_times = Array{Integer,1}(undef,n)
c_created = Array{Integer,1}(undef,n)
c_equals = Array{Integer,1}(undef,n)
c_leq = Array{Integer,1}(undef,n)
cu_plus = Array{Integer,1}(undef,n)
cu_times = Array{Integer,1}(undef,n)
cu_created = Array{Integer,1}(undef,n)
cu_equals = Array{Integer,1}(undef,n)
cu_leq = Array{Integer,1}(undef,n)
depth = Array{Integer,1}(undef,n)

i = 0
println("calculating products")
for i=1:n
    println("   $(examples[i])")
    t[i] = 0
    local b = 0
    local g = 0
    mem = 0
    val = 0
    for j=1:m[i]
        clearcache()
        val, time, b, g, mem = @timed eval(examples[i]);
        t[i] += time
    end
    x[i] = val
    bytes[i] = b
    gctime[i] = g

    c_plus[i] = Count['+']
    c_times[i] = Count['*']
    c_created[i] = Count['n']
    c_equals[i] = Count['=']
    c_leq[i] = Count['≤']
    cu_plus[i] = CountUncached['+']
    cu_times[i] = CountUncached['*']
    cu_created[i] = CountUncached['n'] # not used at the moment
    cu_equals[i] = CountUncached['=']
    cu_leq[i] = CountUncached['≤']

    depth[i] = RecursionDepth["max_depth"]
end
println("calculating stats")
# s = dag_stats.( x )
@static if VERSION < v"0.7.0"
    s = Array{SurrealDAGstats}(n)
else
    s = Array{SurrealDAGstats}(undef, n)
end
for i=1:n
    println("   $(examples[i])")
    s[i] = dag_stats( x[i]; V=false )
end
no = [s[i].nodes for i=1:length(s)]
ed = [s[i].edges for i=1:length(s)] 
mw = [s[i].max_width for i=1:length(s)]
nz = [s[i].n_zeros for i=1:length(s)]
p = [s[i].paths for i=1:length(s)]
d = (ed.+1) ./ no
g = [s[i].generation for i=1:length(s)]
b = [s[i].maxval - s[i].minval for i=1:length(s)]
t = t./m[1:n]
maxv = [s[i].maxval for i=1:length(s)]
minv = [s[i].minval for i=1:length(s)]


t1 = table((i=0:n-1,
                   calc=examples,
                   generation=g,
                   nodes=no,
                   edges=ed,
                   max_width=mw,
                   minval=minv,    # not calculating values
                   maxval=maxv,    # not calculating values
                   n_zeros=nz,     # not calculating values
                   paths=float.(p),
                   density=d,
                   time=t,
                   bytes=bytes,
                   depth=depth,
                   created=c_created,
                   adds=c_plus,
                   adds_new=cu_plus,
                   prods=c_times,
                   prods_new=cu_times,
                   equals=c_equals,
                   equals_new=cu_equals,
                   leq=c_leq,
                   leq_new=cu_leq);
               pkey = [:i])

save(joinpath(@__DIR__, out_dir, "challenging_examples_$(VERSION).csv"), t1)

A = Array{Any}(undef, n+1, 10);
A[1,:] = ["calculation" "generation" "nodes" "edges" "serialisation" "time (s)" "≤"  "+" "x" "depth"]
for j=2:n+1
    A[j,1] = "$(examples[j-1])"
end
A[2:n+1,2] = Int.(g)
A[2:n+1,3] = Int.(no)
A[2:n+1,4] = Int.(ed)
# A[2:n+1,5] = float.(p)
j = 1
for j=2:n+1
    if p[j-1] <= 0
        A[j,5] = " overflow "
    elseif p[j-1]<100000
        A[j,5] = p[j-1]
    else
        A[j,5] = @sprintf("%.1E",  float(p[j-1]))
    end
end
# for j=2:n+1
#     if t[j-1]<0.01
#         A[j,6] = "\$<\$0.01"
#     elseif t[j-1]<1000
#         A[j,6] = @sprintf("%.2f",  t[j-1])
#     else
#         A[j,6] = @sprintf("%.1f",  t[j-1])
#     end
# end
for j=2:n+1
    if t[j-1]<100
        # A[j,6] = round(t[j-1],2)
        A[j,6] = round(t[j-1]; digits=2)
    elseif t[j-1]<1000
        # A[j,6] = round(t[j-1],1)
        A[j,6] = round(t[j-1]; digits=1)
    else
        # A[j,6] = round(t[j-1],0)
        A[j,6] = round(t[j-1]; digits=0)
    end
end
A[2:n+1,7] = c_leq
A[2:n+1,8] = c_plus
A[2:n+1,9] = c_times
A[2:n+1,10] = depth


set_align("r")

file = joinpath(@__DIR__, out_dir, "challenging_examples_table_$(VERSION).tex")
FID = open(file, "w")
L = latex_form(A)
println(FID, L)
close(FID)


