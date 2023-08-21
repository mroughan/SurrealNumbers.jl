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

# old bounds
n2 = 12
n3 = 5 # 3x6 is kiling me on memory -- up to 27G and still growing
n4 = 2 # 4x4 was killed by the box also, prolly too much memory
n_12 = 12
n_14 = 3

# # short time version
# n2 = 11
# n3 = 3
# n4 = 2
# n_12 = 11
# n_14 = 3

clearcache()
ns1 = convert.(SurrealFinite, [0; 1:n2;        2;      1:n3;       1:n4;          1:n_12;          1:n_14; 1//2; 1//2; 1//4; a; b])
ns2 = convert.(SurrealFinite, [2; 2*ones(n2); -2; 3*ones(n3); 4*ones(n4); 1//2*ones(n_12); 1//4*ones(n_14); 1//2; 1//4; 1//4; 2; 2])
n = length(ns1)

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
    cu_plus[i] = CountUncached['+']
    cu_times[i] = CountUncached['*']
    cu_created[i] = CountUncached['c'] # not used at the moment
    cu_equals[i] = CountUncached['=']
    cu_leq[i] = CountUncached['≤']
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
                   n1= convert.(Rational, ns1),
                   n2= convert.(Rational, ns2),
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

save(joinpath(@__DIR__, out_dir, "multiplication_table_$(VERSION).csv"), t1)

A = Array{Any}(undef, n+1,9);
A[1,:] = ["\$x \\times y\$" "g(xy)" "s(xy)" "e(xy)" "p(xy)" "time (s)"  "+" "*" "created"]
for j=2:n+1
    if isinteger(ns1[j-1])
        # tmp1 = Int(ns1[j-1])
        tmp1 = convert(Int, ns1[j-1])
    else
        # tmp1 = Rational(ns1[j-1])
        tmp1 = convert(Rational, ns1[j-1])
    end
    if isinteger(ns2[j-1])
        # tmp2 = Int(ns2[j-1])
        tmp2 = convert(Int, ns2[j-1])
    else
        # tmp2 = Rational(ns2[j-1])
        tmp2 = convert(Rational, ns2[j-1])
    end
    A[j,1] = "\$\\overline{$(tmp1)} \\times \\overline{$(tmp2)}\$"
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
A[2:n+1,7] = c_plus
A[2:n+1,8] = c_times
A[2:n+1,9] = c_created

# A[2:n+1,5] = Int.(x)
# A[:,5] = 3^Int.(ns1.+ns2)
set_align("r")
# tabular(STDOUT, A,"rr|rr")

file = joinpath(@__DIR__, out_dir, "multiplication_table_$(VERSION).tex")
FID = open(file, "w")
L = latex_form(A)
println(FID, L)
close(FID)

print("one example value = ")
surreal2tex(x[3], level=2)


for i=1:4
    file1 = joinpath(@__DIR__, fig_dir, "multiplication_ex_$i.dot")
    println("outputting $file1")
    FID1 = open(file1, "w")
    surreal2dag(FID1, x[i]; direction="back")
    close(FID1)
    run(`dot -Tpdf -O $file1`)
end


# look at this in more detail
y = convert(SurrealFinite, 1) + convert(SurrealFinite, 2)*convert(SurrealFinite,2)
float.(y.L)

# also should do   2^k * 2^(-k), and see how far we get

