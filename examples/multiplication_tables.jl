using IndexedTables
using FileIO
using LatexPrint
using SurrealNumbers 
out_dir = "Data/"
  
convert(SurrealFinite, 1) * convert(SurrealFinite, 2) 
a = convert(SurrealFinite, 2) * convert(SurrealFinite, 2) 
b = convert(SurrealFinite, 1) - convert(SurrealFinite, 1) 

n2 = 12
n3 = 4
n4 = 2
n_12 = 12
n_14 = 3

# n2 = 8 
# n3 = 2
# n4 = 2
# n_12 = 4
# n_14 = 2

ns1 = convert.(SurrealFinite, [0; 1:n2;        2;      1:n3;       1:n4;          1:n_12;          1:n_14; 1//2; 1//2; 1//4; a; b])
ns2 = convert.(SurrealFinite, [2; 2*ones(n2); -2; 3*ones(n3); 4*ones(n4); 1//2*ones(n_12); 1//4*ones(n_14); 1//2; 1//4; 1//4; 2; 2]) 
# ns1 = convert.(SurrealFinite, [0:n2-1;           1:n3;       1:n4;         1:n_12;           1:n_14; 1//2; 1//2])
# ns2 = convert.(SurrealFinite, [2*ones(n2); 3*ones(n3); 4*ones(n4); 1//2*ones(n_12); 1//4*ones(n_14); 1//2; 1//4]) 
# ns1 = convert.(SurrealFinite, [0; 1; 2; 3; 4; 5; 6; 7; 8; 9; 3; a]) 
# ns2 = convert.(SurrealFinite, [2; 2; 2; 2; 2; 2; 2; 2; 2; 2; 3; 2])
# ns1 = [2; 4]
# ns2 = [1/2; 1/4] 
n = length(ns1) 
# m = [10; 10; 10; 10; 10; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1]
m = ones(Integer, n)
x = Array{SurrealFinite,1}(n) 
t = Array{Float64,1}(n) 
bytes = Array{Float64,1}(n)  
gctime = Array{Float64,1}(n)
c_plus = Array{Integer,1}(n) 
c_times = Array{Integer,1}(n) 
c_created = Array{Integer,1}(n) 
c_equals = Array{Integer,1}(n) 
c_leq = Array{Integer,1}(n) 
i = 0
println("calculating products")
for i=1:n 
    println("   ", ns1[i], " x ", ns2[i])
    t[i] = 0 
    b = 0 
    g = 0
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
    c_leq[i] = Count['≦'] 
end
println("calculating stats")
# s = dag_stats.( x )
s = Array{SurrealDAGstats}(n)
for i=1:n
    println("   ", ns1[i], " x ", ns2[i])
    s[i] = dag_stats( x[i]; V=false ) 
end 
no = [s[i].nodes for i=1:length(s)]
ed = [s[i].edges for i=1:length(s)]
p = [s[i].paths for i=1:length(s)]
d = (ed+1) ./ no 
g = [s[i].generation for i=1:length(s)]
b = [s[i].maxval - s[i].minval for i=1:length(s)]
t = t./m[1:n]
bytes = bytes

# add in 2x2x2

t1 = table(@NT(i=0:n-1,
               n1=Rational.(ns1), 
               n2=Rational.(ns2), 
               generation=g,
               nodes=no,
               edges=ed,
               paths=float.(p),
               density=d,
               time=t,
               bytes=bytes,
               adds=c_plus,
               prods=c_times,
               created=c_created,
               equals=c_equals,
               leq=c_leq);
           pkey = [:i])

save("$(out_dir)multiplication_tables.csv", t1)



A = Array{Any}(n+1,9);
A[1,:] = ["\$x \\times y\$" "g(xy)" "s(xy)" "e(xy)" "p(xy)" "time (s)"  "+" "*" "created"]
for j=2:n+1
    if isinteger(ns1[j-1])
        tmp1 = Int(ns1[j-1])
    else
        tmp1 = Rational(ns1[j-1])
    end
    if isinteger(ns2[j-1])
        tmp2 = Int(ns2[j-1])
    else
        tmp2 = Rational(ns2[j-1])
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
        A[j,6] = round(t[j-1],2)
    elseif t[j-1]<1000
        A[j,6] = round(t[j-1],1)
    else
        A[j,6] = round(t[j-1],0)
    end
end
A[2:n+1,7] = c_plus
A[2:n+1,8] = c_times
A[2:n+1,9] = c_created

# A[2:n+1,5] = Int.(x)
# A[:,5] = 3^Int.(ns1.+ns2)
set_align("r")
# tabular(STDOUT, A,"rr|rr")

file = "$(out_dir)multiplication_table.tex"
FID = open(file, "w")
L = latex_form(A) 
println(FID, L)
close(FID)


surreal2tex(x[3], level=2)

for i=1:4
    file = "$(out_dir)multiplication_ex_$i.dot"
    println("outputting $file")
    FID = open(file, "w")
    surreal2dag(FID, x[i])
    close(FID)
    run(`dot -Tpdf -O $file`)
end


# look at this in more detail 
y = convert(SurrealFinite, 1) + convert(SurrealFinite, 2)*convert(SurrealFinite,2)
float.(y.L)

# also should do   2^k * 2^(-k), and see how far we get

