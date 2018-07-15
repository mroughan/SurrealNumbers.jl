using IndexedTables
using FileIO
using LatexPrint
using SurrealNumbers 
out_dir = "Data/"
  
convert(SurrealFinite, 1) * convert(SurrealFinite, 2) 
a = convert(SurrealFinite, 2) * convert(SurrealFinite, 2) 

ns1 = convert.(SurrealFinite, [0; 1; 2; 3; 4; 5; 6; 7; 8; 9; 3; a])
ns2 = convert.(SurrealFinite, [2; 2; 2; 2; 2; 2; 2; 2; 2; 2; 3; 2])
# ns1 = [2; 4]
# ns2 = [1/2; 1/4]
n = length(ns1) 
m = [10; 10; 10; 10; 10; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1]
m = [ 1;  1;  1;  1;  1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1]
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
for i=1:n 
    println(ns1[i], " x ", ns2[i])
    t[i] = 0 
    b = 0
    g = 0
    mem = 0
    val = 0 
    for j=1:m[i]
        clearcache()
        val, time, b, g, mem = @timed ns1[i] * ns2[i]
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
s = dag_stats.( x )
no = [s[i].nodes for i=1:length(s)]
ed = [s[i].edges for i=1:length(s)]
p = [s[i].paths for i=1:length(s)]
d = (ed+1) ./ no 
g = [s[i].generation for i=1:length(s)]
t = t./m[1:n]
bytes = bytes

# add in 2x2x2

t1 = table(@NT(i=0:n-1,
               n1=Int.(ns1), 
               n2=Int.(ns2),
               nodes=no,
               edges=ed,
               paths=float.(p),
               density=d,
               generation=g,
               time=t,
               bytes=bytes,
               adds=c_plus,
               prods=c_times,
               created=c_created,
               equals=c_equals,
               leq=c_leq);
           pkey = [:i])

save("$(out_dir)multiplication_tables.csv", t1)



A = Array(Any,(n+1,9));
A[1,:] = ["x \times y" "g(xy)" "s(xy)" "e(xy)" "p(xy)" "time (s)"  "+" "*" "created"]
for j=2:n+1 
    A[j,1] = "\$\\overline{$(Int(ns1[j-1]))} \\times \\overline{$(Int(ns2[j-1]))}\$"
end
A[2:n+1,2] = Int.(g)
A[2:n+1,3] = Int.(no)
A[2:n+1,4] = Int.(ed)
A[2:n+1,5] = float.(p)
A[2:n+1,6] = t
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
    FID = open(file, "w")
    surreal2dag(FID, x[i])
    close(FID)
    run(`dot -Tpdf -O $file`)
end


# look at this in more detail 
y = convert(SurrealFinite, 1) + convert(SurrealFinite, 2)*convert(SurrealFinite,2)
float.(y.L)

# also should do   2^k * 2^(-k), and see how far we get

