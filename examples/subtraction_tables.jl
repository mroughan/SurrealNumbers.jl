using IndexedTables
using FileIO
using SurrealNumbers 
out_dir = "Data/" 

# subtraction table
x = convert.(SurrealFinite, [0, 1/2, 1, 2])
m = length(x)
A = zeros(SurrealFinite, m, m)    
B = zeros(Int, m, m)    
for i=1:m
    for j=1:m
        A[i,j] = x[i] - x[j]
        B[i,j] = generation(x[i]) + generation(x[j])
     end
end
float.(A)
A
generation.(A)
# A .== canonicalise.(A)
iscanonical.(A)
B
C = dag_stats.(A)
n = [C[i,j].nodes for i=1:size(C)[1], j=1:size(C)[2] ]

# subtraction table 2
ns1 = [4; 4; 4; 4; 4; 4; 1; 2; 3]
ns2 = [0; 1; 2; 3; 4; 5; 1; 2; 3]
n = length(ns1)
m = [10; 10; 10; 10; 10; 10; 10; 10; 10]
x = Array{SurrealFinite,1}(undef, n) 
t = Array{Float64,1}(undef, ) 
bytes = Array{Float64,1}(undef, n) 
gctime = Array{Float64,1}(undef, n)
c_plus = Array{Integer,1}(undef, n) 
c_times = Array{Integer,1}(undef, n) 
c_created = Array{Integer,1}(undef, n) 
c_equals = Array{Integer,1}(undef, n) 
c_leq = Array{Integer,1}(undef, n) 
i = 0
for i=1:n 
    println(ns1[i], " - ", ns2[i])
    t[i] = 0 
    b = 0
    g = 0
    mem = 0
    val = 0 
    for j=1:m[i]
        clearcache()
        val, time, b, g, mem = @timed convert(SurrealFinite, ns1[i]) - convert(SurrealFinite, ns2[i])
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
d = (e+1) ./ no
g = [s[i].generation for i=1:length(s)]
t = t./m[1:n]
bytes = bytes

t1 = table(@NT(i=0:n-1,
               n1=ns1, 
               n2=ns2,
               nodes=no,
               edges=ed,
               paths=p,
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

save("$(out_dir)subtraction_tables.csv", t1)




# examples with diagrams
a21 = convert(SurrealFinite, 1) - convert(SurrealFinite, 1)
file = "$(out_dir)subtraction_ex_1.dot"
FID = open(file, "w")
surreal2dag(FID, a21)
close(FID)
run(`dot -Tpdf -O $file`)
 
a22 = convert(SurrealFinite, 2) - convert(SurrealFinite, 2)
file = "$(out_dir)subtraction_ex_2.dot"
FID = open(file, "w")
surreal2dag(FID, a22)
close(FID)
run(`dot -Tpdf -O $file`)


 
a22 = convert(SurrealFinite, 2) - convert(SurrealFinite, 2)
file = "$(out_dir)subtraction_ex_2.dot"
FID = open(file, "w")
surreal2dag(FID, a22)
close(FID)
run(`dot -Tpdf -O $file`)
 
 
a23 = convert(SurrealFinite, 1/2) - convert(SurrealFinite, 1/2)
file = "$(out_dir)subtraction_ex_3.dot"
FID = open(file, "w")
surreal2dag(FID, a23)
close(FID)
run(`dot -Tpdf -O $file`)
 
 

 
a41 = convert(SurrealFinite, 4) - convert(SurrealFinite, 1)
file = "$(out_dir)subtraction_ex_41.dot"
FID = open(file, "w")
surreal2dag(FID, a41)
close(FID)
run(`dot -Tpdf -O $file`) 
run(`dot -Tsvg -O $file`)

a42 = convert(SurrealFinite, 4) - convert(SurrealFinite, 2)
file = "$(out_dir)subtraction_ex_42.dot"
FID = open(file, "w")
surreal2dag(FID, a42)
close(FID)
run(`dot -Tpdf -O $file`)
run(`dot -Tsvg -O $file`)

a43 = convert(SurrealFinite, 4) - convert(SurrealFinite, 3)
file = "$(out_dir)subtraction_ex_43.dot"
FID = open(file, "w")
surreal2dag(FID, a43)
close(FID)
run(`dot -Tpdf -O $file`)
run(`dot -Tsvg -O $file`)


a5 = convert(SurrealFinite, 3) - (convert(SurrealFinite, 1)-convert(SurrealFinite, 1))
file = "$(out_dir)subtraction_ex_5.dot"
FID = open(file, "w")
surreal2dag(FID, a5)
close(FID)
run(`dot -Tpdf -O $file`)
run(`dot -Tsvg -O $file`)



 
 
