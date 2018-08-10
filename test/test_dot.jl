using SurrealNumbers 
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

out_dir = "Data/"

#### PRESUMES THAT DOT (from GraphVis) is installed

# recursive construction and one and zero
x0 = SurrealFinite("0", ϕ,   ϕ)  
x1 = SurrealFinite("1", [x0], ϕ)
x11 = SurrealFinite("-1", ϕ,  [x0])
x4 = SurrealFinite( [x11],  [x1,x0])
x5 = SurrealFinite( [x11,x0,x4],  [x1])
x21  = SurrealFinite("2", [x1],  ϕ)
x22  = SurrealFinite( [x0, x1],  ϕ)
x23  = SurrealFinite( [x11,x1],  ϕ)
x24  = SurrealFinite( [x11,x0,x1],  ϕ)
x25 = [x11,x0,x1] ≀ ϕ
z = zero(x1)
z = zero(SurrealFinite)
z = one(x1)
s1 = convert(SurrealFinite, 1//2)
s2 = convert(SurrealFinite, 3//4)
s3 = convert(SurrealFinite, 1//2) + convert(SurrealFinite, 1//4) 

x3 = convert(SurrealFinite, 3)
x41 = convert(SurrealFinite, 4)
x42 = SurrealFinite( [x3],  ϕ)
x43 = convert(SurrealFinite, 2) * convert(SurrealFinite, 2)

x6 = convert(SurrealFinite, 3) * convert(SurrealFinite, 2)
x9 = convert(SurrealFinite, 3) * convert(SurrealFinite, 3)

x03 = convert(SurrealFinite, 3) - convert(SurrealFinite, 3)

x00 = SurrealFinite("0", [-x1],   [x1])  
x000 = x00 + x00

surreal2dot(STDOUT, x1)

file = "$(out_dir)test_dot_x41.dot"
FID = open(file, "w")
surreal2dot(FID, x41)
close(FID)
run(`dot -Tpdf -O $file`) 
run(`dot -Tsvg -O $file`) 

file = "$(out_dir)test_dag_x41.dot"
FID = open(file, "w")
surreal2dag(FID, x41)
close(FID)
run(`dot -Tpdf -O $file`) 
run(`dot -Tsvg -O $file`) 

file = "$(out_dir)test_dot_s2.dot"
FID = open(file, "w")
surreal2dot(FID, s2)
close(FID)
run(`dot -Tpdf -O $file`) 
run(`dot -Tsvg -O $file`) 

file = "$(out_dir)test_dag_s2.dot"
FID = open(file, "w")
surreal2dag(FID, s2)
close(FID)
run(`dot -Tpdf -O $file`) 
run(`dot -Tsvg -O $file`) 

file = "$(out_dir)test_dot_s3.dot"
FID = open(file, "w")
surreal2dot(FID, s3)
close(FID)
run(`dot -Tpdf -O $file`) 
run(`dot -Tsvg -O $file`) 

file = "$(out_dir)test_dag_s3.dot"
FID = open(file, "w")
surreal2dag(FID, s3)
close(FID)
run(`dot -Tpdf -O $file`) 
run(`dot -Tsvg -O $file`) 

file = "$(out_dir)test_dot_x43.dot"
FID = open(file, "w")
surreal2dot(FID, x43)
close(FID)
run(`dot -Tpdf -O $file`) 
run(`dot -Tsvg -O $file`) 

file = "$(out_dir)test_dag_x43.dot"
FID = open(file, "w")
surreal2dag(FID, x43)
close(FID)
run(`dot -Tpdf -O $file`) 
run(`dot -Tsvg -O $file`) 

file = "$(out_dir)test_dot_x1.dot"
FID = open(file, "w")
surreal2dot(FID,  convert(SurrealFinite, 1))
close(FID)
run(`dot -Tpdf -O $file`) 
run(`dot -Tsvg -O $file`) 

file = "$(out_dir)test_dot_x00.dot"
FID = open(file, "w")
surreal2dot(FID,  x00)
close(FID)
run(`dot -Tpdf -O $file`) 
run(`dot -Tsvg -O $file`) 

temp1 = values(ExistingSurreals)
temp1 = sort(collect(temp1))

file = "$(out_dir)test_dot_x5.dot"
FID = open(file, "w")
surreal2dot(FID, x5)
close(FID)
run(`dot -Tpdf -O $file`) 
run(`dot -Tsvg -O $file`) 

file = "$(out_dir)test_dag_x5.dot"
FID = open(file, "w")
surreal2dag(FID, x5) 
close(FID)
run(`dot -Tpdf -O $file`) 
run(`dot -Tsvg -O $file`) 
 
# file = "$(out_dir)test_dot_x6.dot"
# FID = open(file, "w")
# surreal2dot(FID, x6)
# close(FID)
# run(`dot -Tpdf -O $file`) 
# run(`dot -Tsvg -O $file`) 

file = "$(out_dir)test_dag_x6.dot"
FID = open(file, "w")
surreal2dag(FID, x6)
close(FID)
run(`dot -Tpdf -O $file`) 
run(`dot -Tsvg -O $file`) 

file = "$(out_dir)test_dag_x9.dot"
FID = open(file, "w")
surreal2dag(FID, x9) 
close(FID) 
run(`dot -Tpdf -O $file`) 
run(`dot -Tsvg -O $file`) 

file = "$(out_dir)test_dag_x1.dot"
FID = open(file, "w")
surreal2dag(FID,  convert(SurrealFinite, 1))
close(FID)
run(`dot -Tpdf -O $file`) 
run(`dot -Tsvg -O $file`) 

file = "$(out_dir)test_dag_x00.dot"
FID = open(file, "w")
surreal2dag(FID,  x00)
close(FID)
run(`dot -Tpdf -O $file`) 
run(`dot -Tsvg -O $file`) 

file = "$(out_dir)test_dag_x000.dot"
FID = open(file, "w")
surreal2dag(FID,  x000)
close(FID)
run(`dot -Tpdf -O $file`) 
run(`dot -Tsvg -O $file`) 

temp2 = values(ExistingSurreals)
temp2 = sort(collect(temp2))
