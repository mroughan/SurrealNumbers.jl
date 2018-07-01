using SurrealNumbers 
out_dir = "Data/" 

# subtraction table
x = convert.(SurrealFinite, [0, 1/2, 1, 2])
A = zeros(SurrealFinite, length(x), length(x))    
B = zeros(Int, length(x), length(x))    
for i=1:length(x)
    for j=1:length(x)
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
 
 
