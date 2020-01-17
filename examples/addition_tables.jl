using SurrealNumbers 
out_dir = "Data/"

# addition table
x = convert.(SurrealFinite, [0, 1/2, 1, 2])
A = zeros(SurrealFinite, length(x), length(x))    
B = zeros(Int, length(x), length(x))    
for i=1:length(x)
    for j=1:length(x)
        A[i,j] = x[i] + x[j]
        B[i,j] = generation(x[i]) + generation(x[j])
    end
end
float.(A)
A
generation.(A)
# A .== canonicalise.(A)
iscanonical.(A)
B

a10 = convert(SurrealFinite, 3/2)
a11 = convert(SurrealFinite, 1/2) + convert(SurrealFinite, 1)
a12 = convert(SurrealFinite, 3/4) + convert(SurrealFinite, 3/4)

file = "$(out_dir)addition_ex_10.dot"
FID = open(file, "w")
surreal2dag(FID, a10)
close(FID)
run(`dot -Tpdf -O $file`)

file = "$(out_dir)addition_ex_11.dot"
FID = open(file, "w")
surreal2dag(FID, a11; direction="back")
close(FID)
run(`dot -Tpdf -O $file`)

file = "$(out_dir)addition_ex_12.dot"
FID = open(file, "w")
surreal2dag(FID, a12)
close(FID)
run(`dot -Tpdf -O $file`)


