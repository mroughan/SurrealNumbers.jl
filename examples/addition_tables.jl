using SurrealNumbers 
fig_dir = "Figs"

# addition table
x = convert.(SurrealFinite, [0, 1/2, 1, 2])
n = length(x)
A = zeros(SurrealFinite, n, n)    
B = zeros(Int, n, n)    
for i=1:n
    for j=1:n
        A[i,j] = x[i] + x[j]
        B[i,j] = generation(x[i]) + generation(x[j])
    end
end
S = dag_stats.(A)
B = [ S[i,j].generation for i=1:n, j=1:n] 
# A .== canonicalise.(A)
C = iscanonical.(A)
A = float.(A)
print("additions = "); display(A)
println()
print("generation = "); display(B)
println()
print("iscanonical = "); display(C)
println()


### some specific examples
a10 = convert(SurrealFinite, 3/2)
a11 = convert(SurrealFinite, 1/2) + convert(SurrealFinite, 1)
a12 = convert(SurrealFinite, 3/4) + convert(SurrealFinite, 3/4)

file = joinpath(@__DIR__, fig_dir, "addition_ex_10.dot")
FID = open(file, "w")
surreal2dag(FID, a10; direction="back")
close(FID)
run(`dot -Tpdf -O $file`)

file = joinpath(@__DIR__, fig_dir, "addition_ex_11.dot")
FID = open(file, "w")
surreal2dag(FID, a11; direction="back")
close(FID)
run(`dot -Tpdf -O $file`)

file = joinpath(@__DIR__, fig_dir, "addition_ex_12.dot")
FID = open(file, "w")
surreal2dag(FID, a12; direction="back")
close(FID)
run(`dot -Tpdf -O $file`)


0
