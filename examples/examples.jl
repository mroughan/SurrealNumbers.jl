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
A .== canonicalise.(A)
iscanonical.(A)
B

a10 = convert(SurrealFinite, 3/2)
a11 = convert(SurrealFinite, 1/2) + convert(SurrealFinite, 1)
file = "$(out_dir)addition_ex_10.dot"
FID = open(file, "w")
surreal2dag(FID, a10)
close(FID)
run(`dot -Tpdf -O $file`)
file = "$(out_dir)addition_ex_11.dot"
FID = open(file, "w")
surreal2dag(FID, a11)
close(FID)
run(`dot -Tpdf -O $file`)


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
A .== canonicalise.(A)
iscanonical.(A)
B



n = 2
x = convert.(SurrealFinite, (0:2^n)/2^n)
A = zeros(SurrealFinite, 2^n+1, 2^n+1)    
for i=1:2^n +1
    for j=1:2^n +1
        A[i,j] = x[i] + x[j]
    end
end
float.(A)
A
generation.(A)
A .== canonicalise.(A)
iscanonical.(A)

n = 1
x = convert.(SurrealFinite, 0:n)
A = zeros(SurrealFinite, n+1, n+1)    
for i=1:n+1
    for j=1:n+1
        A[i,j] = x[i] * x[j]
    end
end
float.(A)
generation.(A)
A .== canonicalise.(A)
iscanonical.(A)






# analyse.(A)


# function analyse(s::SurrealFinite)
#     println("depth = ", generation(s) )
    

# end
