using SurrealNumbers 

n = 3
x = convert.(SurrealFinite, 0:n)
A = zeros(SurrealFinite, n+1, n+1)    
for i=1:n+1
    for j=1:n+1
        A[i,j] = x[i] + x[j]
    end
end
float.(A)
A
generation.(A)
A .== canonicalise.(A)

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


n = 3
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



# analyse.(A)


# function analyse(s::SurrealFinite)
#     println("depth = ", generation(s) )
    

# end
