# do table views of generations for multiplication, and compare to estimated version
#    assumes "multiplication_tables.jl" has run to create its datafile

# Notes: see
#  + https://oeis.org/A047662
#  + https://oeis.org/A008288
#  + https://oeis.org/A047665
#  + Fredman, "On the complexity of maintaining an array and its partial sums"
#  + tribonacci triangle (modified, because we have inhomogenous version)
# plus a couple below
using IndexedTables
using FileIO
# using LatexPrint
# using SurrealNumbers 
out_dir = "Data"


#####
# parse the data from the multiplication table
filename = joinpath(@__DIR__, out_dir, "multiplication_tables_$(VERSION).csv")
t1 = load(filename)
t2 = table(t1)

n = 6
G = zeros(Int,n,n)
a = 0 
i = 0
j = 0
G[1,:] = 1:n
G[:,1] = 1:n
for i=1:n
    for j=1:n
        a = filter(p -> eval(parse(p.n2)) == j, filter(p -> eval(parse(p.n1)) == i, t2))
        if !isempty(a)
            G[i,j] = columns(a, :generation)[1]
            G[j,i] = columns(a, :generation)[1]
        end
    end
end
print("G = "); display(G)
println()

# now do the estimates
m = 6
F = zeros(Int128,m,m)
F[1,:] = 1:m
F[:,1] = 1:m
for i=2:m
    for j=2:m
        F[i,j] = F[i-1,j] + F[i,j-1] + F[i-1,j-1] + 1
    end
end
print("F = "); display(F)
println()

#####
# approximate size for square
N = zeros(Float64,m)
N2 = zeros(Float64,m)
# λ = 3 + sqrt(2) # Mistake in Fredman
# λ = 5.83        # approximate value
λ = 3 + 2*sqrt(2) 
c = sqrt(λ/π) / 2^(9/4)
c2 = sqrt(8+6*sqrt(2))/(8*sqrt(π)) 
for k=1:m
    N[k] = c * λ^k / sqrt(k)
end
 
C = [ [A[k,k], N[k]] for k=1:m ]
print("C = "); display(C)
println()

#####
# generation for powers of 2
n = 0:6
theta = 1.597910218031873178338070118157
g_2_n = floor.(Int, theta.^(2.^n) ) 
# or f_n = f_n*(f_n + 1)
#    g(2^n) = f_n
# https://oeis.org/A007018
#

