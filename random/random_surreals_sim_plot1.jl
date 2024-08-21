# do some plots from the output of random_surreals_sim.jl
#  mostly testing ideas, so may not work now
using DataFrames
using CSV
using Printf 

using Random
using Distributions

using Plots
using PyPlot
pyplot()

path = "Data/Raw"
# list_files(; path=path, pattern=r".*_gen.csv")


α = 0.9
λ = 0.5
g_0 = 1
n = 250
Nseeds = 10
PopulationSize = n * Nseeds
pattern = Regex(@sprintf("random_surreal_sim_%0.2f_%0.2f_%05d_%02d_%07d_(\\d+).csv", λ,  α , n, g_0, PopulationSize))
files = list_files(; path=path, pattern=pattern, join=true)

DF = Vector{DataFrame}(undef, length(files))
for (i,file) in enumerate(files)
    println("  reading $file")
    global DF[i] = CSV.read(file, DataFrame)
end 

df = compile_dataframes( DF )
PyPlot.figure(1)
p1 = Plots.plot(df[:, :Iteration], df[:, Symbol("Mean Generation")]; ylim=(0, 35) )
# p1 = Plots.plot(df[:, :Iteration], df[:, Symbol("Mean Nodes")]; ylim=(0, 35) )
for (i,file) in enumerate(files)
    println("  reading $file")
    global p1
    Plots.plot!(p1, DF[i][:, :Iteration], DF[i][:, Symbol("Mean Generation")] ; alpha=0.2, label="" )
    # Plots.plot!(p1, DF[i][:, :Iteration], DF[i][:, Symbol("Mean Nodes")] ; alpha=0.2, label="")
    # Plots.plot!(p1, DF[i][:, :Iteration], DF[i][:, Symbol("Mean Edges")] ; alpha=0.2, label="" )
    
    # Plots.plot!(p1, DF[i][:, :Iteration], DF[i][:, Symbol("Max Nodes")] ; alpha=0.2, label="" )
    
end 
x = df[:, :Iteration][1:14]
y = df[:, Symbol("Mean Nodes")][1] .+ x
Plots.plot!(p1, x, y ; label="")
    
p1

