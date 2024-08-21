# do some plots for hash functions
using DataFrames
using CSV
using Printf

using Random
using Distributions
using Statistics
using HypothesisTests
# using GLM

using Plots
using PyPlot
using StatsPlots
using LaTeXStrings
pyplot()

input_dir = "Data/Raw/"
plot_dir = "Plots/"

include("random_surreals.jl")
include("utilities.jl")
include("parameters.jl") 
include("plot_parameters.jl") 

α = 0.5
λ = 3.0
g_0 = 5
n = 10000
Nseeds = 30

global DF = Vector{DataFrame}(undef, 0)
global input_file_pattern = Regex(@sprintf("random_surreal_sim_hash_%s_%0.2f_%0.2f_%05d_%02d_%03d_(\\d+).csv", distributions, λ,  α, n, g_0, m))
global files = list_files(; path=input_dir, pattern=input_file_pattern, join=true)
if !isempty(files)
    global PopulationSize = n * length(files)
    for (i,file) in enumerate(files)
        println("    reading $file")
        global DF
        global T = Dict( "Hash"=>UInt64, "Nodes" => Int64, "Edges" => Int64, "Value" => Float64)
        if i==1
            DF = CSV.read(file, DataFrame; types=T)
        else    
            append!(DF, CSV.read(file, DataFrame))
        end
    end
    

end

# count unique hash values at each generation
hash_values = DF[:, Symbol("Hash")]
gen =  DF[:, :Generation]
m_g = maximum(gen)
U = Vector{Float64}(undef,m_g)
N = Vector{Float64}(undef,m_g)
for i=1:m_g
    N[i] = length( hash_values[ gen .== i ] )
    U[i] = length(unique( hash_values[ gen .== i ] ))
end
prop_unique = U ./ N
p1 = Plots.plot( 1:m_g, prop_unique; label="", color=:darkblue, xlabel="generation", ylabel="proportion unique" )
Plots.scatter!( 1:m_g, prop_unique; label="", color=:darkblue, alpha=0.99, )
total_unique = length( unique( hash_values ) ) 
total = length( hash_values )
println("   total proportion unique (all include) = $(total_unique/total)")

# remove the smaller generations, as these are largely non-random
threshold = 3
filter!( row -> row.Generation > threshold, DF)
hash_values = DF[:, Symbol("Hash")]
total_unique = length( unique( hash_values ) ) 
total = length( hash_values )
println("   total proportion unique (g>3) = $(total_unique/total)")

# first consider uniformity
a = 4
b = 2^a   # b=32 has many zeros
# b = 20    # spikes in regular pattern 
# looks like there are a preponderance of 0 mod 4 entries
h_a = mod.( hash_values, b)
h1 = histogram( h_a ; normalize=true, bins=b, label="b=$b")

global plot_file = joinpath( plot_dir, @sprintf("random_surreal_sim_%s_%0.2f_%0.2f_%02d_hash_uniformity_%03d.pdf", 
                                                    distributions, λ, α, g_0, b))
Plots.savefig(h1, plot_file)
println("  output $plot_file")

# b1 = @df DF groupedbar(:method, [:julia :matlab], size=(400,250),
#          width=0.5, 
#          alpha = [0.8  0.9],
#          # alpha = 0.8,
#          # fillcolor = :match,
#          # seriescolor = :lightblue,
#          # color = [:tomato, :crimson, :firebrick, :darkred],
#          # fill = [:tomato, :crimson, :firebrick, :darkred],
#          color = [:darkblue :darkred],
#          label = ["Julia" "Matlab"],
#          xlabel="Algorithm", ylabel="time (ns)",
#          stacked=false,
#          legend=(:topleft),
#          )


# formal test for uniformity isn't needed as it cleary doesn't work. 


# correlations
p = Dict{AbstractString, Any}("Generation"=>0, "Nodes"=>0, "Edges"=>0)
correlation = Dict{AbstractString, Any}()
pvals = Dict{AbstractString, Any}()
for current in keys(p)
    global p
    global S = Symbol(current)
    global x = DF[:, S]
    global mn = Vector{Float64}(undef,maximum(x))
    for i=1:maximum(x)
        mn[i] = mean( hash_values[ x .== i ] )
    end
    global p[current] = Plots.scatter( x, hash_values;  label="", alpha=0.1, marker = marker[1], markersize = sms, markercolor = mc[1], markerstrokecolor = mc[1])
    Plots.plot!(p[current], 1:maximum(x), mn; label="mean")
    Plots.plot!(p[current]; xlabel=current, ylabel="Hash Value", legend=:topleft)
    Plots.plot!(p[current]; yticks=:none )

    # Pearson correlation coefficient and test
    global correlation[current] = Statistics.cor( x, hash_values )
    global pvals[current] = HypothesisTests.pvalue( CorrelationTest(x,hash_values) )

    global plot_file = joinpath( plot_dir, @sprintf("random_surreal_sim_%s_%0.2f_%0.2f_%02d_hash_%s.pdf", 
                                                     distributions, λ, α, g_0, current))
    Plots.savefig(p[current], plot_file)
    println("  output $plot_file")
end
# have to do something different for, "Value"=>0,

h1
p["Generation"]

DF_out = DataFrame( :Variable => ["Generation", "Nodes"], 
                    :Correlation => [0.0, 0.0], 
                    Symbol("P-value") => [0.0, 0.0])
for i=1:size(DF_out,1)
    global DF_out[i,:Correlation] = correlation[ DF_out[i,:Variable] ]
    global DF_out[i,Symbol("P-value")] = pvals[ DF_out[i,:Variable] ]
end
DF_out
