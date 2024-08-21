# take compiled stats, and do summary as a function of α and λ
#   A - generation, nodes, edges
#   B - convergence
#   C - coefficient of n^2 in relationship between g and n
# 
using DataFrames
using CSV
using Printf

using Random
using Distributions
using GLM

include("random_surreals.jl")
include("utilities.jl")
include("parameters.jl") 

input_dir = "Data/Compiled"
output_dir = "Data/"

n = 4000
g_0 = 1
Nseeds = 30
PopulationSize = n * Nseeds

distributions = "P-U"
measures = ["Generation", "Nodes", "Edges"]

DF = DataFrame( Alpha=Float64[], Lambda=Float64[], Generation=Float64[], Nodes=Float64[], Edges=Float64[] )
row = []
len = 0

# A -- means for generation, nodes and edges
for (i1,α) in enumerate(Alphas)
    for (i2,λ) in enumerate(Lambdas)
        global row
        global distributions, n, g_0, PopulationSize, len
        len += 1

        row = zeros(5)
        row[1] = α
        row[2] = λ
        base_file_name = joinpath( input_dir, @sprintf("random_surreal_sim_%s_%0.2f_%0.2f_%05d_%02d_%07d", distributions, λ,  α , n, g_0, PopulationSize))

        for (j,measure) in enumerate(measures)
            global DF, row, base_filename

            println(" α=$α, λ=$λ, measure=$measure")
            global input_file = "$(base_file_name)_$measure.csv" 
            println("    reading $input_file")

            global DF1 = CSV.read(input_file, DataFrame)
            global x = DF1[:, Symbol("$measure")]
            global y = DF1[:, Symbol("f($measure)")] ./ DF1[:, Symbol("Population Size") ]
        
            row[j+2] = sum( x .* y )
        end

        push!( DF, row )

    end
end


#   B - convergence
convergence_pt = Vector{Float64}(undef, len)
len = 0
measure = "Nodes"
thestat = "Mean"
for (i1,α) in enumerate(Alphas)
    for (i2,λ) in enumerate(Lambdas)
        global convergence_pt, measure, thestat
        global distributions, n, g_0, PopulationSize, len
        len += 1
        global file_name = joinpath( input_dir, @sprintf("random_surreal_sim_%s_%0.2f_%0.2f_%05d_%02d_%07d_proc_Nodes_Mean.csv", distributions, λ,  α , n, g_0, PopulationSize))
        println("    reading $input_file")
        global DF2 = CSV.read(file_name, DataFrame)
        
        global m = size(DF2,1)
        global y_symbol_1 = Symbol("Mean $thestat $measure")
        global y_symbol_2 = Symbol("Std $thestat $measure")
        global stop_state = mean( DF2[m-m_final_states+1 : m, y_symbol_1] )
        global ss = repeat( [stop_state], m )
        global rel_difference = abs.( ss .- DF2[:, y_symbol_1]) ./ stop_state
        global k = findfirst( rel_difference .< ϵ )
        if !isnothing(k) && k <= m-m_final_states+1
            convergence_pt[len] = k
        else
            convergence_pt[len] = Inf
        end

    end
end
s = Symbol("Convergence Iteration")
DF[!, s] = convergence_pt


#   C - coefficient of n^2 in relationship between g and n
coeff = Vector{Float64}(undef, len)
len = 0
for (i1,α) in enumerate(Alphas)
    for (i2,λ) in enumerate(Lambdas)
        global coeff
        global distributions, n, g_0, PopulationSize, len
        len += 1

        global input_file = joinpath( input_dir, @sprintf("random_surreal_sim_%s_%0.2f_%0.2f_%05d_%02d_%07d.csv", distributions, λ,  α , n, g_0, PopulationSize))
        println("    reading $input_file")
        global DF3 = CSV.read(input_file, DataFrame)
        x_gen = DF3[:, Symbol("Generation")]
        y_nodes = DF3[:, Symbol("Nodes") ]
        
        global df = DataFrame( x=x_gen, y=y_nodes )
        global ft1 = glm( @formula( y ~ 1 + x + x^2), df, Normal())
        # ft2 = glm( @formula( y ~ 1 + x^2), df, Normal())
        # ft = fit( LinearModel, @formula( y ~ x^2), df)
        global a1 = coef(ft1)
        coeff[len] = a1[3]
    end
end
s = Symbol("Quadratic Coefficient")
DF[!, s] = coeff

output_file = joinpath( output_dir, @sprintf("random_surreal_comp2_%s_%05d_%02d_%07d.csv", distributions, n, g_0, PopulationSize))
CSV.write(output_file, DF)

