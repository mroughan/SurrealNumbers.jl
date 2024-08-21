# do some plots from the output of random_surreals_sim.jl
#  this one shows the evolutions through the iterations of average (and maybe max) stats
using DataFrames
using CSV
using Printf 

using Random
using Distributions

using Plots
using PyPlot
pyplot()

input_dir = "Data/Compiled/"
plot_dir = "Plots/"

include("random_surreals.jl")
include("utilities.jl")
include("parameters.jl") 
include("plot_parameters.jl") 


α = 0.8; λ = 3.5
# α = 0.4; λ = 1.5

g_0 = 1
n = 4000
Nseeds = 30
PopulationSize = n * Nseeds

distributions = "P-U"
base_file_name = joinpath( input_dir, @sprintf("random_surreal_sim_%s_%0.2f_%0.2f_%05d_%02d_%07d_proc", 
                                                distributions, λ,  α , n, g_0, PopulationSize))

p1 = Matrix{Any}(undef,2,3)
p2 = Vector{Any}(undef,2)

# PyPlot.figure(1)
for (i,stat) in enumerate(["Mean"; "Max"])
    global p2
    p2[i]  = Plots.plot( [0], [0]; xlabel="Iteration", ylabel="$stat", legend=:right, label="")

    if stat == "Mean"
        ls = ["Generation", "Nodes", "Edges"]
    else    
        ls = ["Generation", "Nodes"]
    end 

    for (j,measure) in enumerate(ls)
        global p1, p2, DF, filename, x, y_symbol_1, y_symbol_2, stop_state, n, m, m_final_states, k, ss, rel_difference

        filename = "$(base_file_name)_$(measure)_$(stat).csv"
        println("  reading $measure $stat $filename")
        DF = CSV.read(filename, DataFrame)
        
        m = size(DF,1)
        x = 0 : m-1
        y_symbol_1 = Symbol("Mean $stat $measure")
        y_symbol_2 = Symbol("Std $stat $measure")
        p1[i,j] = Plots.plot( x, DF[:, y_symbol_1] ; label="", xlabel="iteration", ylabel="$stat -- $measure")

        Plots.plot!( p2[i], x, DF[:, y_symbol_1] ; label="$measure")

        stop_state = mean( DF[m-m_final_states+1 : m, y_symbol_1] )
        ss = repeat( [stop_state], length(x) )
        Plots.plot!( p2[i], x, ss ; label="", linestyle=:dot, linewidth=lw)

        rel_difference = abs.( ss .- DF[:, y_symbol_1]) ./ stop_state
        k = findfirst( rel_difference .< ϵ )
        if !isnothing(k)
            Plots.plot!( p2[i], [k-1; k-1], [0; stop_state] ; label="", linestyle=:dot, linewidth=2*lw)
        end

        Plots.plot!( p2[i]; xlim=(0,m-1), ylim=(0,Inf), legend=:outerright )
    end 

end 

global plot_file = joinpath( plot_dir, @sprintf("random_surreal_sim_%s_%0.2f_%0.2f_%05d_%02d_%07d_proc_3.pdf", 
                                                 distributions, λ,  α , n, g_0, PopulationSize))
Plots.savefig(p2[1], plot_file)
println("  writing file $plot_file")

global plot_file = joinpath( plot_dir, @sprintf("random_surreal_sim_%s_%0.2f_%0.2f_%05d_%02d_%07d_proc_3_max.pdf", 
                                                 distributions, λ,  α , n, g_0, PopulationSize))
Plots.savefig(p2[2], plot_file)
println("  writing file $plot_file")

p2[1]

# do a version b that compares g_0
# do a version c that plots Max and Mean on samples
# do a version d that works out convergence for all cases