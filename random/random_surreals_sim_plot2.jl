# do some plots from the output of random_surreals_sim.jl
#   generation vs nodes and edges plots
using DataFrames
using CSV
using Printf

using Random
using Distributions
using GLM

using Plots
using PyPlot
pyplot()

input_dir = "Data/Compiled/"
plot_dir = "Plots/"

include("random_surreals.jl")
include("utilities.jl")
include("parameters.jl") 
include("plot_parameters.jl") 

α = 0.8
λ = 3.5
g_0 = 1
n = 4000
Nseeds = 30
PopulationSize = n * Nseeds

base_file_name = joinpath( input_dir, @sprintf("random_surreal_sim_%s_%0.2f_%0.2f_%05d_%02d_%07d", distributions, λ,  α , n, g_0, PopulationSize))

# scatter plot of nodes vs gen and edges vs gen
global input_file = "$(base_file_name).csv" 
global DF = CSV.read(input_file, DataFrame)
x_gen = DF[:, Symbol("Generation")]
y_nodes = DF[:, Symbol("Nodes") ]
PyPlot.figure("Gen v Nodes")
global p12 = Plots.scatter( x_gen, y_nodes;  label="", alpha=ma, marker = marker[1], markersize = sms, markercolor = mc[1], markerstrokecolor = mc[1])
Plots.plot!(p12; xlabel="Generation (g)", ylabel="Nodes (n)", legend=:topleft)

df = DataFrame( x=x_gen, y=y_nodes )
ft1 = glm( @formula( y ~ 1 + x + x^2), df, Normal())
ft2 = glm( @formula( y ~ 1 + x^2), df, Normal())
# ft = fit( LinearModel, @formula( y ~ x^2), df)

mg = maximum(x_gen)
xx = collect(0:mg)
a1 = coef(ft1)
a2 = coef(ft2)

average_nodes = Vector{Float64}(undef, length(xx))
for (i,x) in enumerate(xx)
    average_nodes[i] = mean( y_nodes[ x_gen .== x] )
end 
Plots.scatter!(p12, xx, average_nodes; label="mean", alpha=1.0, marker = marker[2], markercolor = mc[2], markerstrokecolor = mc[2])
# Plots.plot!( p12, xx, repeat([ n*(1 - exp(-λ)) ], length(xx)); linestyle=:dot, label="maximum size")

# f1 = 0.5
# f2 = 1.0
# f3 = 0.25
# Plots.plot!(p12, 0:mg, f1*(0:mg).^2; )
#Plots.plot!(p12, 0:mg, f2*(0:mg).^2; )
# Plots.plot!(p12, 0:mg, f3*(0:mg).^2; )
Plots.plot!(p12, xx, @. a1[1] + a1[2]*xx + a1[3]*xx^2; label="quadratic")
# Plots.plot!(p12, xx, @. a2[1]*xx^2; label="pure quad")
plot_file = joinpath( plot_dir, @sprintf("random_surreal_sim_%s_%0.2f_%0.2f_%05d_%02d_%07d_gen_node.pdf", distributions, λ,  α , n, g_0, PopulationSize))
Plots.savefig(p12, plot_file)
println("   output $plot_file")

##########################################################################33
y_edges = DF[:, Symbol("Edges") ]
PyPlot.figure("Gen v Edges")
global p13 = Plots.scatter( x_gen, y_edges;  label="", alpha=ma, marker = marker[1], markersize = sms, markercolor = mc[1], markerstrokecolor = mc[1] )
Plots.plot!(p13; xlabel="Generation (g)", ylabel="Edges (e)", legend=:topleft)

average_edges = Vector{Float64}(undef, length(xx))
for (i,x) in enumerate(xx)
    average_edges[i] = mean( y_edges[ x_gen .== x] )
end  

df = DataFrame( x=x_gen, y=y_edges )
ft3 = glm( @formula( y ~ 0 + x + x^2), df, Normal())
ft4 = glm( @formula( y ~ 0 + x^2), df, Normal())
a3 = coef(ft3)
a4 = coef(ft4)
a4_cm = (a4[1] - 1.96*stderror(ft4)[1])
a4_cp = (a4[1] + 1.96*stderror(ft4)[1])
Plots.scatter!(p13, xx, average_edges; label="mean", alpha=1.0, marker = marker[2], markercolor = mc[2], markerstrokecolor = mc[2])
Plots.plot!(p13, xx, @. a3[1]*xx + a3[2]*xx^2; label="quadratic")
# Plots.plot!(p13, xx, @. a4[1]*xx^2; label="pure quad")
# Plots.plot!(p13, xx, @. a4_cm*xx^2; linestyle=:dash, label="")
# Plots.plot!(p13, xx, @. a4_cp*xx^2; linestyle=:dash, label="")

plot_file = joinpath( plot_dir, @sprintf("random_surreal_sim_%s_%0.2f_%0.2f_%05d_%02d_%07d_gen_edges.pdf", distributions, λ,  α , n, g_0, PopulationSize))
Plots.savefig(p13, plot_file)
println("   output $plot_file")

##########################################################################33
y_edges = DF[:, Symbol("Edges") ]
PyPlot.figure("Nodes v Edges")
global p23 = Plots.scatter( y_nodes, y_edges;  label="", alpha=ma, marker = marker[1], markersize = sms, markercolor = mc[1], markerstrokecolor = mc[1] )
Plots.plot!(p23; xlabel="Nodes (n)", ylabel="Edges (e)")

max_nodes = maximum(y_nodes)
yy = collect(0:max_nodes)

# Plots.plot!(p23, yy, yy*λ; linestyle=:dash, label=L"\lambda n")
Plots.plot!(p23, yy, yy .- 1; linestyle=:dash, label=L"e = n-1")
Plots.plot!(p23; legend=:topleft, xlim=(0,250), ylim=(0,1000) ) # α = 0.8, λ = 3.5

df = DataFrame( x=y_nodes, y=y_edges )
ft6 = glm( @formula( y ~ 1 + x ), df, Normal())
slope = round( coef(ft6)[2]; digits = 2)
zz = coef(ft6)[1] .+ yy*coef(ft6)[2]
Plots.plot!(p23, yy, zz; linestyle=:dot, label=LaTeXString("e = $(slope) \$n\$"), linewidth=2*lw, color=:black)

plot_file = joinpath( plot_dir, @sprintf("random_surreal_sim_%s_%0.2f_%0.2f_%05d_%02d_%07d_nodes_edges.pdf", distributions, λ,  α , n, g_0, PopulationSize))
Plots.savefig(p23, plot_file)
println("   output $plot_file")


PyPlot.figure("Nodes v Edges (residual)")
global p23b = Plots.scatter( y_nodes, y_edges .- coef(ft6)[1] .- y_nodes*coef(ft6)[2];  label="", alpha=ma, marker = marker[1], markersize = sms, markercolor = mc[1], markerstrokecolor = mc[1] )
Plots.plot!(p23b; xlabel="Nodes (n)", ylabel="Residual Edges (minus linear model)")

p12

