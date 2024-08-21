# summary plots of output of comp2
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

using Plots
using PyPlot
pyplot()

input_dir = "Data/"
plot_dir = "Plots/"

include("random_surreals.jl")
include("utilities.jl")
include("parameters.jl") 
include("plot_parameters.jl") 

n = 4000
g_0 = 1
Nseeds = 30
PopulationSize = n * Nseeds
distributions = "P-U"
file_name = joinpath( input_dir, @sprintf("random_surreal_comp2_%s_%05d_%02d_%07d.csv", 
                                            distributions, n, g_0, PopulationSize))
DF = CSV.read(file_name, DataFrame)

sym_a = Symbol("Alpha")
sym_l = Symbol("Lambda")
sym_g = Symbol("Generation")
sym_n = Symbol("Nodes")
sym_e = Symbol("Edges")
sym_c = Symbol("Convergence Iteration")
sym_q = Symbol("Quadratic Coefficient")

filter!(row -> isfinite(row."Convergence Iteration"), DF)

DF_a = groupby(DF, sym_a) 
DF_l = groupby(DF, sym_l) 

Alphas = unique( DF[:,sym_a])
Lambdas = unique( DF[:,sym_l])

##########################################################33
# A -  
pA1 = Plots.plot( ; xlabel="λ", ylabel="Average Generation"  )
for (i,α) in enumerate(Alphas)
    global x = DF_a[i][:,sym_l]
    Plots.plot!( pA1, x, DF_a[i][:,sym_g]; label="α=$α" )
    expected = [generation_exp.( α, a)[1] for a in x]
    Plots.plot!( pA1, x, expected; label="", linestyle=:dash )
end
Plots.plot!( pA1; legend=:outerright)
plot_file = joinpath( plot_dir, @sprintf("random_surreal_comp2_%s_%05d_%02d_%07d_gen_a.pdf", 
                                            distributions, n, g_0, PopulationSize))
Plots.savefig(pA1, plot_file)
println("   output $plot_file")

pA2 = Plots.plot( ; xlabel="α", ylabel="Average Generation"  )
for (i,λ) in enumerate(Lambdas)
    global x = DF_l[i][:,sym_a]
    Plots.plot!( pA2, x, DF_l[i][:,sym_g]; label="λ=$λ" )
    expected = [generation_exp.( a, λ )[1] for a in x]
    Plots.plot!( pA2, x, expected; label="", linestyle=:dash )
end
Plots.plot!( pA2; legend=:outerright)
plot_file = joinpath( plot_dir, @sprintf("random_surreal_comp2_%s_%05d_%02d_%07d_gen_l.pdf", 
                                            distributions, n, g_0, PopulationSize))
Plots.savefig(pA2, plot_file)
println("   output $plot_file")

##########################################################33
# B -  
pB1 = Plots.plot( ; xlabel="λ", ylabel="Converged Iteration"  )
for (i,α) in enumerate(Alphas)
    Plots.plot!( pB1, DF_a[i][:,sym_l], DF_a[i][:,sym_c]; label="α=$α" )
end
Plots.plot!( pB1; legend=:outerright)
plot_file = joinpath( plot_dir, @sprintf("random_surreal_comp2_%s_%05d_%02d_%07d_converge_a.pdf", 
                                            distributions, n, g_0, PopulationSize))
Plots.savefig(pB1, plot_file)
println("   output $plot_file")


pB2 = Plots.plot( ; xlabel="α", ylabel="Converged Iteration", yscale=:log10  )
for (i,λ) in enumerate(Lambdas)
    Plots.plot!( pB2, DF_l[i][:,sym_a], DF_l[i][:,sym_c]; label="λ=$λ" )
end
x_a = ( DF[:, sym_a] )
y_c = log.( DF[:, sym_c] )
global df = DataFrame( x=x_a, y=y_c )
global ft1 = glm( @formula( y ~ 1 + x ), df, Normal())
global a1 = coef(ft1)
fitted_line(x, a) = exp(a[1]) .* exp.(a[2] .* x)

# x = [Alphas; 0.9]
x = [0.05; Alphas; 0.85]
Plots.plot!( pB2, x, fitted_line(x, a1); linestyle=:dash, 
                       label=@sprintf("%.1f exp(%.1f α)", exp(a1[1]), a1[2]), linewidth=2*lw, alpha=0.7)
# annotate!( pB2, x[end], fitted_line(x[end], a1), 
#                 Plots.text( @sprintf("%.1f", fitted_line(x[end], a1)), 16, :left) )

Plots.plot!( pB2; legend=:outerright)
plot_file = joinpath( plot_dir, @sprintf("random_surreal_comp2_%s_%05d_%02d_%07d_converge_l.pdf", 
                                            distributions, n, g_0, PopulationSize))
Plots.savefig(pB2, plot_file)
println("   output $plot_file")

##########################################################33
# C - 
pC1 = Plots.plot( ; xlabel="λ", ylabel="Quadratic Coefficient", yscale=:log10  )
for (i,α) in enumerate(Alphas)
    Plots.plot!( pC1, DF_a[i][:,sym_l], DF_a[i][:,sym_q]; label="α=$α" )
end
Plots.plot!( pC1; legend=:outerright)
plot_file = joinpath( plot_dir, @sprintf("random_surreal_comp2_%s_%05d_%02d_%07d_quad_a.pdf", 
                                            distributions, n, g_0, PopulationSize))
Plots.savefig(pC1, plot_file)
println("   output $plot_file")


pC2 = Plots.plot( ; xlabel="α", ylabel="Quadratic Coefficient", yscale=:log10  )
for (i,λ) in enumerate(Lambdas)
    Plots.plot!( pC2, DF_l[i][:,sym_a], DF_l[i][:,sym_q]; label="λ=$λ" )
end
Plots.plot!( pC2; legend=:outerright)
plot_file = joinpath( plot_dir, @sprintf("random_surreal_comp2_%s_%05d_%02d_%07d_quad_l.pdf", 
                                            distributions, n, g_0, PopulationSize))
Plots.savefig(pC2, plot_file)
println("   output $plot_file")

