using Plots
using PyPlot
using Printf
pyplot()

include("random_surreals.jl")

n = 1000; # population size per iteration
gen = 6;

# λ = 3.0; α = 0.8; m = 40;
λ = 3.0; α = 0.5; m = 20;

Ns = 10
Mean_stats = Array{MeanDAGstats}(undef, Ns, m+1)
Max_stats = Array{MeanDAGstats}(undef, Ns, m+1)
Min_stats = Array{MeanDAGstats}(undef, Ns, m+1)

output = []
stats = []
max_stats = []
min_stats = []
mean_stats = []
for seed = 1:Ns
    rng = StableRNG(seed);
    global output, stats, mean_stats, max_stats, min_stats = 
              m_step_generation( n, m, gen, α, λ, seed)
    global Mean_stats[seed,:] = mean_stats
    global Max_stats[seed,:] = max_stats
    global Min_stats[seed,:] = min_stats
end

out_vals = Float64.(convert.(Rational, output))
out_g = generation.( output) 
mn = mean(out_vals)
mx = maximum(out_vals)
v = var(out_vals)
sd = sqrt(v)
si = sum(isinteger.(out_vals))
pi = si / n

# mg = mean(out_g)
cn,cx = count_n(out_g)

# Dout = Poisson( mean_stats[end].generation )
# pn = Int64.(round.(n*pdf.( Dout, cx )))

geometric_appox = λ .* (1 .- α) .* α .^ cx
pn = Int64.(round.(n .* geometric_appox  ))

pmf = generation_pmf.( cx; α=α, λ=λ)
pm = Int64.(round.(n .* pmf))

[cn'; pn'; pm']'

x = collect(cx)

pnx = 8
p1 = Plots.plot( x, collect(cn); label="Empirical" )
Plots.plot!(p1, x[pnx:end], collect(pn)[pnx:end]; label="Geometric approximation" )
Plots.plot!(p1, x, collect(pm); label="PMF" )
Plots.xlabel!(p1, "k")
Plots.ylabel!(p1, "P(X = k)")

filename = @sprintf( "Data/generation_%0.2f_%0.2f_%02d_%05d_%02d.svg", λ,  α , m, n, gen)
Plots.savefig(p1, filename)


mx = [0; collect(1:m)]
p2 = Plots.plot( mx,  [ ms.generation for ms in mean_stats]; label="generation" )
Plots.plot!( p2, mx, [ ms.nodes for ms in mean_stats]; label="nodes" )
Plots.plot!( p2, mx, [ ms.edges for ms in mean_stats]; label="edges" )
Plots.xlabel!(p2, "iteration")
Plots.ylabel!(p2, "mean (population)")
# Plots.plot!( p2, mx, [ ms. for ms in mean_stats] )
