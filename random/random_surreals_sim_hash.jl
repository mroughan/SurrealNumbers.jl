using DataFrames
using CSV
using Printf

include("random_surreals.jl")
include("parameters.jl")

# most of the globals don't have to be, but its helpful for debugging at least
output = []
stats = []
max_stats = []
min_stats = []
mean_stats = []
DF = []
gens = []
nodes = [] 
edges = []
integer_freq = 0.0

α = 0.5
λ = 3.0
g_0 = 5
n = 10000
# Nseeds = 1

for seed = (s0+1) : (s0+Nseeds)
	println("   seed=$seed")
	global filename = @sprintf("Data/Raw/random_surreal_sim_hash_%s_%0.2f_%0.2f_%05d_%02d_%03d_%03d.csv", distributions, λ,  α , n, g_0, m, seed)
	if isfile(filename)
			continue  # don't redo a sim that's already done
	end
	rng = StableRNG(seed);
	global output, stats, mean_stats, max_stats, min_stats = 
	m_step_generation( n, m, g_0, α, λ, seed; V=false, parent_dist = parent_dist,  split_dist = split_dist)
	
	# remove the zeros
	global nonzero = map(!iszero, output)
	output = output[ nonzero ]

	# calculate the hashes with some other stats, so we can do correlations
	stats = dag_stats.(output; V=true)
	global DF = DataFrame( 
			"Hash" => hash.(output),
			"Generation" => [s.generation for s in stats],
			"Nodes" => [s.nodes for s in stats],
			"Edges" => [s.edges for s in stats],
			"Value" => [Float64(s.value) for s in stats],
	)
	CSV.write(filename, DF)
	println("   wrote $filename")
end       

