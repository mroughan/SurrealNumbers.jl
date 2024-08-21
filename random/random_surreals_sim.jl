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

for (i1,α) in enumerate(Alphas)
    for (i2,λ) in enumerate(Lambdas)
        for (i3,n) in enumerate(Ns)
            for (i4,g_0) in enumerate(Gs)
                clearcache() # avoid overcommitting everything to memory
                println(" α=$α, λ=$λ, n=$n, g_0=$g_0")
                
                for seed = (s0+1) : (s0+Nseeds)
                    println("   seed=$seed")
                    global filename = @sprintf("Data/Raw/random_surreal_sim_%s_%0.2f_%0.2f_%05d_%02d_%03d_%03d_proc.csv", distributions, λ,  α , n, g_0, m, seed)
                    if isfile(filename)
                        continue  # don't redo a sim that's already done
                    end
                    rng = StableRNG(seed);
                    global output, stats, mean_stats, max_stats, min_stats = 
                            m_step_generation( n, m, g_0, α, λ, seed; V=false, parent_dist = parent_dist,  split_dist = split_dist)

                    # show the evolution for each sim
                    global DF = DataFrame( 
                            "Iteration" => collect(0:m),
                            "Mean Generation" => [s.generation for s in mean_stats],
                            "Max Generation" => [s.generation for s in max_stats],
                            "Mean Nodes" => [s.nodes for s in mean_stats],
                            "Max Nodes" => [s.nodes for s in max_stats],
                            "Mean Edges" => [s.edges for s in mean_stats],
                            "Max Edges" => [s.edges for s in max_stats],
                            # "Mean Value" => [s.value for s in mean_stats],
                            # "Max Value" => [s.value for s in max_stats],
                            # "Min Value" => [s.value for s in min_stats],
                    )
                    CSV.write(filename, DF)

                    # show the final distribution
                    stats = dag_stats.(output; V=true) # redo stats for final generation to get values
                    global filename2 = @sprintf("Data/Raw/random_surreal_sim_%s_%0.2f_%0.2f_%05d_%02d_%03d_%03d_final.csv", distributions, λ,  α , n, g_0, m, seed)
                    global DF2 = DataFrame( 
                            "Generation" => [s.generation for s in stats],
                            "Nodes" => [s.nodes for s in stats],
                            "Edges" => [s.edges for s in stats],
                            "Value" => [Float64(s.value) for s in stats],
                            "Integer" => [Int64(isinteger(s.value)) for s in stats],
                            "Gap" => [Float64(s.gap) for s in stats],
                            "Breadth" => [Float64(s.breadth) for s in stats],
                            "Parents" => [length(parents(o)) for o in output ], # wasn't included when I rain most of the existing sims
                            "Max parents" => [s.max_parents for s in stats]
                    )
                    CSV.write(filename2, DF2)
                    # I would like to add to this
                    #    1 = number of unique forms (other than canonical-zero)
                    #              but ideally, would do this by writing out a hash for each ???
                    #              adding a thing in here would require redoing the syms ...

                end
            end
        end
    end
end
