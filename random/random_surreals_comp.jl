# compile the Raw simulation data from all seeds
#    # final state distributions, record across all seeds
#    # average state as number of iterations progress (across all seeds)
using DataFrames
using CSV
using Printf

include("random_surreals.jl")
include("utilities.jl")
include("parameters.jl") 

input_dir = "Data/Raw"
output_dir = "Data/Compiled"


Ns = [4000]
Alphas = [0.4]

for (i1,α) in enumerate(Alphas)
    for (i2,λ) in enumerate(Lambdas)
        for (i3,n) in enumerate(Ns)
            for (i4,g_0) in enumerate(Gs)
                println(" α=$α, λ=$λ, n=$n, g_0=$g_0")
                
                # output statistics
                global input_file_pattern = 
                    Regex(@sprintf("random_surreal_sim_%s_%0.2f_%0.2f_%05d_%02d_%03d_(\\d+)_final.csv", 
                                     distributions, λ,  α, n, g_0, m))
                global files = list_files(; path=input_dir, pattern=input_file_pattern, join=true)
                if isempty(files)
                    println("    no output data files")
                else    

                    global PopulationSize = n * length(files)
                    
                    global DF = Vector{DataFrame}(undef, 0)
                    for (i,file) in enumerate(files)
                        println("    reading $file")
                        if i==1
                            DF = CSV.read(file, DataFrame)
                        else    
                            append!(DF, CSV.read(file, DataFrame))
                        end
                    end
                    DF[:,"Log2_Denominator"] = [ Int64.(round.(log2.(Rational(r).den))) for r in DF[:,"Value"] ]
                    DF[:,"Canonical-Zeros"] = [ Int64(n==1) for n in DF[:,"Nodes"]]
                    global output_file = @sprintf("%s/random_surreal_sim_%s_%0.2f_%0.2f_%05d_%02d_%07d.csv", output_dir, distributions, λ,  α , n, g_0, PopulationSize)
                    CSV.write(output_file, DF)

                    for k in ["Generation","Nodes","Edges","Log2_Denominator","Canonical-Zeros"]
                        global var = DF[:,k]
                        global output_file = @sprintf("%s/random_surreal_sim_%s_%0.2f_%0.2f_%05d_%02d_%07d_%s.csv", output_dir, distributions, λ,  α , n, g_0, PopulationSize, k)
                        global cn,cx = count_n(var)
                        global DF2 = DataFrame(
                            "$k" => collect(cx),
                            "f($k)" => collect(cn),
                            "Population Size" => repeat( [PopulationSize], length(cx) )
                        )
                        CSV.write(output_file, DF2)
                    end 
                end

                # time-process statistics
                global input_file_pattern = Regex(@sprintf("random_surreal_sim_%s_%0.2f_%0.2f_%05d_%02d_%03d_(\\d+)_proc.csv", distributions, λ,  α , n, g_0, m))
                global files = list_files(; path=input_dir, pattern=input_file_pattern, join=true)
                if isempty(files)
                    println("    no process data files")
                else    
                    
                    global DF3 = Vector{DataFrame}(undef, length(files))
                    for (i,file) in enumerate(files)
                        println("    reading $file")
                        global DF3[i] = CSV.read(file, DataFrame)
                    end 
                    global DF_mean = compile_dataframes( DF3 )
                    global DF_std = compile_dataframes( DF3; fn=std )

                    for k in ["Generation","Nodes","Edges"]
                        for j in ["Max", "Mean"]
                            tmp = "$j $k"
                            global DF5 = DataFrame("Mean $tmp" => DF_mean[!,tmp], "Std $tmp" => DF_std[!,tmp])
                            global output_file = @sprintf("%s/random_surreal_sim_%s_%0.2f_%0.2f_%05d_%02d_%07d_proc_%s_%s.csv", output_dir, distributions, λ,  α , n, g_0, PopulationSize, k, j)
                            CSV.write(output_file, DF5)
                        end
                    end
                end

            end
        end
    end
end
