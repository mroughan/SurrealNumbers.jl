Code for generating ensembles of random surreals and their summary
stats.

## Generation

+ parameters.jl           -- parameter ranges to test
+ random_surreals.jl      -- code to create the ensemble
+ random_surreals_test.jl -- examples of generation code
+ random_surreals_sim.jl  -- main simulation over parameters
+ random_surreals_sim_hash.jl -- generate data for testing hash functions

## Processing

+ utilities.jl                  -- some utility functions for summarising (and plotting)  
+ random_surreals_comp.jl  -- compile the raw sim data and create distributions and averages
+ random_surreals_comp2.jl -- further summarize stats as functions of parameters

## Plotting

+ plot_parameters.jl               -- parameters for plots
+ random_surreals_sim_plot2.jl     -- generation vs number nodes and edges
+ random_surreals_sim_plot3.jl     -- show evolution of averages (maxes) of various stats
+ random_surreals_sim_plot3a.jl    -- show evolution of averages, and compares size of population
+ random_surreals_sim_plot3b.jl    -- show evolution of averages, and compares initial population
+ random_surreals_sim_plot4b.jl    -- plots over the final distribution
+ random_surreals_sim_plot5.jl     -- summaries of random_surreals_comp2.jl, showing various features as function of parameters
+ random_surreals_sim_hash_plot.jl -- hash function plots

NB there are a lot of output files, so they are omitted, so

+ output data should go in a directory ```Data/``` that needs creation (with subdirs ```Raw/``` and ```Compiled/```)
+ output plots should go in a directory called ```Plot/s``` which needs creation


