

# utility function to do quick and dirty histograms over the integers
#  could use StatsBase, but this way I avoid one more dependency
function count_n(X::Vector{T}) where T <: Integer
    # assume array of ints from 1--max_x
    #d = Dict{T, Int}()
    max_x = 0
    for x in X
        if x>max_x
            max_x = x
        end
        # if haskey(d,x)
        #    d[x] += 1
        #else
        #    d[x] = 1
        #end
    end 
    c = zeros(Int64,max_x+1)
    for x in X
        c[x+1] += 1
    end
    cx = 0:max_x;
    co = OffsetArray(c, cx)
    return co, cx
end

# compose results over a set of dataframes holding the results
function compile_dataframes( DF::Vector{DataFrame}; fn = mean) #assumes each DF in vector is the same
    df = similar(DF[1])
    for name in names(df)
        df[:,name] = fn( [d[:,name] for d in DF] )
    end
    return df
end

# create a list of files based on an input pattern
function list_files(; path::AbstractString=".", pattern::Regex=r"", join::Bool=false)::Vector{String}
    raw"""
        list_files( ; path::AbstractString=".", pattern::Regex=r"", join::Bool=false) -> Vector{String}

    Returns an array of file names in the specified directory path matching the regular expression if given.

    # Arguments:

    * `path::AbstractString`: Directory path to list files from, by default te current directory.
    * `pattern::Regex`: Only file names which match the regular expression will be returned.
    * `join::Bool=false`: When join is true, it returns joinpath(dir, name) for each name so that the returned strings are full paths

    """
    files = [file for file in readdir(path; join=join) if occursin(pattern, file)]
    return files
end
