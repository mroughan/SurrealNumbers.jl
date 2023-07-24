# generate the positive quadrant of the dyadic tree
using SurrealNumbers

out_dir = "Figs"
direction = "back" # version used in Mathematics magazine
# direction = "forward" # this is the original and default version, used in SoftwareX

n=5
x_scale = 4.0

file2 = joinpath(@__DIR__, out_dir, "dyadic_tree_quad_$n.dot")
io = open(file2, "w")
println(io, "digraph dyadic_tree {")

name_map2 = Dict{SurrealFinite, Integer}()
color = Dict{SurrealFinite, Integer}()
color_map = Dict(1=>"red3", 2=>"blue3", 3=>"green3")
C = Dict(1=>"R", 2=>"B", 3=>"G")
shape = "shape=polygon,sides=4,"
shape = "shape=ellipse,width=0.5,height=0.5,"
shape = "shape=circle,"
shape = "shape=box,"
shape = "shape=plaintext,"

z = zero(SurrealFinite)
m = 0
name_map2[z] = m
var_args = "pos=\"0,0!\""
cn = 1
println(io, "   node_$m [label=0, $shape  fontsize=36, color=\"$(color_map[cn])\", $var_args];")
x = z
color[x] = cn
    
for k=1:n
    global m
    global fs
    cn = mod(k,2) + 1
    
    fs = 36 - 4*k 

    var_args = "pos=\"-$(x_scale*k),$(k)!\""
    x = convert(SurrealFinite, -k)
    m += 1
    # println(io, "   node_$m [label=$(-k),fontsize=$fs, color=\"$(color_map[cn])\", $var_args];")
    # name_map2[x] = m
    # color[x] = cn
    # println(" k = $k, x = $x")
    # tmp = name_map2[x.R[1]]
    # println(io, "   node_$m -> node_$tmp;")
    
    var_args = "pos=\"$(x_scale*k),$(k)!\""  
    x = convert(SurrealFinite, k)
    m += 1
    println(io, "   node_$m [label=$(k), $shape  fontsize=$fs, color=\"$(color_map[cn])\", $var_args];")
    name_map2[x] = m
    color[x] = cn
    println(" k = $k, x = $x")
    tmp = name_map2[x.L[1]] 
    println(io, "   node_$m -> node_$tmp;")

    for f=2:k
        for j=(f-2)*2^(k-f+1) .+ (1:2:2^(k-f+1))
            v = j // 2^(k-f+1)
            s = string(v.num) * "/" * string(v.den) 
            x = convert(SurrealFinite, v)
            cn = setdiff([1,2,3], [color[x.L[1]], color[x.R[1]]] )[1]
            color[x] = cn
            var_args = "pos=\"$(x_scale*convert(Float64, x)),$(k)!\""
            m += 1
            println(io, "   node_$m [label=\"$(s)\", $shape  fontsize=$fs, color=\"$(color_map[cn])\", $var_args];")
            name_map2[x] = m 
            println(" k = $k, f = $f, x = $x")
            tmp = name_map2[x.L[1]]
            println(io, "   node_$m -> node_$tmp;")
            tmp = name_map2[x.R[1]]
            println(io, "   node_$m -> node_$tmp;")

            # v = -j // 2^(k-f+1)
            # s = string(v.num) * "/" * string(v.den)
            # x = convert(SurrealFinite, v)
            # cn = setdiff([1,2,3], [color[x.L[1]], color[x.R[1]]] )[1]
            # color[x] = cn
            # var_args = "pos=\"$(x_scale*convert(Float64, x)),$(k)!\""
            m += 1
            # println(io, "   node_$m [label=\"$(s)\", fontsize=$fs, color=\"$(color_map[cn])\", $var_args];")
            # name_map2[x] = m 
            # println(" k = $k, f = $f, x = $x")
            # tmp = name_map2[x.L[1]]
            # println(io, "   node_$m -> node_$tmp;")
            # tmp = name_map2[x.R[1]]
            # println(io, "   node_$m -> node_$tmp;") 
        end 
    end
end

println(io, "}")
close(io)
# run(`dot -n -Kfdp -Tpdf -O $file`) 
run(`neato -Tpdf -O $file2`) 
run(`neato -Tsvg -O $file2`) 

 
X = sort( collect(keys(color)))
for x in X
    print(" , $(C[color[x]])")
end
println()


