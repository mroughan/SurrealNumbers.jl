# generate the dyadic tree, suitable for output to DOT
using SurrealNumbers
# var_args = "color=red"
var_args = ""
n=3; x_scale = 2.25
# n=4; x_scale = 5.00
out_dir = "Data/"

name_map = Dict{SurrealFinite, Integer}()

file = "$(out_dir)dyadic_tree_$n.dot"

io = open(file, "w")
println(io, "digraph dyadic_tree {")
z = zero(SurrealFinite)
m = 0
name_map[z] = m
var_args = "pos=\"0,0!\""
surreal2node(io, z, 0, extra_args=var_args)
println(io, "node_day_0 [shape=none, pos=\"$(x_scale*(n+0.5)),$(0)!\", label=\"Day 0\", fontsize=18]")

for k=1:n
    println(io, "node_day_$k [shape=none, pos=\"$(x_scale*(n+0.5)),$(k)!\", label=\"Day $k\", fontsize=18]")
    
    var_args = "pos=\"-$(x_scale*k),$(k)!\""
    x = convert(SurrealFinite, -k)
    m += 1
    surreal2node(io, x, m, extra_args=var_args)
    name_map[x] = m
    println(" k = $k, x = $x")
    tmp = name_map[x.R[1]]
    println(io, "   node_$m:R -> node_$tmp [color=\"blue3\"];")

    var_args = "pos=\"$(x_scale*k),$(k)!\""
    x = convert(SurrealFinite, k)
    m += 1
    surreal2node(io, x, m, extra_args=var_args)
    name_map[x] = m
    println(" k = $k, x = $x")
    tmp = name_map[x.L[1]] 
    println(io, "   node_$m:L -> node_$tmp [color=\"red3\"];")

    for f=2:k
        for j=(f-2)*2^(k-f+1) + (1:2:2^(k-f+1))
            x = convert(SurrealFinite, j // 2^(k-f+1))
            var_args = "pos=\"$(x_scale*convert(Float64, x)),$(k)!\""
            m += 1
            surreal2node(io, x, m, extra_args=var_args)
            name_map[x] = m
            println(" k = $k, f = $f, x = $x")
            tmp = name_map[x.L[1]]
            println(io, "   node_$m:L -> node_$tmp [color=\"red3\"];")
            tmp = name_map[x.R[1]]
            
            println(io, "   node_$m:R -> node_$tmp [color=\"blue3\"];")

            x = convert(SurrealFinite, -j // 2^(k-f+1))
            var_args = "pos=\"$(x_scale*convert(Float64, x)),$(k)!\""
            m += 1
            surreal2node(io, x, m, extra_args=var_args)
            name_map[x] = m
            println(" k = $k, f = $f, x = $x")
            tmp = name_map[x.L[1]]
            println(io, "   node_$m:L -> node_$tmp [color=\"red3\"];")
            tmp = name_map[x.R[1]]
            println(io, "   node_$m:R -> node_$tmp [color=\"blue3\"];") 
        end
    end
end

println(io, "}")
close(io)
# run(`dot -n -Kfdp -Tpdf -O $file`) 
run(`neato -Tpdf -O $file`) 
run(`neato -Tsvg -O $file`) 


