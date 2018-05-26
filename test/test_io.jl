using SurrealNumbers

out_dir = "Data/"
 

# MISTAKE! a1 = convert(SurrealFinite, "{1//2| \phi, 1}")

a2 = convert(SurrealFinite, "{1//2 | 2,  1}")
a3 = convert(SurrealFinite, "{-1//2 | {|\phi}}")
a4 = convert(SurrealFinite, "{-1//2 | 2.0, 3}")

X = read("$(out_dir)test_surreals.dat", SurrealFinite, 4)


a4 = convert(SurrealFinite, 13//16)
a5 = convert(SurrealFinite, -2//16)
io_test_file = "$(out_dir)test_surreals_io.dat"
io = open(io_test_file, "w")
pf(io, a4)
pf(io, a5)
println(io, expand(a4))
println(io, expand(a5))
close(io)

Y = read(io_test_file, SurrealFinite, 5)

surreal2tex(a4)

io_test_file = "$(out_dir)test_surreals_io.tex"
io = open(io_test_file, "w")
surreal2tex(io, a4; level = 0)
surreal2tex(io, a4; level = 1)
surreal2tex(io, a4; level = 2)
close(io)

x = SurrealFinite([a4], [])
surreal2tex(x; level = 0)
surreal2tex(x; level = 1)
surreal2tex(x; level = 2)

