using SurrealNumbers 
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

@static if VERSION < v"0.7.0"
    const devnull = DevNull
end 
 
clearcache()
 
# recursive construction and one and zero
x0 = SurrealFinite("0", ϕ,   ϕ)  
x1 = SurrealFinite("1", [x0], ϕ)
x11 = SurrealFinite("-1", ϕ,  [x0])
x111  = SurrealDyadic("-1", ϕ,  [x0])
x4 = SurrealFinite( [x11],  [x1,x0])
x5 = SurrealFinite( [x11,x0,x4],  [x1])
x21  = SurrealFinite("2", [x1],  ϕ)
x22  = SurrealFinite( [x0, x1],  ϕ)
x23  = SurrealFinite( [x11,x1],  ϕ)
x24  = SurrealFinite( [x11,x0,x1],  ϕ)
x25 = [x11,x0,x1] ≀ ϕ
z = zero(x1)
z0 = zero(SurrealFinite)
o = one(x1)
o1 = one(SurrealFinite)
x00 = SurrealFinite([x11],[x1])
z11 = x1 - x1
   
x3 = convert(SurrealFinite, 3)
x41 = convert(SurrealFinite, 4)
x42 = SurrealFinite( [x3],  ϕ)
x43 = convert(SurrealFinite, 2) * convert(SurrealFinite, 2)
x44 = SurrealFinite(4)
x6 = convert(SurrealFinite, 2) * convert(SurrealFinite, 3)

print("x1 = ")
pf(x1)
println()

println(" z = ", z)
println("x4 = ", x4)
println("x5 = ", x5)

s1 = convert(SurrealFinite, 1//2)
s2 = convert(SurrealFinite, 3//4)
s2a = convert(SurrealShort, 3//4)
s2b = convert(SurrealDyadic, 3//4) 
print("s1 = ")
pf(s1)
println()
print("s2 = ")
pf(s2)
println()

println("x41 = ", x41, ", x42 = ", x42, ", x43 = ", x43)

pf(x43)
println()
 
spf(x43)
println()

pf( dali(3) * dali(3) ) # this caused hash violations at one point, because not all equal-valued surreals were being deterministically sorted

@testset "constructor test" begin
    @test SurrealFinite("1/2", [0], [1]) ≅ 1/2
    @test_throws ErrorException SurrealFinite("1", [x1], [x0])
    @test_throws DomainError convert(SurrealFinite, NaN )
    @test_throws ErrorException SurrealFinite(1//3)
    @test_throws ErrorException SurrealFinite(100000 + 1//2)
    @test_throws ErrorException SurrealFinite(100000.000012)

    @test convert(SurrealFinite, 1.3 ) == 5854679515581645//4503599627370496
    @test all( convert(Array{SurrealFinite}, [1,2]) .== convert.(SurrealFinite, [1,2]) )
    @test z === zero(SurrealFinite)
    @test z0 === zero(SurrealFinite)
    @test o === one(SurrealFinite) 
    @test o1 === one(SurrealFinite)
    @test s2a === s2 
    @test s2b === s2
    @test x11 == x111
    @test dali(0) === SurrealZero
    @test dali(1) === SurrealOne
    @test dali(-1) == SurrealMinusOne
    @test dali(2) == SurrealTwo
    @test dali(-2) == SurrealMinusTwo
end

@testset "comparisons" begin
    @test x0 <= x1
    @test x0 < x1
    @test x0 <= x0
    @test x1 >= x11
    @test x11 < x0
    @test !( x11 ≅ x1 )
    @test x21 ≅ x22 ≅ x23 ≅ x24 # should be an equivalence class
    @test x41 ≅ x42 ≅ x43
    @test x41 == x42
    @test x44 == x41
    @test x41 != x43

    @test x0 ≺ x5
    @test x4 ≺ x5
    @test x5 ≻ x4
    @test !(x5 ≺ x4)
    @test !(x4 ≺ x4)
    @test x0 ⪯ x5
    @test x5 ⪯ x5
    @test x41 ⪰ x3
    @test s1 ⪰ s1
    @test x1 ⪯ s2
    @test x21 ≺ x43
    @test isancestor(x21, x43)

    @test all( parents(x24) .== [x11,x0,x1] )
    @test all( parents(x5) .== sort([x11,x0,x4,x1]) )
end

@testset "basic operators" begin
    @test -x1 == x11
    @test - -x0 == x0
    @test - -x4 == x4
    @test - -x4 === x4   # still flaky for version 0.7, so expect this to break
    @test x1 + x11 ≅ x0 
    @test x1 - x1 ≅ x0 
    @test x1 + x0 == x1
    @test x1*x1 ≅ x1
    @test x0*x1 ≅ x0

    @test x4*x0 ≅ x0 # this case broke things before

    @test x11/one(x11) ≅ -x1
    @test x22 / x22 ≅ one(x22)
    @test x22/x1 ≅ x22
    @test x4/x1 ≅ x4
    @test float(x1/x4) == float(x1)/float(x4)
    @test_throws ErrorException x1/x0
    @test_throws ErrorException x1/x3
    
    @test convert(SurrealFinite, 6)/ convert(SurrealFinite, 3) ≅ convert(SurrealFinite, 2)
end

@testset "convert to surreal" begin
    @test convert(SurrealFinite,2)*convert(SurrealFinite,2) ≅ convert(SurrealFinite,4)
    @test convert(SurrealFinite,2)*convert(SurrealFinite,-3) ≅ convert(SurrealFinite,-6)
    @test convert(SurrealFinite,2)*x0 ≅ x0 
    @test convert(SurrealFinite,2)*x1 ≅ convert(SurrealFinite,2)
    @test x4 ≅ -s1
    @test convert(SurrealFinite, 0.5) == convert(SurrealFinite, 1//2 )
    @test convert(SurrealFinite, 1.000000000001) == 1.000000000001
      # tests promotion
      #  float(562949953421875//562949953421312)
    
    # f =  1122342342.23422522
    # convert(SurrealFinite, 1122342342.23422522)
    @test convert(SurrealFinite, 0.5) == dali(0.5)
end

A = [x0, x1] 
B = [x11, x4]
X = [ x0, x1, x1, x11, x1, x11, x23]
Y = copy(X)
unique2!(Y) 
@testset "array operations" begin
    @test all( sort(X) .≅ convert.(SurrealFinite, [-1, -1, 0, 1, 1, 1, 2] ))
    @test all( unique(X) .≅ convert.(SurrealFinite, [0, 1, -1, 2] ))
    # @test all( Y .≅ convert.(SurrealFinite, [-1, 0, 1, 2] ))

    @test all( sort(float.(A + B)) .== [-1.0, -0.5, 0.0, 0.5])
    @test all( sort(float.(A - B)) .== [0.5, 1.0, 1.5, 2.0])
    @test A + ϕ == ϕ
    @test ϕ - B == ϕ

    @test ϕ <= A # can see from this that <= is not a proper ordering on 'sets' of surreals
    @test A <= ϕ # its just a convenience comparison operator
    @test ϕ < A
    @test A < ϕ
    @test B <= A
    @test B < A
    @test !(B > A)
    @test !(B >= A)
 
    @test all( x1*A .≅ A )
    @test all( x23*A .≅ A*x23 )

    # @test hash(A) == 0x54b75ebf9d6abfa4 # this is a bad test because the hash values are version dependent (and many OS dep as well)
end
 
V = [-6.25, -4.0, -3.25, -2.5, -2.0, -1.625, -1.0, -0.25, 0.0, 1.0, 7.0, 0.5, 1.625, 2.25, 4.5, 6.75, 8.0, 12.125]
@testset "simple functions" begin
    @test sign(x0) == x0
    @test sign(x1) == x1
    @test sign(x4) == -one(x4)
    @test sign(-x4) == one(x4)

    @test round(x22) ≅ x22 
    @test round(x23) ≅ x23
    @test round(x4) ≇ x4 
    @test round(-x4) ≇ -x4
    @test round( convert(SurrealFinite, -0.75) ) ≅ round( -0.75 )
    @test round( convert(SurrealFinite, -1.25) ) ≅ round( -1.25 )
 
    @test floor(zero(SurrealFinite)) ≅ 0.0
    @test floor(one(SurrealFinite)) ≅ 1.0
    @test floor(x22) ≅ x22 
    @test floor(x4) ≅ -one(x4)
    @test floor(x41) ≅ 4
    @test floor(x4) ≇ x4
    @test floor(-x4) ≇ -x4
    @test floor(-x4) == zero(x4)
    for v in V
        println("    testing value v = $v")
        @test floor( convert(SurrealFinite, v) ) ≅ floor( v )
        @test floor( Int, convert(SurrealFinite, v) ) == floor( v )
        @test ceil( convert(SurrealFinite, v) ) ≅ ceil( v )
        @test round( convert(SurrealFinite, v) ) ≅ round( v, RoundNearestTiesUp )
        @test trunc( convert(SurrealFinite, v) ) ≅ trunc( v )
        @test isinteger( convert(SurrealFinite, v) ) == isinteger( v )
    end
    @test_throws ErrorException floor( convert(SurrealFinite, 10000) )
    @test_throws ErrorException floor( convert(SurrealFinite, -10000) )
 
    # @test floor2(zero(SurrealFinite)) ≅ 0.0
    # @test floor2(one(SurrealFinite)) ≅ 1.0
    # @test floor2(x22) ≅ x22
     
    @test ceil(zero(SurrealFinite)) ≅ 0.0
    @test ceil(one(SurrealFinite)) ≅ 1.0
    @test ceil(x22) ≅ x22 
    @test ceil(x4) ≅ zero(x4) 
    @test ceil(x4) ≇ x4
    @test ceil(-x4) ≇ -x4
    @test ceil(-x4) == one(x4)

    @test mod(zero(SurrealFinite), x21) == 0
    @test mod(one(SurrealFinite), x21) == 1
    @test mod(x22, x21) ≅ 0.0 
    @test mod(x41, x21) ≅ zero(x4)
    @test mod(x3, x43) ≅ 3
    @test mod(-x43, x3) ≅ 2
    @test mod(s1 + 2, x21) ≅ s1
    @test_throws ErrorException mod( x3, s2 )

    @test iszero( z )
    @test !iszero( z11 )
    @test equivtozero( z )
    @test !equivtozero( x11 )

    @test isinteger(x0) == true
    @test isinteger(x22) == true
    @test isinteger(s1) == false
    @test isinteger(x43) == true
    @test isinteger(x41) == true
    @test isinteger( convert(SurrealFinite,-4) ) == true
    @test isinteger( convert(SurrealFinite,-1//4) ) == false

    @test iseven(x0) == true
    @test iseven(x22) == true
    @test iseven(s1) == false
    @test iseven(x3) == false
    @test isodd(s1) == false
    @test isodd(x3) == true   
    @test iseven( convert(SurrealFinite,-4) ) == true
    @test iseven( convert(SurrealFinite,-1) ) == false
    
    @test isinf(s1) == false
    @test isnan(s1) == false
    @test isfinite(s1)  
end


@testset "simplifications" begin
    @test canonicalise(x22) ≅ x22
    @test canonicalise(x43) ≅ convert(SurrealFinite,4)
    @test canonicalise(x43) != convert(SurrealFinite,2)*convert(SurrealFinite,2)
    @test canonicalise(x43) ≅ convert(SurrealFinite,2)*convert(SurrealFinite,2)
    @test !iscanonical(x22)
    @test iscanonical(x21)
end

g3 = SurrealFinite( [ 7//16 ], [ 15//16 ] )
g4 = SurrealFinite( [convert(SurrealFinite,3)], [convert(SurrealFinite,17)] )
# g4 ≅ convert(SurrealFinite,4)

@testset "conversion from surreal back to real" begin
    @test convert(Int64, convert(SurrealFinite, -2)) == -2
    @test convert(Int32, convert(SurrealFinite, 6)) == 6
    @test convert(Int, convert(SurrealFinite, 0)) == 0 
    # @test Int64( convert(SurrealFinite, 0) ) == 0 
    @test convert(Int64, convert(SurrealFinite, -2)) == -2
    @test convert(UInt32, convert(SurrealFinite, 6)) == 6
    @test convert(Int8, convert(SurrealFinite, 0)) == 0
    @test_throws InexactError convert(UInt8, convert(SurrealFinite, -1))
    
    @test convert(Rational, convert(SurrealFinite, 1//8)) == 1//8
    @test convert(Rational, convert(SurrealFinite, 9//8)) == 9//8
    @test convert(Rational, convert(SurrealFinite, -9//8)) == -9//8
    @test convert(Rational, convert(SurrealFinite, -27//16) ) == -27//16
    @test convert(Rational, x4) == -1//2
    @test convert(Rational, convert(SurrealFinite,2)*convert(SurrealFinite,2) ) == 4

    @test g3 ≅ convert(SurrealFinite, 0.5) 
    @test convert(Rational, SurrealFinite( [ 7//16 ], [ 15//16 ] ) ) == 1//2
    @test convert(Rational{Int32}, SurrealFinite( [ 7//16 ], [ 15//16 ] ) ) == 1//2

    @test convert(Float64, x43 ) == 4.0
    @test convert(Float64, x4) == -0.5
    @test convert(AbstractFloat, x43 ) == 4.0
    @test convert(AbstractFloat, x4) == -0.5
    @test float( x4 ) == -0.5
    @test float( s2 ) == 0.75
    @test float( g4 ) == 4.0 

    @test all(convert.(SurrealFinite, [-1, 0, 1, 2] ) == [ convert(SurrealFinite, i) for i=-1:2 ])

    @test convert(String, x1) == "1"
    @test convert(String, x4) == "-1/2"
end

@testset "promotion" begin
    # lots of this has already been implicitly tested
    # note promoted type will be in canonical form, but other side might not
    @test x0 <= 1.0
    @test x1 == 1//1
    @test 1.0 + x1 == 2.0
    @test canonicalise(1//2 + x1) == 1.5
    @test 1//2 + x1 ≅ 1.5
    @test 1//2 + x1 ≇ 2.5
end

out_dir = "Data"
a4 = convert(SurrealFinite, 13//16)
a5 = convert(SurrealFinite, -2//16)
s = "{ {{{{{|0}|{0|}}|{{0|}|}}|{{{0|}|}|}}|{{{{0|}|}|}|}} | }"    # from http://goodmath.scientopia.org/2006/08/21/arithmetic-with-surreal-numbers/
ss = convert(SurrealFinite, s)
@testset "I/O and string parsing" begin
    @test_throws Meta.ParseError convert(SurrealFinite, "{1//2| \\phi, 1}")
    @test convert(SurrealFinite, "{1//2 | 2,  1}") == SurrealFinite( [1//2], [2,1] )
    @test convert(SurrealFinite, raw"{-1//2 | {|\phi}}") ==  SurrealFinite( [-1//2], [zero(SurrealFinite)] )
    @test convert(SurrealFinite, "{-1//2 | 2.0, 3}") ==  SurrealFinite( [-1//2], [2,3] )
    
    io_test_file_str = joinpath(@__DIR__, out_dir, "test_surreals.dat")
    R = read(io_test_file_str, SurrealFinite, 4)
    @test all( float(R[1:3]) .== [0.0, -1.0, 0.0] )
    @test_throws UndefRefError R[4]

    io_test_file = joinpath(@__DIR__, out_dir, "test_surreals_io.dat")
    io = open(io_test_file, "w")
    pf(io, a4)
    pf(io, a5)
    println(io, expand(a4))
    println(io, expand(a5))
    close(io)
    YY = read(io_test_file, SurrealFinite, 2)
    @test YY[1] == a4
    @test YY[2] == a5

    io_test_file_tex = joinpath(@__DIR__, out_dir, "test_surreals_io.tex")
    io = open(io_test_file_tex, "w")
    surreal2tex(io, a4; level = 0)
    surreal2tex(io, a4; level = 1)
    surreal2tex(io, a4; level = 2)
    surreal2tex(a4; level = 2)
    close(io)
    # should compare this to a calibration file
    #  but not sure if this might introduce potential for system dependencies that aren't real errors

    @test ss ≅ 4.0
    @test ss == convert(SurrealShort, 2.0) * convert(SurrealDyadic, 2.0)
end

@testset "structured output" begin
    @test surreal2dot(devnull, x0) == 1
    @test surreal2dot(x0) == 1 
    @test surreal2dot(devnull, x23) == 5
    # @test surreal2dot(devnull, x43) == tree_nodes(x43)
    @test surreal2dag(devnull, x23) == 4
    # @test surreal2dag(devnull, x43) == nodes(x43)
    # @test surreal2dag(devnull, -x43) == nodes(x43) 
    @test surreal2dag(devnull, -s2) == nodes(s2) 

    # should compare these against a reference, but sort could get non-unique order, so potential for difference
end

a24 = dali(2) * dali(4);
(P,Q,U,V) = uniqueness(a24);   
V2 = uniqueness_failure( U, V)
U2 = filter( ((k,v),) -> v>1 , U) 
P2 = filter( ((k,v),) -> v==0xa09c3927cec6b030, P )
 
x32 = convert(SurrealFinite, 3/4) + convert(SurrealFinite, 3/4)
@testset "DAG statistics" begin
    @test dag_stats(x1) == SurrealDAGstats(2,2,1,1,ϕ,1,1,0,1,1,1)
    @test dag_stats(x41) == SurrealDAGstats(5,5,4,4,ϕ,1,4,0,4,1,1)
    @test dag_stats(x43) == SurrealDAGstats(11,21,14,6,ϕ,5,4,-1,4,2,2)
    @test dag_stats(x00) == SurrealDAGstats(4,5,4,2,ϕ,2,0,-1,1,2,2) 
    @test dag_stats(x5) == SurrealDAGstats(5,12,9,3,ϕ,6,1//2,-1,1,1,4)
    @test dag_stats(s2) == SurrealDAGstats(4,7,5,3,ϕ,3,3//4,0,1,1,2)    
    @test dag_stats(x6) == SurrealDAGstats(45,2574,82,12,ϕ,625,6,-3,7,4,3) 
    @test dag_stats(x1, Dict{SurrealFinite,SurrealDAGstats}())[2] == Dict(x0=>SurrealDAGstats(1,1,0,0,ϕ,1,0,0,0,1,0),
                                                                          x1=>SurrealDAGstats(2,2,1,1,ϕ,1,1,0,1,1,1)    )
    @test dag_stats(x43).nodes == 11
    @test dag_stats(x43).nodes == nodes(x43)
    @test dag_stats(x43).tree_nodes == 21
    @test dag_stats(x43).tree_nodes == tree_nodes(x43)
    @test dag_stats(x43).edges == 14
    @test dag_stats(x43).edges == edges(x43)
    @test dag_stats(x43).generation == 6
    @test dag_stats(x43).paths == 5
    @test dag_stats(x43).value == 4
    @test dag_stats(x43).minval == -1
    @test dag_stats(x43).maxval == 4
    @test dag_stats(x43).n_zeros == 2
    @test dag_stats(x43).max_width == 2
    @test generation(x0) == 0
    
    @test generation(x1) == 1
    @test generation(x21) == 2
    @test generation(x22) == 2
    @test generation(x4) == 2
    @test generation(x5) == 3
    @test generation( x43 ) == 6
    @test generation( convert(SurrealFinite,4) ) == 4

    @test nodes(x0) == 1
    @test nodes(x1) == 2
    @test nodes(x41) == 5
    @test nodes(x43) == 11
    @test nodes(x5) == 5
    @test nodes(s2) == 4
 
    @test breadth(x1) == 1
    @test breadth(x43) == 5

    @test paths(x43) == 5

    @test dag_stats(x32).generation == 6
    @test dag_stats(x32).generation == generation(x32)

    # @test dag_stats(x41; LP=true) == SurrealDAGstats(5, 4, 4, [0, 1, 2, 3, 4], 1,4,0,4)
    @test dag_stats(x41; LP=true).longest_path == convert.(SurrealFinite, [0, 1, 2, 3, 4])
  
    @test uniqueness_max(U) == 1 # this test causes an error if we do say 2*4, without the caching in +
end
 
d33 = dali(33//2)
@testset "memory allocation and caching" begin
    # @test Base.summarysize(  ) == 943
    # 1775 was the old size before we made "dali" reuse cached numbers
    # 1743 comes up in v0.7, instead of the number above, not sure why versions are so different,
    # but if it is implementation dependent it shouldn't be a test
    
    # depends on previous results -- could be changed by earlier operations
    @test SurrealNumbers.ExistingSurreals[ SurrealNumbers.ExistingCanonicals[-1//2] ] == -1//2
    # don't use hardcoded hash values: @test SurrealNumbers.ExistingSurreals[ SurrealNumbers.ExistingProducts[0x013b60c22ac142fa][0x84037d1bb0afff04] ] == -1
    clearcache()
    @test isempty(SurrealNumbers.ExistingSurreals)
    @test isempty(SurrealNumbers.ExistingCanonicals)
    @test isempty(SurrealNumbers.ExistingProducts)
    @test isempty(SurrealNumbers.ExistingSums)
    @test isempty(SurrealNumbers.ExistingNegations)
    Count['+'] == 0
    Count['*'] == 0 

    m = dali(2) * dali(2)
    Count['+'] == 4
    Count['+'] == 13
    Count['*'] == 6 
    Count['c'] == 20 
end


@testset "new identity test" begin
    c1 = dali(1) + dali(3)
    c2 = dali(2) + dali(2)
    @test c1 == c2
    @test c1 === c2

end

