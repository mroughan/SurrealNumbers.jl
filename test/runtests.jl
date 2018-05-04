using SurrealNumbers 
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

# recursive construction and one and zero
x0 = SurrealFinite("0", ϕ,   ϕ)  
x1 = SurrealFinite("1", [x0], ϕ)
x11 = SurrealFinite("-1", ϕ,  [x0])
x4 = SurrealFinite( [x11],  [x1,x0])
x5 = SurrealFinite( [x11,x0,x4],  [x1])
x21  = SurrealFinite("2", [x1],  ϕ)
x22  = SurrealFinite( [x0, x1],  ϕ)
x23  = SurrealFinite( [x11,x1],  ϕ)
x24  = SurrealFinite( [x11,x0,x1],  ϕ)
x25 = [x11,x0,x1] ≀ ϕ
z = zero(x1)
z = zero(SurrealFinite)
z = one(x1)

x3 = convert(SurrealFinite, 3)
x41 = convert(SurrealFinite, 4)
x42 = SurrealFinite( [x3],  ϕ)
x43 = convert(SurrealFinite, 2) * convert(SurrealFinite, 2)

print("x1 = ")
pf(x1)
println()

println(" z = ", z)
println("x4 = ", x4)
println("x5 = ", x5)

s1 = convert(SurrealFinite, 1//2)
s2 = convert(SurrealFinite, 3//4)
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

@testset "constructor test" begin
    @test SurrealFinite("1/2", [0], [1]) ≅ 1/2
    @test_throws ErrorException SurrealFinite("1", [x1], [x0])
    @test all( convert(Array{SurrealFinite}, [1,2]) .== convert.(SurrealFinite, [1,2]) )
    
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
    @test x41 != x43
end
 
@testset "basic operators" begin
    @test -x1 == x11
    @test - -x0 == x0
    @test - -x4 == x4 
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


@testset "surreal specific functions" begin
    @test generation(x0) == 0
    @test generation(x1) == 1
    @test generation(x21) == 2
    @test generation(x22) == 2
    @test generation(x4) == 2
    @test generation(x5) == 3
    @test generation( convert(SurrealFinite,2)*convert(SurrealFinite,2) ) == 6
    @test generation( convert(SurrealFinite,4) ) == 4

    @test size(x0) == 1
    @test size(x1) == 2
    @test size(x41) == 5
    @test size(x43) == 21
    @test size(x5) == 12
    @test size(s2) == 7

    @test n_zeros(x0) == 1
    @test n_zeros(x1) == 1
    @test n_zeros(x41) == 1
    @test n_zeros(x43) == 5
    @test n_zeros(x5) == 6
    @test n_zeros(s2) == 3

    @test depth_max(x0) == 0
    @test depth_max(x1) == 1
    @test depth_max(x41) == 4
    @test depth_max(x43) == 6
    @test depth_max(x5) == 3
    @test depth_max(s2) == 3

    @test depth_av(x0) == 0.0
    @test depth_av(x1) == 1.0
    @test depth_av(x41) == 4.0
    @test depth_av(x43) == 6.0
    @test depth_av(x5) == 13/6
    @test depth_av(s2) == 7/3

    @test all( sort(count_n(x0)) .== [0] )
    @test all( sort(count_n(x1)) .== [0, 1] )
    @test all( sort(count_n(x41)) .== [0,1,2,3,4] )
    @test all( sort(count_n(x43)) .== [-1,0,0,0,0,0,0,1,1,1,1,1,2,2,2,2,3,3,3,4,4] )
    @test all( sort(count_n(x5)) .== [-1, -1, -1//2, 0, 0, 0, 0, 0, 0, 1//2, 1, 1] )
    @test all( sort(count_n(s2)) .== [0, 0, 0, 1//2, 3//4, 1, 1] )
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

    @test B <= A
    @test B < A
    @test !(B > A)
    @test !(B >= A)

    @test all( x1*A .≅ A )
    @test all( x23*A .≅ A*x23 )
end

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
    @test floor(x4) ≇ x4
    @test floor(-x4) ≇ -x4
    @test floor(-x4) == zero(x4)

    @test ceil(zero(SurrealFinite)) ≅ 0.0
    @test ceil(one(SurrealFinite)) ≅ 1.0
    @test ceil(x22) ≅ x22 
    @test ceil(x4) ≅ zero(x4)
    @test ceil(x4) ≇ x4
    @test ceil(-x4) ≇ -x4
    @test ceil(-x4) == one(x4)

     # add some for ceil as well
    
    @test isinteger(x0) == true
    @test isinteger(x22) == true
    @test isinteger(s1) == false
    @test isinteger( convert(SurrealFinite,2)*convert(SurrealFinite,2) ) == true
    @test isinteger( convert(SurrealFinite,-4) ) == true
    @test isinteger( convert(SurrealFinite,-1//4) ) == false

    @test isinf(s1) == false
    @test isnan(s1) == false
    @test isfinite(s1)  
end


@testset "simplifications" begin
    @test canonicalise(x22) ≅ x22 
    @test canonicalise( convert(SurrealFinite,2)*convert(SurrealFinite,2) ) ≅ convert(SurrealFinite,4)
    @test canonicalise( convert(SurrealFinite,2)*convert(SurrealFinite,2) ) != convert(SurrealFinite,2)*convert(SurrealFinite,2)
    @test canonicalise( convert(SurrealFinite,2)*convert(SurrealFinite,2) ) ≅ convert(SurrealFinite,2)*convert(SurrealFinite,2)
end

g3 = SurrealFinite( [ 7//16 ], [ 15//16 ] )
g4 = SurrealFinite( [convert(SurrealFinite,3)], [convert(SurrealFinite,17)] )
# g4 ≅ convert(SurrealFinite,4)

@testset "conversion from surreal back to real" begin
    @test convert(Integer, convert(SurrealFinite, -2)) == -2
    @test convert(Integer, convert(SurrealFinite, 6)) == 6
    @test convert(Integer, convert(SurrealFinite, 0)) == 0
    
    @test convert(Rational, convert(SurrealFinite, 1//8)) == 1//8
    @test convert(Rational, convert(SurrealFinite, 9//8)) == 9//8
    @test convert(Rational, convert(SurrealFinite, -9//8)) == -9//8
    @test convert(Rational, convert(SurrealFinite, -27//16) ) == -27//16
    @test convert(Rational, x4) == -1//2
    @test convert(Rational, convert(SurrealFinite,2)*convert(SurrealFinite,2) ) == 4

    @test g3 ≅ convert(SurrealFinite, 0.5) 
    @test convert( Rational, SurrealFinite( [ 7//16 ], [ 15//16 ] ) ) == 1//2

    @test convert(AbstractFloat, convert(SurrealFinite,2)*convert(SurrealFinite,2) ) == 4.0
    @test convert(AbstractFloat, x4) == -0.5
    @test float( x4 ) == -0.5
    @test float( s2 ) == 0.75
    @test float( g4 ) == 4.0 

    @test all(convert.(SurrealFinite, [-1, 0, 1, 2] ) == [ convert(SurrealFinite, i) for i=-1:2 ])

    @test convert(String, x1) == "1"
    @test convert(String, x4) == "-1/2"
end

@testset "promotion" begin
    # note promoted type will be in canonical form, but other side might not
    @test x0 <= 1.0
    @test x1 == 1//1
    @test 1.0 + x1 == 2.0
    @test canonicalise(1//2 + x1) == 1.5
    @test 1//2 + x1 ≅ 1.5
    @test 1//2 + x1 ≇ 2.5
end

@testset "output" begin
    @test surreal2dot(DevNull, x0) == 1
    @test surreal2dot(DevNull, x23) == 5
    @test surreal2dot(DevNull, x43) == size(x43)
end
