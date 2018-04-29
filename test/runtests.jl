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


@testset "comparisons" begin
    @test x0 <= x1
    @test x0 <= x0
    @test x1 >= x11
    @test !( x11 ≅ x1 )
    @test x21 ≅ x22 ≅ x23 ≅ x24 # should be an equivalence class
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
end
    
@testset "generation function" begin
    @test generation(x0) == 0
    @test generation(x1) == 1
    @test generation(x21) == 2
    @test generation(x22) == 2
    @test generation(x4) == 2
    @test generation(x5) == 3
    @test generation( convert(SurrealFinite,2)*convert(SurrealFinite,2) ) == 6
    @test generation( convert(SurrealFinite,4) ) == 4
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

X = [ x0, x1, x1, x11, x1, x11, x23]
Y = copy(X)
unique2!(Y) 
@testset "array operations" begin
    @test all( sort(X) .≅ convert.(SurrealFinite, [-1, -1, 0, 1, 1, 1, 2] ))
    @test all( unique(X) .≅ convert.(SurrealFinite, [0, 1, -1, 2] ))
    # @test all( Y .≅ convert.(SurrealFinite, [-1, 0, 1, 2] ))
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

    @test floor(x22) ≅ x22 
    @test floor(x4) ≅ -one(x4)
    @test floor(x4) ≇ x4
    @test floor(-x4) ≇ -x4
    @test floor(-x4) == zero(x4)

    # add some for ceil as well
    
    @test isinteger(x0) == true
    @test isinteger(x22) == true
    @test isinteger(s1) == false
    @test isinteger( convert(SurrealFinite,2)*convert(SurrealFinite,2) ) == true
end

@testset "simplifications" begin
    @test simplify(x22) ≅ x22 
    @test simplify( convert(SurrealFinite,2)*convert(SurrealFinite,2) ) ≅ convert(SurrealFinite,4)
    @test simplify(convert(SurrealFinite,2)*convert(SurrealFinite,2) ) != convert(SurrealFinite,2)*convert(SurrealFinite,2)
    @test simplify(convert(SurrealFinite,2)*convert(SurrealFinite,2) ) ≅ convert(SurrealFinite,2)*convert(SurrealFinite,2)
end

g3 = SurrealFinite( [ 7//16 ], [ 15//16 ] )
g4 = SurrealFinite( [convert(SurrealFinite,3)], [convert(SurrealFinite,17)] )
# g4 ≅ convert(SurrealFinite,4)

@testset "conversion from surreal back to ordinary numbers" begin
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

end


@testset "promotion" begin
    @test 1.0 + x1 == 2.0
end



