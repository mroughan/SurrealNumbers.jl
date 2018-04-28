using SurrealNumbers 
  
x0 = SurrealFinite("0", ϕ,   ϕ)  
x1 = SurrealFinite("1", [x0], ϕ)
x11 = SurrealFinite("-1", ϕ,  [x0])
pf(x1)
println()

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
 
println(" z = ", z)
println("x4 = ", x4)
println("x5 = ", x5)


x0 <= x1
x0 <= x0
x11 >= x1
x11 ≅ x1

x21 ≅ x22 ≅ x23 ≡ x24 # should be an equivalence class
 
-x1
- -x0
-x4


x1 + x11 
x1 - x1
x1 + x0

generation(x0)
generation(x1)
generation(x21)
generation(x22)
generation(x4)
generation(x5)

x1*x1 ≅ x1
x0*x1 ≅ x0
convert(SurrealFinite,2)*convert(SurrealFinite,2) ≅ convert(SurrealFinite,4)
convert(SurrealFinite,2)*convert(SurrealFinite,-3) ≅ convert(SurrealFinite,-6)
convert(SurrealFinite,2)*x0 ≅ x0
convert(SurrealFinite,2)*x1 ≅ convert(SurrealFinite,2)

X = [ x0, x1, x1, x11, x1, x11, x23]
sort(X)
unique(X)

s1 = convert(SurrealFinite, 1//2)
s2 = convert(SurrealFinite, 3//4)
x4 ≅ -s1
pf(s1)
println()
pf(s2)
println()

pf(convert(SurrealFinite, 0.5))
println()
pf(convert(SurrealFinite, 1.000000000001))
println()


f =  1122342342.23422522
# convert(SurrealFinite, 1122342342.23422522)

simplify(x22) ≅ x22 
simplify( convert(SurrealFinite,2)*convert(SurrealFinite,2) ) ≅ convert(SurrealFinite,4)

sign(x0) == x0
sign(x1) == x1

generation( convert(SurrealFinite,2)*convert(SurrealFinite,2) ) == 6
generation( convert(SurrealFinite,4) ) == 4

round(x22) ≅ x22 
round(x23) ≅ x23
round(x4) ≇ x4
round(-x4) ≇ -x4 

convert(Rational, convert(SurrealFinite, 1//8))
convert(Rational, convert(SurrealFinite, 9//8))
convert(Rational, convert(SurrealFinite, -9//8))
convert(Rational, convert(SurrealFinite, -27//16) )
convert(Rational, x4)
convert(Rational, convert(SurrealFinite,2)*convert(SurrealFinite,2) )

convert(AbstractFloat, convert(SurrealFinite,2)*convert(SurrealFinite,2) )
convert(AbstractFloat, x4)
float( x4 )
float( s2 )

isinteger(x0) == true
isinteger(x22) == true
isinteger(s1) == false
isinteger( convert(SurrealFinite,2)*convert(SurrealFinite,2) ) == true


g4 = SurrealFinite( [convert(SurrealFinite,3)], [convert(SurrealFinite,17)] )
g4 ≅ convert(SurrealFinite,4)
float( g4 )


simplify(convert(SurrealFinite,2)*convert(SurrealFinite,2) ) != convert(SurrealFinite,2)*convert(SurrealFinite,2)
simplify(convert(SurrealFinite,2)*convert(SurrealFinite,2) ) ≅ convert(SurrealFinite,2)*convert(SurrealFinite,2)

# test promotion
1.0 + x1 == 2.0


SurrealFinite( [1,2], [3,4] )

SurrealFinite( [ 7//16 ], [ 15//16 ] ) ≅ convert(SurrealFinite, 0.5)
convert( Rational, SurrealFinite( [ 7//16 ], [ 15//16 ] ) )


