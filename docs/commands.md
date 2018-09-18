# Quick Doc 

This is just a list of commands, grouped in some rough fashion to list
what we can do with this package with a few examples.

Note that I started out calling the type `SurrealFinite` before
finding the better names `SurrealShort` or `SurrealDyadic`, so
everything is definite in terms of SurrealFinite, but there are
aliases. 

### Creating Surreals

+ Constructors

  + SurrealFinite(shorthand::String, L::Array{SurrealFinite}, R::Array{SurrealFinite}, h::UInt64)
  + SurrealFinite(shorthand::String, L::Array, R::Array )
  + SurrealFinite(L::Array, R::Array )
  + ≀(L::Array, R::Array)
  + zero( )
  + one( )

+ Predefines (use these in preference to manual constructions)

  + ∅ and ϕ empty array of SurrealFinite
  + SurrealZero
  + SurrealOne
  + SurrealMinusOne
  + SurrealTwo
  + SurrealMinusTwo
  + SurrealThree

+ Conversions

  + convert(::Type{SurrealFinite}, n::Int )
  + convert(::Type{SurrealFinite}, r::Rational )
  + convert(::Type{SurrealFinite}, f::Float64 ) 

+ Aliases to constructors and conversion routines

  + dali(x) = convert(SurrealFinite, x)
  + SurrealFinite(x) = convert(SurrealFinite, x)
  + SurrealShort
  + SurrealDyadic

### Working with Surreals

+ Operations (that return surreals by default)

  + +, -, * 
  + / (only for values that result in dyadics)
  + floor(), ceil(), trunc(), round()
  + mod()
  + sign()

+ Comparisons

  + === which tests if they are the same object in memory
  + == which tests identity (but ignores shorthand) 
  + ≅ which tests if values are equal (and conversely ≇)
  + <=, >=, <, > standard comparisons of surreal value
  + and array variants of these
  
+ Query

  + iszero( )
  + iscanonical( )
  + sign()
  + isinteger() 
  + isdivisible(), isodd(), iseven()
  + isnan(), isinf(), isfinite()

+ Surreal functions -- special functions unique to surreal forms

  + canonicalise
  + parents
  + ⪯ (is a parent of), ⪰, ≺, ≻
  + see also `Information functions` below

+ convert to standard types

  + convert(::Type{Rational}, s::SurrealFinite )
  + similarly for other numerical types

### I/O

+ read

  + read(io::IO, ::Type{SurrealFinite}, n::Int=1)
  + convert(::Type{SurrealFinite}, s::AbstractString ) 

+ write

  + expand -- see internal manual information
  + show(io::IO, x::SurrealFinite)
  + convert(::Type{String}, s::SurrealFinite )  
  + surreal2tex(io::IO, x::SurrealFinite; level=0)
  + pf, spf -- debugging routines, not really part of the external interface

+ graph formats

  + surreal2dag -- output to a DAG in DOT format (obsolete)
  + surreal2dot -- output to a tree in DOT format (obsolete)
  + other routines are really part of the external interface


### Informational functions

+ General
 
  + generation

+ DAG stats

  + dag_stats
  + nodes
  + edges
  + paths
  + tree_nodes
  + breadth
  + uniqueness, uniqueness_max, uniqueness_failure


### Misc

+ hash 

+ unique2!

+ caches and caching -- not sure I want to document these -- they
  aren't really part of the user interface



