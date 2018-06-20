# SurrealNumbers

[![Build Status](https://travis-ci.org/mroughan/SurrealNumbers.jl.svg?branch=master)](https://travis-ci.org/mroughan/SurrealNumbers.jl)

[![Coverage Status](https://coveralls.io/repos/mroughan/SurrealNumbers.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/mroughan/SurrealNumbers.jl?branch=master)

[![codecov.io](http://codecov.io/github/mroughan/SurrealNumbers.jl/coverage.svg?branch=master)](http://codecov.io/github/mroughan/SurrealNumbers.jl?branch=master)


## Intro 

This is a slightly crazy package implementing some parts of the
Surreal Number system invented by John Horton Conway, and explained by
Knuth in ["Surreal Numbers: How Two Ex-Students Turned on to Pure
Mathematics and Found Total Happiness."](http://www.math.harvard.edu/~knill/teaching/mathe320_2015_fall/blog15/surreal1.pdf)

It isn't intended to be useful, so much as educational. It was
educational for me to write, in terms of learning Julia. I wanted a
task that would be painful to do in Matlab but "easy" to do in
Julia. It might also be helpful, I hope, for someone trying to learn
about Surreal numbers. I certainly did learn a lot about them that
would have been easily brushed to the side if I had only gone through
the theorems.

## Background: the Surreal Numbers

Surreal numbers aren't numbers as we are taught, but they have many of
the same properties. The tricky thing is that they are defined
recursively from the very start. The nice part is that they use only
set operations.

The definition is as follows: a surreal number $x$ is an ordered pair of
sets of surreal numbers (call them $X_L$ and $X_R$) such that every
member of the left set is `<` all of the members of the right set.

There is a starting point -- we can always use empty sets -- as so the
first surreal number (usually denoted zero, because it will turn out
to be the additive identity) is $\{ | \}$, where I will use this
 bracket and pipe notation to denote $x = \{ X_L | X_R\}$. Also, to
 make the empty spaces clearer (and the coding more efficient) I have
 defined $ ϕ = \{ \} $ the empty set. 

Then on the "first day" a new generation of surreals can be created in
terms of the initial case. On the second day we create a new generation
and so on. Each has a meaning corresponding to traditional numbers in
order to place a consistent interpretation on standard mathematical
operators defined on the surreals. 

Just to reiterate, the tricky thing is that everything is
recursive. Even comparitors like `<`, and hence, we can't even know
if something is a valid surreal in its own right, but only through
recursively investigating its component sets. 

In any case, this is not a survey or tutorial on surreal
numbers. There are many out there better than I can write, e.g.,

+ https://en.wikipedia.org/wiki/Surreal_number
+ https://conservancy.umn.edu/bitstream/handle/11299/174778/Hebert_umn_0130M_15912.pdf?sequence=1
+ https://math.stackexchange.com/questions/816540/proof-of-conways-simplicity-rule-for-surreal-numbers
+ http://www.math.harvard.edu/~knill/teaching/mathe320_2015_fall/blog15/surreal1.pdf
+ https://www.cut-the-knot.org/WhatIs/Infinity/SurrealNumbers.shtml
+ https://www.whitman.edu/Documents/Academics/Mathematics/Grimm.pdf
+ https://www.tondering.dk/download/sur16.pdf

The purpose here is simply to introduce the key reasons for
implementing this in Julia -- it enables users to define
high-performance types of their own, and those types can be
recursively defined.

I know you can do this in other languages, and in some cases also
achieve high performance. But it's *so* easy in Julia. 

On the other hand, while surreals use sets, and Julia has a Set type,
I didn't find that easy to work with. As a type, it doesn't seem to
have all the bells and whistles I would expect from sets.

I could have added these, but Julia creates a whole suite of array
functions automagically when you define scalars so using arrays was a
low pain way to get things working. That is, although the surreals use
sets, i.e., order of the elements is not important, almost all texts
do write these sets in order. So I felt justified in putting these
elements into sorted arrays. It works nicely, as long as you make sure
to only put unique entries into the array, and this is a little more
tricky for surreals for reasons we will get to in a moment.

Someone may correct me about the best way to do this, but for the
moment, using arrays seemed like a low pain way to get surreals working. 

### An example

Let's start off with some small examples of using the package. After
installing the package we can create basic surreals using (i) the
constructor, (ii) conversion from another real number, or (iii) a couple of
special functions, e.g., `zero` and `one`. There are two constructors,
one includes an extra string we'll call the *shorthand* for the
surreal. It's use in printing out numbers. The second constructor, and
many other operators leave this blank. The empty set is designated by
ϕ, which you can get in Julia by typing `\phi TAB`. 

    julia> using SurrealNumbers
    julia> z = zero(SurrealFinite)
    julia> x1 = SurrealFinite("1", [z], ϕ)
    julia> x2 = SurrealFinite("2", [x1], ϕ)
    julia> x_something = SurrealFinite([z, x1], ϕ)
    julia> x_half = convert(SurrealFinite, 1 // 2)

These commands create several surreals, starting at zero, then 1, and
2, showing the recursive construction. Then a surreal whose value we
might not know (to start with), and then we convert the rational value
1/2 into a surreal.

Note that the output of these varies: in the case where a shorthand
was defined it just outputs (in bold) the shorthand, but otherwise it
will show the bracket format of the components. To see the bracket
format even when there is a shorthand defined use the command `pf`, e.g.

    julia> println(" x_half = ", x_half, " = ")
    julia> pf(x_half)
    julia> println()

We have most of the standard arithmetic operators (division has some
restrictions -- see below), so you can do things such as

    julia> x2 + x_half 

which will produce quite a long sequence. To see what it is, convert
back to a real number, 

    julia> float( x2 + x_half )

or do a picture of the recursion using DOT (which needs to be
separately installed from [GraphViz](https://www.graphviz.org/)), e.g.,

    julia> file = "test_dot.dot"
    julia> FID = open(file, "w")
    julia> surreal2dot(FID, x2 + x_half)
    julia> close(FID)
    julia> run(`dot -Tsvg -O $file`) 

with the following result. Note that the red notations were added
manually. Each box is a surreal number, designated by the number at
the top of the box, and its left and right sets are in the
corresponding boxes below. The recursion for each is show below, down
to the point where each recursion stops at zero.
 
![3/2](/test/Data/test_dot.svg)

That seems like enough to get you started, so now a little about the
implementation. 

## Icky bits -- implementation details

The surreals include all real numbers (and infinity and epsilon and
others). However, many of these require infinitely large sets $X_L$
and/or $X_R$. I have some ideas about how to do that (using lazy
evaluation), but they aren't fleshed out yet so for the moment, I will
restrict myself to finite surreals, i.e., surreals with finite sets
$X_L$ and $X_R$. 

So the type `Surreal` is an abstract type (a subtype of `Real`) with
at present only one useful subtype `SurrealFinite`, where finite
(here) means that the representation is finite, not that the number
itself is finite. 

This is a pretty rich set by itself, but it doesn't cover even the
entire set of Rationals; only the *dyadic* rationals. So a few words
on them seem in order (see below).

The plan is to add a `SurrealTrans` for the transfinites and other
surreals with infinite representations. But that is a little harder to
do, so I leave that for the moment.

The actual type structure (minus the constructor -- see the code) is
just

    struct SurrealFinite <: Surreal
        shorthand::String
        L::Array{SurrealFinite,1} 
        R::Array{SurrealFinite,1} 
    end 

Note the addition of a `shorthand` string, which isn't necessary, but
carries a little bit of extra information to make pretty printing a
little easier.


### Dyadic numbers

The dyadic rational numbers are those that have a denominator that is
an exact power of 2, that is, numbers of the form (note seems that
GitHub doesn't let me use Mathjax to interpret maths). 

     \[ x = \frac{ n }{ 2^k } \]

It turns out that every dyadic has a finite representation as a
surreal, and every finite surreal represents a dyadic. The easiest way
to represent this is as a DAG (Directed Acyclic Graph) of the constructions as shown below:

![dyadics](/examples/Data/dyadic_tree_3.dot.svg)
 
It might seem a little limiting to be restricted to this set, but
remember that floating point numbers are dyadics. They are a (binary)
integer (the *significand*) multiplied by 2 to the power of a (binary)
exponent (just called the *exponent*). Thus all (finite) floating
point numbers have a finite surreal representation (though it may be
very, very long).

I did start writing a type for dyadics, but I'm not sure what it would
be useful for, so it didn't get far (I realise the irony here -- the
surreals themselves aren't exactly useful).

### Converting a number to a surreal

The description of a surreal given above generates a *form*. The forms
satisfy the rules of arithmetic ($+$, $-$, $\times$ and $/$), and so
we can identify these with the real numbers. However, there are many
*forms* that equate to the same real number. Thus there are sets of
*forms* that are equivalent in the sense that $x \equiv y$ if and only
if $x \leq y$ and $y \leq x$.

I think of this loosely by analogy to the rationals, e.g., we can have

    2 // 4 == 1 // 2

However, it seems important to distinguish equality from equivalence
here. In programming terms two "things" are equal when they are the
same, not some airy-fairy notion of equivalence, so equality and
equivalence have different meanings and uses. Hence, here we have the
relation ≅ defined to test equivalence separate from ==. 

BTW, here we hit one of the weirdnesses of Julia; 99\% of the time,
you can redefine operators and comparators to do whatever you like on
your new type. But you can't redefine $===$ or $\equiv$. The blog I
read suggests that this is because this is a core operation, that
might cause problems if a user broke it.

 + https://docs.julialang.org/en/stable/devdocs/functions/#Builtins-1
 + https://discourse.julialang.org/t/overload-for-custom-type/4898

Anyway, I have defined *congruence* $\cong$ or `≅` to do the same
thing, check for equivalence. However, as there are many possible
surreal forms we could create to represent any given real number, we
have to chose one. Call that the *canonical* form. We could define it
in several ways, but the standard is

+ zero => $0 = \{ | \}$

+ integers n => $n+1 = \{ n | \}$

+ dyadic fraction => $ x = \frac{ n }{ 2^k } $ (for n odd) becomes

     \[  \{ x - 1/2^{k} | x + 1/2^{k} \} \]

+ negative number => use the identity that $-x = \{-X_R | -X_L\}$

This is the set of conversions implemented, and it includes all
dyadics, and hence floats. Note the recursion inherent in the
construction/conversion. So some care (particularly) with floats
should be taken not to create a variable that exhausts the stack. The
current code isn't very clever in checking for this another place some
work is needed. 

Examples:

    julia> pf( convert(SurrealFinite, 1//2) )
    { 0 | 1 }
    julia> pf( convert(SurrealFinite, 3//4) )
    { 1//2 | 1 } 
    julia> pf( convert(SurrealFinite, -11//8) )
    { -3//2 | -5//4 }

The `pf` function used here is a "print-in-full", which prints the
left and right sets of the surreal, not the real equivalent
shorthand. 

### Converting a surreal back to a real

The conversion of canonical forms is relatively easy, but remember
that non-canonical forms are possible, and can be quite
counter-intuitive. 

For instance, naively, you might expect that the form $\{ 3 |
17 \}$ could be mapped to the mean of the two sets, i.e., 10.
However, $x = \{ X_L | X_R \}$ is the simplest number such that $X_L < x
< X_R$, so, in fact, this form is equivalent to 4.

Note, often in texts it is written $X \not \leq Y$ whereas I am
writing $X > Y$. The original definition is intended (I think) to take
care of the case where one or the other is the empty set, but we can
equally define $>$ to be synonomous with $\not \leq$, and just move
on. 

Conway defines a surreal $x = \{ X_L | X_R \}$ to be the simplest
$x$ that satisfies $X_L < x < X_R$, where simplest means comes from
the earliest generation.

The secret of the conversion is again to use recursion, but that is
not quite enough in this case. We use several tricks along the way:

+ If $x$ is equivalent to a known surreal such as 0 or 1, we convert
directly
+ If $x$ is negative, we use the identity that $-x = \{ -X_R | -X_L \}$
+ And, most importantly, if $0 < x < 1$ we know $x$ will be the
"simplest" dyadic rational number such that $x_L < x < x_R$. In the
interval $(0,1)$ simplest means having the denominator with the lowest
power, i.e., in order of preference we would like the denomator to the
$1,2,3, \ldots$.

We can find the latter though a simple modification of the standard
binary search a simplified version of which is shown below.
 
     a = 0; b = 1
     while true
       d = (a + b) / 2   
       if xl < c < xr
         return d
       elseif c <= xl
         a = d 
       elseif c >= xr 
         b = d
       end
     end
  
The first and last rules allow us to convert any number $x \in
[0,1]$. To convert numbers above this range, i.e., $x > 1$, we
subtract 1 (the surreal additive identity), convert the result
(recursively), and then add back 1 (this time a real). To convert
negative numbers we apply the second rule.

The result is not the world's most beautiful code -- I'm sure it can
be improved, but there are so many other inefficiency's here, I am not
sure it warrants it.


### Implementations of operators and standard functions

Most of the operators follow standard surreal definitions and defining
them in Julia is easy. They are all recursive, as you might guess,
and so very inefficient -- I wouldn't want to do demanding
computations this way, but they are easy to understand, for instance

     +(x::SurrealFinite, y::SurrealFinite) = SurrealFinite([x.X_L .+ y; x .+ y.X_L], [x.X_R .+ y; x .+ y.X_R] )

Notice that we are exploiting here Julia's natural extension of
operators from scalar to vectors (this is one of the reasons that
using arrays instead of sets for $X_L$ and $X_R$ is appealling). Thus we can write

     x .+ [y_1, y_2]

once addition is defined on the scalar surreal type, without any additional definitions, and this is particular appealling here as the scalar addition operator is recursively defined in terms of the vector+scalar addition `.+`. 

Mutiplication was a bit of a bugbear to get right because I
misinterpreted the definition. The definition of multiplication of
$x= \{ X_L |X_R \}$ and $y=\{ Y_L | Y_R \}$ had terms like

    X_L y + x Y_L - X_L Y_L

I assumed the `+` operator in this definition was the standard
(surreal) addition of sets of surreals, but it isn't. Instead the
entire expression should be interpreted as "for all pairs of elements
from X_L and Y_L perform the above computation, and form a set from
the results." That's easy enough to do once you work out what you are
trying to do. There is another underlying problem which is that
canonical forms all have sets with 1 or 0 elements, and my broken
multiplication worked for them. It was only when I created products of
non-canonical forms that I saw problems, and then only when I had
cleaned out another misunderstanding from the code.

The one interesting thing to note is that even simple computations often
generate non-canonical forms. Part of the aim of this package was to
let people experiment and see such things. It is otherwise far to
laborious to calculate anything but the very simplest cases (none of
the above texts do any but the simplest).

One really useful example (for me to understand the surreals) is the
following, which was obtained by multiplying the canonical form of 2
with itself. You might naively think that operators applied to
canonical forms resulted in canonical forms, but this is not the
case. The form is quite complex. Moreover, superficially it recurse
into itself, i.e., the surreal number 4 appears twice in the
tree. However, the form is defined in terms of the tree, so the
dependent "4" is not the same surreal as the parent "4". That is, they
are not  equal, but they are equivalent. Thus the parent "4" is
defined in terms of simpler forms (even though one of these is
equivalent). 

![2x2 = 4](/test/Data/test_dot_x43.svg)

The other pieces of the toolkit are the standard things you expect to
be able to do with numbers, e.g. round, sign, isinteger, ...  I
haven't implemented all of these, but hopefully enough that any others
can be easily added.

There are two approaches: one is to use intrinsic surreal arithmetic,
e.g. `sign` and `abs` are implemented using native surreal arithmetic and operators. The result is that they looks almost exactly like it would for any other number.

    sign(x::SurrealFinite) = x<zero(x) ? -one(x) : x>zero(x) ? one(x) : zero(x)
    abs(x::SurrealFinite) = x<zero(x) ? -x : x

Actually, I don't even have to define `abs` as I get this for free
because Julia has a similar operation defined `abs(x::Real)` for all
real numbers. This is a great feature of Julia's type hierarchy and
multiple dispatch function selection. 

The other approach to function definition is somewhat of a cheat. It
involves converting the number to a real, and then using the
appropriate operation on that field. I have tried to avoid that
approach when possible. At the moment, only the `simplify` function
uses this approach. But if I am going to do more complex math
functions, e.g., logs or trig functions, I think I will have to take
this approach. The mathematical definitions of such in pure surreal
terms are obscure. 

The cheat is used as part of one of the pieces of this that is hard to
implement, namely division. Division is well-defined on the surreals,
but, as even simple divisions such as 1/3 result in non-dyadic rational
numbers, and hence have a infinite representations as a surreal, the
class of "finite" surreals defined here is not *closed* under
division. Eventually, this can be solved by outputting a non-finite
surreal, but for the moment I have only implemented division when
dividing by $2^k$ or by a divisor of the numerator of the dyadic. The
first is easy (it can be implemented in terms of multiplication). The
second (at the moment) requires the cheat (see below) of converting
back into a rational number, but in either case the check to see which
case requires the conversion. Some work is needed here.

The `show` command is designed to show the full structure unless there
is a "shorthand" defined for a surreal (most of the simple conversions
will set this up). This aids in viewing the surreals succintly, but
sometimes we want to see deeper. In the case where shorthand is
defined we can use the command `pf` to see deeper, but it will stop at
the first layer below with a shorthand.

To see the full recursion I have implemented output of a surreal form
into the DOT syntax from [GraphViz](https://www.graphviz.org/). The
function `surreal2dot` can output a `.dot` file, and then this can be
parsed (assuming you have GraphVis installed) by commands such as 

    using SurrealNumbers

    s2 = convert(SurrealFinite, 3//4)
    file = "test_dot_s2.dot"
    FID = open(file, "w")
    surreal2dot(FID, s2)
    close(FID)
    run(`dot -Tsvg -O $file`)

    x5 = SurrealFinite( convert.(SurrealFinite, [-1, -1//2, 0]), [one(SurrealFinite)] )
    file = "test_dot_x5.dot"
    FID = open(file, "w")
    surreal2dot(FID, x5)
    close(FID)
    run(`dot -Tsvg -O $file`)

Which produces the figures like those below, illustrating the
recursive definition of the two surreal numbers given. 

![s2](/test/Data/test_dot_s2.dot.svg) <span>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span>
![x5](/test/Data/test_dot_x5.dot.svg). 

The thing to note about these tree representations is that they are inefficient,
the same surreal forms are repeated. The code includes a function `surreal2dag`
that you use in exactly the same way to generate a DAG (Directed Acyclic Graph)
representation of the surreal as you see below.

![s2](/test/Data/test_dag_s2.dot.svg) <span>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span>
![x5](/test/Data/test_dag_x5.dot.svg). 


See the code `test/test_dot.jl` for other examples. 

###  Uniquely surreal functions

There are two pieces that are unique to surreals:

**Generation:** the generation of a surreal is 1 + the generation of
the surreals used to construct it. Again this is easy to implement
recursively.  Generation comes from Knuth's story where it's called
the "birth day" of the surreal. Generation can be thought of as a
function $\rho(\cdot)$

     \[ \rho(x) = \sup_{i \in X_L \cup X_R} \big[ \rho(i) + 1 \big] \]

Implicitly, generation is the tree *depth*, if viewed as the recursive
tree shown above. 

Generation is important as $x = \{ X_L | X_R \}$ is defined to be the
simplest number such that $X_L < x < X_R$, i.e., the number with the
lowest birth day. 

Note that within an equivalence class of forms, not all forms have the
same generation. 

**Canonicalise:** convert a surreal into its equivalent canonical
form. The easiest way to implement this was to use a similar cheat to
that above, i.e., convert to a real, and then convert back to the
equivalent surreal in canonical form.

### Arrays of Surreals

We mentioned arrays of surreals earlier -- they are used in the
implementation instead of sets. One of the nice things about Julia us
that you get many of the array operators and functions for free when
you create scalar operators. So, for instance, you can write

    convert.(SurrealFinite, [-1, 0, 1, 2] )

which will broadcasts the convert function across the array of
integers to create an array containing the respective
surreals. Likewise, we can use *comprehensions* to construct arrays,
e.g.,

    [ convert(SurrealFinite, i) for i=-1:2 ]

Or we could construct and iterator for the same thing (once `floor`
and some promotion rules are defined), e.g.., the iterator from -1 to 2 is  

    convert(SurrealFinite, -1 ):convert(SurrealFinite, 2)

However, in order to use arrays as sets we need, in the constructor
for a surreal to reduce the "set" to a sorted array containing unique
elements. Julia has nice sort and unique functions, but they rely on
hash functions, so we have to add a hash. These need to be recursive,
and work for arrays of surreals as well. The hash function help says that
we should implement such for new types such that `isequal(x,y)`
implies `hash(x)==hash(y)`, with a second argument to be mixed in the
results. This is linked to the idea that we should separate the == and
the *congruent* comparisons -- my hash is based on two forms being 
equal, not being equivalent. In any case, Julia's syntax is again
simple and concise for specifying the hashes need.  The approach I
adoped was the following (which may well not be the best, but seems to
work).

    hash(x::SurrealFinite, h::UInt) = hash( hash(x.X_L, h) * hash(x.X_R, h), h )
    hash(X::Array{SurrealFinite}, h::UInt) = isempty(X) ? hash(0) : hash( hash(X[1]) * hash( X[2:end]), h )

Having set this up, I realised that when specifying the component sets
of a surreal, I actually did want my unique function to only include
one element from each equivalence class, i.e., I could have used
congruences as equality, but in general these are short lists (with
lots of recursion internally, but the arrays themselves are short), so
I just wrote my own unique function, `unique2!`. Yes I know the name
isn't great, but for newbies note that the ! is not a `not` symbols
here. In Julia's idiom this signifies that the function modifies its
arguments.

https://docs.julialang.org/en/stable/manual/style-guide/

In any case, `unique2!` also uses sort, but then crudely eliminates
duplicates based on congruence. It is then called as part of the
constructor for a surreal, which also checks the condition that $X_L <
X_R$, i.e., no element of $X_R$ is $\leq$ an element of $X_L$.

## Other comments

Don't bother to tell me that this is horribly inefficient. It will
never be anything but. Surreals were not created with numerical
computing in mind. They are about as inefficient a way to do
calculations as I can reasonably think of (excluding the addition of
nullops everywhere).

And note that this isn't really exploring an entire chunk of the
surreals, i.e., the transfinite numbers that can be represented this
way. I'll get to that one day, time gods willing.

There are other implementations of the surreals. For instance
 
+ Coq https://dl.acm.org/citation.cfm?id=2150520
+ Coq https://github.com/pirapira/surreal
+ Haskell https://github.com/Lacaranian/surreal (only integers)
+ Haskell https://github.com/serialhex/Surreal-Numbers (not quite working)
+ Haskell https://github.com/elfakyn/Haskell-surreals (only integers)
+ Mathematica http://demonstrations.wolfram.com/GeneratingTheSurrealNumbers/
              (only generation)
+ Python https://github.com/codeinthehole/python-surreal (not much
         implemented)
+ Python https://github.com/314eter/surreal-numbers
+ Ruby http://raganwald.com/2009/03/07/elegance-and-the-surreals.html


And some of these languages might be more appropriate in some ways for
this task. But I wanted to learn Julia, and see how far I could take
it here. Moreover, most of these are at least as incomplete as the
code here. For instance, none (as far as I am aware) implement
division (the Julia toolkit here has a very limited version of
division included). 


### More information about Surreals

+ https://www.ics.uci.edu/~eppstein/cgt/surreal.html


### More Examples

Some examples of code are included here to make a little of this more real. 

#### Creating surreals

Let's start by creating some surreals

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
    s1 = convert(SurrealFinite, 1//2)
    s2 = convert(SurrealFinite, 3//4)

You can check the values using `float` to convert them back to
standard floating point numbers. Note that ϕ is shorthand for an empty array of surreal numbers, which is quite helpful in many places.

Try printing them out in various forms:

    print("s1 = ", s1, " = " )
    pf(s1)
    println()


#### Operations on surreals and arrays of surreals

Test out some operations on the above variables with expected results

     # comparisons
     x0 <= x1
     x0 <= x0
     x1 >= x11
     !( x11 ≅ x1 )
     x21 ≅ x22 ≅ x23 ≅ x24

     # addition and subtraction
     -x1 == x11
     - -x0 == x0
     - -x4 == x4 
     x1 + x11 ≅ x0
     x1 - x1 ≅ x0
     x1 + x0 == x1

     # multiplication
     x1*x1 ≅ x1
     x0*x1 ≅ x0
     x4*x0 ≅ x0 
     convert(SurrealFinite,2)*convert(SurrealFinite,2) ≅ convert(SurrealFinite,4)
     
     # division
     x11/one(x11) ≅ -x1
     x22 / x22 ≅ one(x22)
     x22/x1 ≅ x22
     x4/x1 ≅ x4
     float(x1/x4) == float(x1)/float(x4)
     convert(SurrealFinite, 6)/ convert(SurrealFinite, 3) ≅ convert(SurrealFinite, 2)

#### Simple functions

Many of the simple numerical functions will work, but not the more
advanced ones such as `log`. The two specific to surreals are
`generation` which returns the generation or birth day of a surreal,
and `canonicalise` which reduces it to the canonical form of the
surreal. 

     generation( zero(SurrealFinite) ) == 0
     canonicalise( convert(SurrealFinite,2)*convert(SurrealFinite,2) ) == convert(SurrealFinite,4)
     
#### Converting back to standard real

Examples of how to convert a number back to rationals or floats. 

     convert(Rational, convert(SurrealFinite, 1//8)) == 1//8
     float( x4 ) == -0.5
     convert( Rational, SurrealFinite( [ 7//16 ], [ 15//16 ] ) ) == 1//2
     convert(String, x4) == "-1/2"

And of how conversion is automatically applied by promotion rules

     x0 <= 1.0
     x1 == 1//1
     1.0 + x1 == 2.0
     canonicalise(1//2 + x1) == 1.5
     1//2 + x1 ≅ 1.5

## Final notes

The total implementation here is less than 1,000 lines of code. No doubt
an expert in Julia could make it a good deal tighter -- I have
concentrated on making the code easy (for me) to understand rather
than super concise. 

This little project would have been a lot easier if I knew more Julia,
or more about the surreals. Trying to build something to learn about
two moving parts at once wasn't a brilliant idea :)

But I was correct in thinking that (at least for me) this would have
been almost impossible to build in Matlab. And at the least, it would
have required a good deal more work without all the automagical pieces
of Julia helping.




