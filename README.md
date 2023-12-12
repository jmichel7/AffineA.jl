
<a id='AffineA'></a>

<a id='AffineA-1'></a>

# AffineA

- [`AffineA`](index.md#AffineA)
- [`AffineA.PPerm`](index.md#AffineA.PPerm)
- [`AffineA.PPerm`](index.md#AffineA.PPerm-Tuple{Integer, Vararg{Any}})
- [`Chevie.Garside.DualBraidMonoid`](index.md#Chevie.Garside.DualBraidMonoid-Tuple{AffineA.Atilde})
- [`AffineA.coxeter_PPerm_group`](index.md#AffineA.coxeter_PPerm_group)
- [`Chevie.PermRoot.refls`](index.md#Chevie.PermRoot.refls-Tuple{AffineA.Atilde, Integer})
- [`PermGroups.Perms.cycles`](index.md#PermGroups.Perms.cycles-Tuple{PPerm})

<a id='AffineA' href='#AffineA'>#</a>
**`AffineA`** &mdash; *Module*.



This package implements:

  * the type `PPerm`, periodic permutations of the integers.
  * the  function  `coxeter_PPerm_group(n)`  (or `coxPPerm(n)`), the Coxeter group `Ãₙ₋₁` as a group of `PPerm` of period `n`.
  * the function `DualBraidMonoid` for such groups.

It is based on the papers

  * [Digne, F.] "Presentations duales pour les groupes de tresses de type affine A", Comment. Math. Helv. 81 (2006) 23–47
  * [Shi] The Kazhdan-Lusztig cells in certain affine Weyl groups  Springer LNM 1179 (1986)

©  2007 François Digne for the  mathematics, François Digne and Jean Michel for the code.

**Installing**

To install this package, at the Julia command line:

  * enter package mode with ]
  * do the command

```
(@v1.7) pkg> add "https://github.com/jmichel7/AffineA.jl"
```

  * exit package mode with backspace and then do

```
julia> using Chevie, AffineA
```

and you are set up.

To update later to the latest version, do

```
(@v1.7) pkg> update AffineA
```

An example:

```julia-repl
julia> W=coxPPerm(3) # The group Ã₂ as periodic permutations of period 3.
coxeter_PPerm_group(3)

julia> l=elements(W,3) # elements of Coxeter length 3
9-element Vector{PPerm}:
 (2,3₋₁)
 (1,2)₁(3)₋₁
 (1)₋₁(2,3)₁
 (1,3)
 (1,3)₁(2)₋₁
 (1,2₋₁)₋₁(3)₁
 (1)₁(2,3₋₁)₋₁
 (1,3₋₁)₋₁(2)₁
 (1,2₋₁)
```

```julia-rep1
julia> print(l[6]) # the images of 1:3 by the 6th element
PPerm(Int16[-1, 1, 6])
```

```julia-repl
julia> mod1.([-1, 1, 6],3) # the image in 𝔖₃
3-element Vector{Int64}:
 2
 1
 3

julia> l[6] # printed as product of cycles with shifts -1 and 1 (sum 0). 2₋₁ is 2-3
PPerm(3): (1,2₋₁)₋₁(3)₁

julia> B=DualBraidMonoid(W) # The Coxeter element is W(1,2,3)
DualBraidMonoid(coxeter_PPerm_group(3),c=[1, 2, 3])

julia> b=prod(B.(refls(W,1:7))) # the product in B of 7 dual atoms
c.4.5.46
```

```julia-rep1
julia> print(b)
B(2,4,3,4,5,4,6)
```

```julia-repl
julia> refls(W,2:5) # the corresponding reflections
4-element Vector{PPerm}:
 (2,3)
 (1,3₋₁)
 (1,3)
 (1,2₋₁)
```


<a target='_blank' href='https://github.com/jmichel7/AffineA.jl/blob/53aa6c6bfd271e9328d5b70440647c34b96ece33/src/AffineA.jl#L1-L91' class='documenter-source'>source</a><br>

<a id='AffineA.PPerm' href='#AffineA.PPerm'>#</a>
**`AffineA.PPerm`** &mdash; *Type*.



a `PPerm` represents a shiftless periodic permutation `f` of the integers

  * periodic of period `n` means `f(i+n)=f(i)+n`
  * then permutation means all `f(i)` are distinct mod `n`.
  * no shift means `sum(f.(1:n))==sum(1:n)`

it is represented in field `d` as the `Vector` `[f(1),…,f(n)]`. The default constructor  takes a  vector of  integers, and  checks its  validity if the keyword `check=true` is given.


<a target='_blank' href='https://github.com/jmichel7/AffineA.jl/blob/53aa6c6bfd271e9328d5b70440647c34b96ece33/src/AffineA.jl#L96-L105' class='documenter-source'>source</a><br>

<a id='AffineA.PPerm-Tuple{Integer, Vararg{Any}}' href='#AffineA.PPerm-Tuple{Integer, Vararg{Any}}'>#</a>
**`AffineA.PPerm`** &mdash; *Method*.



`PPerm(n,c₁,…,cₗ;check=true)`  constructs a `PPerm` of period `n` by giving its decomposition into cycles.

cycles `cᵢ` are given as pairs `(i₁,…,iₖ)=>d` representing the  permutation `i₁↦ i₂↦ …↦ iₖ↦ i₁+d*n`.  `=>d` can be omitted when `d==0` and `(i₁,)=>d` can be abbreviated to `i₁=>d`. 

An `iⱼ` itself may be given as a pair `v=>d` representing `v+n*d`.

The  argument is  tested for  validity if  `check=true`; in  particular the cycles must be disjoint `mod. n`.

```julia-repl
julia> PPerm([-1,1,6])
PPerm(3): (1,2₋₁)₋₁(3)₁

julia> PPerm(3,(1,2=>-1)=>-1,3=>1)
PPerm(3): (1,2₋₁)₋₁(3)₁

julia> PPerm(3,(1,-1)=>-1,3=>1)
PPerm(3): (1,2₋₁)₋₁(3)₁
```


<a target='_blank' href='https://github.com/jmichel7/AffineA.jl/blob/53aa6c6bfd271e9328d5b70440647c34b96ece33/src/AffineA.jl#L128-L151' class='documenter-source'>source</a><br>

<a id='PermGroups.Perms.cycles-Tuple{PPerm}' href='#PermGroups.Perms.cycles-Tuple{PPerm}'>#</a>
**`PermGroups.Perms.cycles`** &mdash; *Method*.



`cycles(a::PPerm)`

The non-trivial cycles of a PPerm; each cycle is returned as a pair `(i₁,…,iₖ)=>d` and  is normalized such that `mod1(i₁,n)` is the smallest of the  `mod1(iⱼ,n)` and `1≤i₁≤n`.

```julia-repl
julia> cycles(PPerm([-1,1,6]))
2-element Vector{Pair{Vector{Int64}, Int64}}:
 [1, -1] => -1
     [3] => 1
```


<a target='_blank' href='https://github.com/jmichel7/AffineA.jl/blob/53aa6c6bfd271e9328d5b70440647c34b96ece33/src/AffineA.jl#L227-L240' class='documenter-source'>source</a><br>

<a id='AffineA.coxeter_PPerm_group' href='#AffineA.coxeter_PPerm_group'>#</a>
**`AffineA.coxeter_PPerm_group`** &mdash; *Function*.



`coxeter_PPerm_group(n::Integer)` returns `W(Ãₙ₋₁)` as a group of `PPerm` of period `n`.


<a target='_blank' href='https://github.com/jmichel7/AffineA.jl/blob/53aa6c6bfd271e9328d5b70440647c34b96ece33/src/AffineA.jl#L412-L414' class='documenter-source'>source</a><br>

<a id='Chevie.PermRoot.refls-Tuple{AffineA.Atilde, Integer}' href='#Chevie.PermRoot.refls-Tuple{AffineA.Atilde, Integer}'>#</a>
**`Chevie.PermRoot.refls`** &mdash; *Method*.



`refls(W::Atilde,i::Integer)`

returns  the `i`-th reflection of  `W`. Reflections `(a,bⱼ)` are enumerated by  lexicographical order of `(j,a,b-a)` with `j` nonnegative; however when `a>b`  this reflection  is printed  `(b,a₋ⱼ)`. `i`  can also be a `Vector`; then the corresponding list of reflections is returned.


<a target='_blank' href='https://github.com/jmichel7/AffineA.jl/blob/53aa6c6bfd271e9328d5b70440647c34b96ece33/src/AffineA.jl#L378-L385' class='documenter-source'>source</a><br>

<a id='Chevie.Garside.DualBraidMonoid-Tuple{AffineA.Atilde}' href='#Chevie.Garside.DualBraidMonoid-Tuple{AffineA.Atilde}'>#</a>
**`Chevie.Garside.DualBraidMonoid`** &mdash; *Method*.



`DualBraidMonoid(W)`

If `W=coxeter_PPerm_group(n)`, constructs the dual braid monoid for `Ãₙ₋₁` and the Coxeter element `c=PPerm([1-n;3:n;2+n])`. If `M=DualBraidMonoid(W)`,  used as a  function, `M(w)` returns  an element of the dual braid monoid if `w` belongs to the interval `[1,c]`, and `nothing` otherwise.

```julia-repl
julia> W=coxPPerm(3);l=elements(W,3)
9-element Vector{PPerm}:
 (2,3₋₁)
 (1,2)₁(3)₋₁
 (1)₋₁(2,3)₁
 (1,3)
 (1,3)₁(2)₋₁
 (1,2₋₁)₋₁(3)₁
 (1)₁(2,3₋₁)₋₁
 (1,3₋₁)₋₁(2)₁
 (1,2₋₁)

julia> B=DualBraidMonoid(W);B.(l)
9-element Vector{Union{Nothing, GarsideElt{PPerm, AffineA.AffaDualBraidMonoid{PPerm, AffineA.Atilde}}}}:
 6
 nothing
 c
 4
 nothing
 nothing
 nothing
 nothing
 5
```


<a target='_blank' href='https://github.com/jmichel7/AffineA.jl/blob/53aa6c6bfd271e9328d5b70440647c34b96ece33/src/AffineA.jl#L443-L476' class='documenter-source'>source</a><br>

