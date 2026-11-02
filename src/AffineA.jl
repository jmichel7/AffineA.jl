"""
This package implements:
  - the type [`PPerm`](@ref), the periodic permutations of the integers.
  - the  function  [`coxeter_PPerm_group`](@ref)`(n)`  (or `coxPPerm(n)`), the
    Coxeter group of type `Ãₙ₋₁`, as a group of `PPerm`s of period `n`.
  - the function [`coxeter_sym_PPerm_group`](@ref)`(n)` (or `coxsymPPerm(n)`), the 
    Coxeter group of type `C̃_n` as a subgroup of `coxPPerm(2n)`.
  - the function [`DualBraidMonoid`](@ref) for such groups.
 It is based on the papers
  - [Digne] Presentations duales pour les groupes de tresses de type
    affine A, Comment. Math. Helv. 81 (2006) 23--47
  - [Digne2] A Garside presentation for Artin-Tits groups of type `C̃ₙ`,
    Ann. Inst. Fourier 62 (2012) 641--666 
  - [Shi] The Kazhdan-Lusztig cells in certain affine Weyl groups,
    Springer LNM 1179 (1986) 
©François Digne (2007) for the  mathematics, François Digne and Jean Michel
for the code.

### Installing
To install this package, start Julia, enter package mode with ], and then do
```
(@v1.12) pkg> add https://github.com/jmichel7/AffineA.jl
```
exit package mode with backspace and then do
```
julia> using Chevie, AffineA
```
and you are set up.

To update later to the latest version, do
```
(@v1.12) pkg> update AffineA
```
An example:
```julia-repl
julia> W=coxPPerm(3) # The group Ã₂ as periodic permutations of period 3.
coxeter_PPerm_group(3)

julia> l=elements(W,3) # elements of Coxeter length 3
9-element Vector{PPerm}:
 (1,3)
 (1)₁(2,3₋₁)₋₁
 (2,3₋₁)
 (1,3)₁(2)₋₁
 (1,2)₁(3)₋₁
 (1,3₋₁)₋₁(2)₁
 (1,2₋₁)
 (1,2₋₁)₋₁(3)₁
 (1)₋₁(2,3)₁
```
```julia-rep1
julia> print(l[8]) # the images of 1:3 by the 8th element
PPerm(Int16[-1, 1, 6])
```
```julia-repl
julia> mod1.([-1, 1, 6],3) # the image in 𝔖₃ of l[8]
3-element Vector{Int64}:
 2
 1
 3

julia> l[8] # printed as product of cycles with shifts -1 and 1 (sum 0). 2₋₁ is 2-3
PPerm(3): (1,2₋₁)₋₁(3)₁
```
other functions on PPerm are `*,/,\\,^,inv,==,copy,hash,isone,cycles,reflength`.
See the help for the last two functions and for `PPerm`.

The group `W` has the expected methods:
```julia-repl
julia> diagram(W)
       3
      / \
Ã₂   1———2


julia> cartan(W)
3×3 Matrix{Int64}:
  2  -1  -1
 -1   2  -1
 -1  -1   2

julia> reflrep(W)
3-element Vector{Matrix{Int64}}:
 [-1 0 0; 1 1 0; 1 0 1]
 [1 1 0; 0 -1 0; 0 1 1]
 [1 0 1; 0 1 1; 0 0 -1]

julia> refls(W,1:5) # the first 5 of the infinitely many reflections of W; see the help
5-element Vector{PPerm}:
 (1,2)
 (2,3)
 (1,3₋₁)
 (1,3)
 (1,2₋₁)

julia> B=BraidMonoid(W) # the Artin braid monid, a locally Garside monoid
BraidMonoid(coxeter_PPerm_group(3))

julia> B(1,2,3)^6
123123123123123123
```
The dual braid monoid is also implemented
```julia-repl
julia> B=DualBraidMonoid(W) # The Coxeter element is W(1,2,3)
DualBraidMonoid(coxeter_PPerm_group(3),c=[1, 2, 3])

julia> b=prod(B.(refls(W,1:7))) # the product in B of 7 dual atoms
c.4.5.46

julia> repr(b)
"B(2,4,3,4,5,4,6)"

julia> leftdescents(B,B.δ) # the infinitely many left descents are encoded (see help)
2-element Vector{Vector{PPerm}}:
 [(2,3), (2,3₋₁)]
 [(1,2), (1,3)]


julia> B(3,-1,2) # one can work in the group
(1)⁻¹6.2
```

Let `σ` be the automorphism of `W` which reverses the `n-1`first generators
```julia-repl
julia> σ.(gens(W))
3-element Vector{PPerm}:
 (2,3)
 (1,2)
 (1,3₋₁)
```
The group of `σ` fixed points in `coxPPerm(2n)` is isomorphic to `C̃ₙ`.
```julia-repl
julia> W1=coxsymPPerm(6)
coxeter_sym_PPerm_group(6)

julia> diagram(W1)
C̃₃ 1 > 2 - 3 < 4

julia> gens(W1)
4-element Vector{PPerm}:
 (1,6₋₁)
 (1,2)(5,6)
 (2,3)(4,5)
 (3,4)

julia> DualBraidMonoid(W1) # and it has also a dual braid monoid
DualBraidMonoid(coxeter_sym_PPerm_group(6),c=[1, 3, 2, 4])
```
"""
module AffineA
using Chevie
export PPerm, coxeter_PPerm_group, coxPPerm

"""
a `PPerm` represents a shiftless periodic permutation `f` of the integers
  - periodic of period `n` means `f(i+n)=f(i)+n`
  - then permutation means all `f(i)` are distinct mod `n`.
  - no shift means `sum(f.(1:n))==sum(1:n)`

it is represented in field `d` as the `Vector` `[f(1),…,f(n)]`. The default
constructor  takes a  vector of  integers, and  checks its  validity if the
keyword `check=true` is given.
"""
struct PPerm
  d::Vector{Int16}
  function PPerm(d::AbstractVector{<:Integer};check=true)
    if check validate(d) end
    new(convert(Vector{Int16},d))
  end
end

period(w::PPerm)=length(w.d)

"check the validity of a `PPerm`"
function validate(d)
  n=length(d)
  if sort(mod.(d,n))!=0:n-1 error(d,": images must be distinct mod ",n) end
  if sum(d)!=sum(1:n) error(d,": sum of shifts is ",sum(d)-sum(1:n)," must be 0") end
end

Base.one(p::PPerm)=PPerm(1:period(p);check=false)
Base.isone(p::PPerm)=p.d==eachindex(p.d)
Base.copy(p::PPerm)=PPerm(copy(p.d);check=false)
Base.broadcastable(p::PPerm)=Ref(p)

perm(a::PPerm)=mod1.(a.d,period(a))

"""
`PPerm(n,c₁,…,cₗ;check=true)`  constructs a `PPerm` of period `n` by giving
its decomposition into cycles.

cycles `cᵢ` are given as pairs `(i₁,…,iₖ)=>d` representing
the  permutation `i₁↦ i₂↦ …↦ iₖ↦ i₁+d*n`.  `=>d` can be omitted when `d==0`
and `(i₁,)=>d` can be abbreviated to `i₁=>d`. 

An `iⱼ` itself may be given as a pair `v=>d` representing `v+n*d`.

The  argument is  tested for  validity if  `check=true`; in  particular the
cycles must be disjoint `mod. n`.

When printed a cycle is normalised as `(i₁,i₂=>j₂,…,iₖ=>jₖ)=>j` with `1≤i₁<iₘ≤n` for all `m`.

```julia-repl
julia> PPerm([-1,1,6])
PPerm(3): (1,2₋₁)₋₁(3)₁

julia> PPerm(3,(1,2=>-1)=>-1,3=>1)
PPerm(3): (1,2₋₁)₋₁(3)₁

julia> PPerm(3,(1,-1)=>-1,3=>1)
PPerm(3): (1,2₋₁)₋₁(3)₁

julia> PPerm(3,(3,2=>1,4))
PPerm(3): (1,3₋₁,2)
```
"""
function PPerm(n::Integer,cc...;check=true)
  if isempty(cc) return PPerm(1:n;check=false) end
  cc=map(cc) do cyc
    c=if cyc isa Int 
      ((cyc,),0)
    elseif cyc isa Tuple 
      (cyc,0)
    elseif cyc isa Pair 
      if first(cyc) isa Int
        ((first(cyc),),last(cyc))
      else cyc
      end
    end
    map(c[1])do v
      if v isa Pair  
        v[1]+n*v[2]
      else v
      end
    end=>c[2]
  end
  if check
    u=collect(Iterators.flatten(map(x->mod.(x[1],n),cc)))
    if length(unique(u))!=length(u)
      error(cc," : the cycles must be disjoint mod ",n)
    end
  end
  p=prod(cc)do (cyc,d)
    perm=collect(1:n)
    for i in eachindex(cyc)
      k=mod1(cyc[i],n)
      perm[k]=cyc[mod1(i+1,length(cyc))]+k-cyc[i]
    end
    perm[mod1(cyc[end],n)]+=d*n
    PPerm(perm;check=false)
  end
  if check validate(p.d) end
  p
end

function Base.:*(x::PPerm,y::PPerm)
  n=period(x)
  res=similar(y.d)
  for i in 1:n 
    u=mod1(y.d[i],n)
    res[i]=x.d[u]+y.d[i]-u
  end
  PPerm(res;check=false)
end

Base.:isless(x::PPerm,y::PPerm)=x.d<y.d

function Base.inv(x::PPerm)
  l=perm(x)
  ll=invperm(l)
  PPerm(ll.-@view (x.d.-l)[ll];check=false)
end

function Base.:^(i::Int,p::PPerm)
  y=mod1(i,period(p))
  i+p.d[y]-y
end

Base.:^(q::PPerm,p::PPerm)= inv(p)*q*p
Base.:^(a::PPerm, n::Integer)=n>=0 ? Base.power_by_squaring(a,n) :
                                     Base.power_by_squaring(inv(a),-n)
Base.:/(a::PPerm,b::PPerm)=a*inv(b)
Base.:\(a::PPerm,b::PPerm)=inv(a)*b

Base.:(==)(a::PPerm,b::PPerm)=a.d==b.d
function Base.hash(a::PPerm, h::UInt)
  for (_,v) in pairs(a.d)
    h=hash(v,h)
  end
  h
end

"""
`cycles(a::PPerm)`

The non-trivial cycles of a PPerm; each cycle is returned as a pair
`(i₁,…,iₖ)=>d` and  is normalized such that `mod1(i₁,n)` is the smallest of the 
`mod1(iⱼ,n)` and `1≤i₁≤n`.

```julia-repl
julia> cycles(PPerm([-1,1,6]))
2-element Vector{Pair{Vector{Int64}, Int64}}:
 [1, -1] => -1
     [3] => 1
```
"""
function Perms.cycles(a::PPerm;trivial=false)
  n=period(a)
  res=Pair{Vector{Int},Int}[]
  l=trues(n)
  while true
    x=findfirst(l)
    if isnothing(x) break end
    cyc=Int[]
    while true
      push!(cyc,x)
      xn=mod1(x,n)
      l[xn]=false
      x+=a.d[xn]-xn
      if mod(x-cyc[1],n)==0 break end
    end
    shift=div(x-cyc[1],n)
    if trivial || length(cyc)>1 || shift!=0 push!(res,cyc=>shift) end
  end
  res
end    

# same as last.(cycles(w;trivial=true)) but faster
function shifts(a::PPerm)
  n=period(a)
  res=Int[]
  l=trues(n)
  while true
    x=findfirst(l)
    if isnothing(x) break end
    cyc=x
    while true
      xn=mod1(x,n)
      l[xn]=false
      x+=a.d[xn]-xn
      if mod(x-cyc,n)==0 break end
    end
    shift=div(x-cyc,n)
    push!(res,shift)
  end
  res
end    

function Base.show(io::IO, ::MIME"text/plain", a::PPerm)
  if !haskey(io,:typeinfo) print(io,"PPerm(",period(a),"): ") end
  show(io,a)
end

function Base.show(io::IO,a::PPerm)
  n=period(a)
  TeX=get(io,:TeX,false)
  function stringshift(d)
    if get(io,:sgn,false)
      if TeX 
        if d==0 return "" end
        res=d>0 ? "+"^d : "-"^(-d)
        if abs(d)>1 res="{"*res*"}" end
        "_"*res
      else d>0 ? "₊"^d : "₋"^(-d)
      end
    else
      d==0 ? "" : stringind(io,d)
    end
  end
  if !get(io,:limit,false) && !TeX
    print(io,"PPerm(",a.d,")");return
  end
  c=cycles(a)
  if isempty(c) print(io,"()")
  else
    for cc in c
      cyc,d=cc
      print(io,"(",join(map(cyc)do y
        x=mod1(y,n)
        string(x,stringshift(div(y-x,n)))
      end,","),")")
      print(io,stringshift(d))
    end
  end
end

##------------------------coxPPerm----------------------------
"""
`reflection_length(w::PPerm)` or `reflength`

returns the minimum number of transpositions of which `w` is a product.
see [2.8, Digne].
"""
function Perms.reflection_length(w::PPerm)
  n=period(w)
  d=shifts(w)
  p0=count(iszero,d)
  res=n+length(d)-2*p0
  if p0==length(d) return res end
  # we now compute κ(w)-p0
  pos=filter(>(0),d)
  neg=-filter(<(0),d)
  v(pp,k)=sort!.(map(x->sum.(x),partitions(pp,k)))
  for k in min(length(pos),length(neg)):-1:1
    if !isempty(intersect(v(pos,k),v(neg,k)))
      return res-2*k
    end
  end
end

function firstdescendingreflection(W,w)
  i=1
  l=reflength(w)
  while true
    s=refls(W,i)
    if reflength(s*w)<l return i end
    i+=1
  end
end

function refword(W,w::PPerm)
  res=Int[]
  while !isone(w)
    s=firstdescendingreflection(W,w)
    push!(res,s)
    w=refls(W,s)*w
  end
  res
end

#--------------- common code of AffaDualBraidMonoid and AffcDualBraidMonoid
abstract type AffDualBraidMonoid{T,TW}<:Garside.GarsideMonoid{T} end

Garside.IntervalStyle(::AffDualBraidMonoid)=Garside.Interval()

Base.show(io::IO, M::AffDualBraidMonoid)=print(io,"DualBraidMonoid(",M.W,
                                                ",c=",word(M.W,M.δ),")")

function CoxGroups.firstleftdescent(M::AffDualBraidMonoid,w::PPerm)
  l1,l2=leftdescents(M,w)
  if !isempty(l1) return whichrefl(M.W,first(l1)) end
  if !isempty(l2) return whichrefl(M.W,first(l2)) end
end
  
function Garside.leftgcdc(M::AffDualBraidMonoid,a::PPerm,b::PPerm)
  x=one(M)
  while true
    t=firstintersectionleftdescents(leftdescents(M,a),leftdescents(M,b))
    if isnothing(t) return x,(a,b) end
    x*=t
    a=t\a
    b=t\b
  end
end

function (M::AffDualBraidMonoid)(r::PPerm)
  if isdualsimple(M,r) GarsideElt(M,[r]) end
end

# a,b are results of leftdescents. This is valid in types Atilde  and Ctilde
function firstintersectionleftdescents(a,b)
  for t in a[1] if t in b[1] || perm(t) in map(x->x.d,b[2]) return t end end
  for t in b[1] if perm(t) in map(x->x.d,a[2]) return t end end
  for t in a[2] if t in b[2] return t end end
end

#--------------- AffaDualBraidMonoid ---------------------------------------
@GapObj struct AffaDualBraidMonoid{T,TW}<:AffDualBraidMonoid{T,TW}
  δ::T
  stringδ::String
  W::TW
end

#test if y divides δ 
function isdualsimple(M::AffaDualBraidMonoid,y::PPerm)
  period(y)==period(M.δ) && reflength(y)+reflength(M.δ/y)==period(y)
end

"""
`leftdescents(M::AffaDualBraidMonoid,a::PPerm)`

descent  sets are encoded as a pair: a list of atoms, and a list of atoms
`(1,u)`  representing all atoms `(1,uᵢ)`. This uses lemma 2.20 of [Digne] and
is valid only if a is a dual simple.
"""
function CoxGroups.leftdescents(M::AffaDualBraidMonoid,a::PPerm)
  n=period(a)
  if !isdualsimple(M,a) error(a," is not a dual simple") end
  res=[PPerm[],PPerm[]]
  for (x,d) in cycles(a)
    if x[1]!=1 || length(x)!=1
      for j in 1:length(x)
        for k in j+1:length(x)
          push!(res[1],PPerm(n,(x[j],x[k]);check=false))
          if d!=0 push!(res[1],PPerm(n,(x[j]+n,x[k]);check=false)) end
        end
        if d!=0 push!(res[2],PPerm(n,(1,mod1(x[j],n));check=false)) end
      end
    end
  end
  res
end

Perms.cycletype(a::PPerm)=sort(map(cyc->length(cyc[1])=>cyc[2],cycles(a)))

##--------------------------------------------------------------------------
@GapObj struct Atilde <: CoxeterGroup{PPerm}
  gens::Vector{PPerm}
end

"""
`refls(W::Atilde,i)`

returns  the  `i`-th  reflection  of  `W`.  Let  `n=ngens(W)`.  Reflections
`(a,bⱼ)`  are enumerated by  lexicographical order of  `(j,b-a,a)` with `j`
nonnegative; `b-a` varies from `1` to `n-1`, and `a` from `1` to `n`.

`i` can also be a vector of integers, the corresponding list of reflections
is then returned.
"""
function PermRoot.refls(W::Atilde,i::Integer)
  n=ngens(W)
  j,r=divrem(i-1,n*(n-1))
  ecart,pos=divrem(r,n) #pos=a-1, ecart=b-a-1
  res=collect(1:n)
  res[1+pos]=2+pos+ecart+j*n
  b=mod1(res[1+pos],n)
  res[b]=b-1-ecart-j*n
  PPerm(res;check=false)
end

"""
`whichrefl(W::Union{Atilde,Ctilde},pp::PPerm)`

returns  `i`  such  that  `pp==refls(W,i)`.  It  is  assumed that `pp` is a
reflection.
"""
function whichrefl(::Atilde,pp::PPerm)
  n=period(pp)
  pos=findfirst(i->pp.d[i]!=i,1:n)-1
  v=pp.d[pos+1]
  if v<pos+1
    pos=mod1(v,n)-1
    v=pp.d[pos+1]
  end
  p,ecart=divrem(v-pos-2,n)
  p*n*(n-1)+ecart*n+pos+1
end

Base.show(io::IO,G::Atilde)=print(io,"coxeter_PPerm_group(",ngens(G),")")
  
"""
`coxeter_PPerm_group(n::Integer)` returns `W(Ãₙ₋₁)` as a group of `PPerm` of period `n`.
"""
function coxeter_PPerm_group(n::Integer)
  if n<2 error(n," should be >=2") end
  gens=map(i->PPerm(n,(i,i+1);check=false),1:n-1)
  push!(gens,PPerm([0;2:n-1;n+1];check=false))
  Atilde(gens,Dict{Symbol,Any}())
end

const coxPPerm=coxeter_PPerm_group

PermRoot.refltype(W::Atilde)=get!(W,:refltype)do
  n=ngens(W)
  [TypeIrred(;series=Symbol("Ã"),indices=1:n,rank=n-1)]
end

Base.one(W::Atilde)=one(first(gens(W)))

"""
`length(::Atilde,w::PPerm)`

computes  the length of a periodic permutation using the formula from [Shi]
Lemma 4.2.2
"""
Base.length(::Atilde,w::PPerm)=
  sum(j->sum(i->abs(fld(j^w-i^w,period(w))),1:j-1;init=0),eachindex(w.d))

CoxGroups.isrightdescent(::Atilde,w::PPerm,i)=
  i==period(w) ? w.d[i]>w.d[1]+period(w) : w.d[i]>w.d[1+i]

CoxGroups.isleftdescent(W::Atilde,w::PPerm,i)=isrightdescent(W,w^-1,i)

Perms.reflength(::Atilde,w)=reflength(w)
Base.isfinite(::Atilde)=false

function Chevie.PermRoot.cartan(W::Atilde)
  r=ngens(W)
  if r==2 return [2 -2;-2 2] end
  res=cartan(:A,r)
  res[1,end]=-1
  res[end,1]=-1
  res
end

function PermRoot.reflrep(W::Atilde)
  get!(W,:matgens)do
    c=cartan(W)
    reflectionMatrix.(eachrow(one(c)),eachrow(c))
  end
end

PermRoot.reflrep(W::Atilde,w)=prod(reflrep(W)[word(W,w)];init=one(cartan(W)))
PermRoot.reflrep(W::Atilde,i::Integer)=reflrep(W)[i]

"""
`DualBraidMonoid(W::Atilde)`

If `W=coxeter_PPerm_group(n)`, constructs the dual braid monoid for `Ãₙ₋₁`
and the Coxeter element `c=PPerm([1-n;3:n;2+n])`. If
`M=DualBraidMonoid(W)`,  used as a  function, `M(w)` returns  an element of
the  dual  braid  monoid  if  `w∈W`  belongs  to  the interval `[1,c]`, and
`nothing` otherwise.
```julia-repl
julia> W=coxPPerm(3);l=elements(W,3)
9-element Vector{PPerm}:
 (1,3)
 (1)₁(2,3₋₁)₋₁
 (2,3₋₁)
 (1,3)₁(2)₋₁
 (1,2)₁(3)₋₁
 (1,3₋₁)₋₁(2)₁
 (1,2₋₁)
 (1,2₋₁)₋₁(3)₁
 (1)₋₁(2,3)₁

julia> B=DualBraidMonoid(W);B.(l)
9-element Vector{Union{Nothing, GarsideElt{PPerm, AffineA.AffaDualBraidMonoid{PPerm, AffineA.Atilde}}}}:
 4
 nothing
 6
 nothing
 nothing
 nothing
 5
 nothing
 c
```
"""
function Garside.DualBraidMonoid(W::Atilde;c=(n=ngens(W);PPerm([1-n;3:n;2+n];check=false)),
  revMonoid=nothing)
# If revMonoid is given, constructs the reversed monoid
# which allows to fill the field M.revMonoid after building the dual monoid
  M=AffaDualBraidMonoid(c,"c",W,Dict{Symbol,Any}())
  if revMonoid===nothing
    M.revMonoid=DualBraidMonoid(W;c=inv(c),revMonoid=M)
  else M.revMonoid=revMonoid
  end
  M
end

function Garside.atom(M::AffaDualBraidMonoid,i::Integer)
  s=refls(M.W,i)
  if !isdualatom(s) error("$s=atom($i) is not a dual atom") end
  s
end

function isdualatom(a::PPerm)
  if count(x->x[1]!=x[2],pairs(a.d))!=2 return false end
  n=period(a)
  i=findfirst(i->i!=a.d[i],1:n)
  i==1 || abs(i-a.d[i])<n
end
include("AffineC.jl")
end
