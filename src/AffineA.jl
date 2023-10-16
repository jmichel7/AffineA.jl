"""
This package implements:

  - the type `PPerm`, periodic permutations of the integers.

  - the type `Atilde`, the Coxeter group `Ãₙ` as a group of `PPerm`.

  - the function `DualBraidMonoid` for such groups.

 It is based on the papers

  - [Digne, F.] "Presentations duales pour les groupes de tresses de type
    affine A", Comment. Math. Helv. 81 (2006) 23--47

  - [Shi] The Kazhdan-Lusztig cells in certain affine Weyl groups 
    Springer LNM 1179 (1986) 

©  2007 François Digne for the  mathematics, François Digne and Jean Michel
for the code.

An example:
```julia-repl
julia> W=Atilde(3) # The group Ã₂ as periodic permutations of period 3.
Atilde(3)

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

julia> l[6] # printed as cycles with shift -1 and 1 (sum 0). 2₋₁ is 2-3
(1,2₋₁)₋₁(3)₁

julia> B=DualBraidMonoid(W) # The Coxeter element is W(1,2,3)
DualBraidMonoid(Atilde(3),c=[1, 2, 3])

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
"""
module AffineA
using Gapjm
export PPerm, Atilde

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
  function PPerm(d::AbstractVector{<:Integer};check=false)
    if check validate(d) end
    new(convert(Vector{Int16},d))
  end
end

"check the validity of a `PPerm`"
function validate(d)
  n=length(d)
  if sort(mod.(d,n))!=0:n-1 error(d,": images must be distinct mod ",n) end
  if sum(d)!=sum(1:n) error(d,": sum of shifts is ",sum(d)-sum(1:n)," must be 0") end
end

Base.one(p::PPerm)=PPerm(1:length(p.d))
Base.isone(p::PPerm)=p.d==eachindex(p.d)
Base.copy(p::PPerm)=PPerm(copy(p.d))
Base.broadcastable(p::PPerm)=Ref(p)

perm(a::PPerm)=mod1.(a.d,length(a.d))

"""
`PPerm(n,c₁,…,cₗ;check=true)`  constructs a `PPerm` of period `n` by giving
its decomposition in cycles.

cycles `cᵢ` are given as pairs `(i₁,…,iₖ)=>d` representing
the  permutation `i₁↦ i₂↦ …↦ iₖ↦ i₁+d*n`.  `=>d` can be omitted when `d==0`
and `(i₁,)=>d` can be abbreviated to `i₁=>d`. 

An `iⱼ` itself may be given as a pair `v=>d` representing `v+n*d`.

The cycles must be disjoint `mod. n`. The argument is tested for validity if
`check=true`; in particular the cycles must be disjoint `mod. n`.

```julia-repl
julia> PPerm([-1,1,6])
(1,2₋₁)₋₁(3)₁

julia> PPerm(3,(1,2=>-1)=>-1,3=>1)
(1,2₋₁)₋₁(3)₁
```
"""
function PPerm(n::Integer,cc...;check=true)
  if isempty(cc) return PPerm(1:n) end
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
    PPerm(perm)
  end
  if check validate(p.d) end
  p
end

function Base.:*(x::PPerm,y::PPerm)
  n=length(x.d)
  res=similar(y.d)
  for i in 1:n 
    u=mod1(y.d[i],n)
    res[i]=x.d[u]+y.d[i]-u
  end
  PPerm(res)
end

function Base.inv(x::PPerm)
  n=length(x.d)
  l=perm(x)
  ll=invperm(l)
  PPerm(ll.-@view (x.d.-l)[ll])
end

function Base.:^(i::Int,p::PPerm)
  y=1+mod(i-1,length(p.d))
  i+p.d[y]-y
end

Base.:^(q::PPerm,p::PPerm)= inv(p)*q*p
Base.:^(a::PPerm, n::Integer)=n>=0 ? Base.power_by_squaring(a,n) :
                                     Base.power_by_squaring(inv(a),-n)
Base.:/(a::PPerm,b::PPerm)=a*inv(b)
Base.:\(a::PPerm,b::PPerm)=inv(a)*b

Base.:(==)(a::PPerm,b::PPerm)=a.d==b.d
function Base.hash(a::PPerm, h::UInt)
  for (i,v) in pairs(a.d)
    h=hash(v,h)
  end
  h
end

"""
`cycles(a::PPerm)`

The non-trivial cycles of a PPerm; each cycle is return as a pair
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
  n=length(a.d)
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
    pp=cyc=>div(x-cyc[1],n)
    if trivial || length(cyc)>1 || pp[2]!=0 push!(res,pp) end
  end
  res
end    

function Base.show(io::IO,a::PPerm)
  n=length(a.d)
  function stringdec(d)
    if get(io,:sgn,false)
      d>0 ? "₊"^d : "₋"^(-d)
    else
      d==0 ? "" : stringind(io,d)
    end
  end
  if !get(io,:limit,false) && !get(io,:TeX,false) 
    print(io,"PPerm(",a.d,")");return
  end
  c=cycles(a)
  if isempty(c) print(io,"()")
  else
    for cc in c
      cyc,d=cc
      print(io,"(",join(map(cyc)do y
        x=mod1(y,n)
        string(x,stringdec(div(y-x,n)))
      end,","),")")
      print(io,stringdec(d))
    end
  end
end

##------------------------AtildeGroup----------------------------
##The following formula is from [Shi] Lemma 4.2.2
Base.length(w::PPerm)=sum(j->sum(i->abs(fld(j^w-i^w,length(w.d))),1:j-1),eachindex(w.d))

isrightdescent(w::PPerm,i)= i==length(w.d) ? w.d[i]>w.d[1]+length(w.d) : w.d[i]>w.d[1+i]

CoxGroups.isleftdescent(w::PPerm,i)=isrightdescent(w^-1,i)

# for this function see [Digne],2.8
function Perms.reflength(w::PPerm)
  n=length(w.d)
  function v(pp,i)
    pp=map(x->sum.(getindex.(Ref(pp),x)),partitions(eachindex(pp),i))
    if pp[1][1]<0 pp=-pp end
    return tally.(pp)
  end
  d=last.(cycles(w))
  p0=count(iszero,d)
  res=n+length(d)-2*p0-count(i->i==i^w,1:n)
  if p0==length(d) return res end
  pos=filter(>(0),d)
  neg=filter(<(0),d)
  m=min(length(pos),length(neg))
  r=m:-1:1
  res-2*r[findfirst(i->length(intersect(v(pos,i),v(neg,i)))>0,r)]
end

function refword(w::PPerm)
  n=length(w.d)
  cnt=0
  function ff(w) local s
    l=reflength(w)
    d=1
    while true
      for i in 1:n
        s=PPerm(n,(i,i+d))
        if reflength(s*w)<l return s end
      end
      d+=1
      if mod(d,n)==0 d+=1 end
    end
  end
  res=PPerm[]
  while !isone(w)
    s=ff(w)
    push!(res,s)
    w=s*w
  end
  res
end

# this is unsufficient!!!
isdualsimple(y::PPerm)=count(c->c[2]!=0,cycles(y)) in [0,2]

# descent  sets are encoded as a pair: a list of atoms, and a list of atoms
# (1,u)  representing all atoms (1,uᵢ). This uses lemma 2.20 of [Digne] and
# is valid only if a is a dual simple.
function dualleftdescents(a::PPerm)
  n=length(a.d)
  if !isdualsimple(a) error(a," is not a dual simple") end
  res=[PPerm[],PPerm[]]
  for (x,d) in cycles(a)
    if x[1]!=1 || length(x)!=1
      for j in 1:length(x)
        for k in j+1:length(x)
          push!(res[1],PPerm(n,(x[j],x[k])))
          if d!=0 push!(res[1],PPerm(n,(x[j]+n,x[k]))) end
        end
        if d!=0 push!(res[2],PPerm(n,(1,mod1(x[j],n)))) end
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
`refls(W::Atilde,i::Integer)`

returns  the `i`-th reflection of  `W`. Reflections `(a,bⱼ)` are enumerated
by  lexicographical order  of `(j,a,b-a)`  with `j`  positive; however when
`a>b`  this reflection  is printed  `(b,a₋ⱼ)`. `i`  can also be a `Vector`;
then the corresponding list of reflections is returned.
"""
function PermRoot.refls(W::Atilde,i::Integer)
  n=ngens(W)
  p,r=divrem(i-1,n*(n-1))
  ecart,pos=divrem(r,n)
  res=collect(1:n)
  res[1+pos]=2+pos+ecart+p*n
  u=mod1(res[1+pos],n)
  res[u]=u-1-ecart-p*n
  PPerm(res)
end

# finds i such that pp=refls(W,i). It is assumed that pp is a reflection.
function whichatom(pp::PPerm)
  n=length(pp.d)
  pos=findfirst(i->pp.d[i]!=i,1:n)-1
  v=pp.d[pos+1]
  if v<pos+1
    pos=mod1(v,n)-1
    v=pp.d[pos+1]
  end
  p,ecart=divrem(v-pos-2,n)
  p*n*(n-1)+ecart*n+pos+1
end

Base.show(io::IO,G::Atilde)=print(io,"Atilde(",ngens(G),")")
  
"""
`Atilde(n::Integer)` returns `W(Ãₙ₋₁)` as a group of `PPerm` of period `n`.
"""
function Atilde(n::Integer)
  if n<2 error(n," should be >=2") end
  gens=map(i->PPerm(n,(i,i+1)),1:n-1)
  push!(gens,PPerm([0;2:n-1;n+1]))
  Atilde(gens,Dict{Symbol,Any}())
end

PermRoot.refltype(W::Atilde)=get!(W,:refltype)do
  n=ngens(W)
  [TypeIrred(Dict(:series=>Symbol("Ã"),:indices=>1:n,:rank=>n-1))]
end

Base.one(W::Atilde)=one(first(gens(W)))
Base.length(W::Atilde,w)=length(w)
CoxGroups.isleftdescent(W::Atilde,w,i)=isleftdescent(w,i)
Perms.reflength(W::Atilde,w)=reflength(w)
Base.isfinite(W::Atilde)=false

@GapObj struct AffaDualBraidMonoid{T,TW}<:Garside.GarsideMonoid{T}
  δ::T
  stringδ::String
  W::TW
end

Garside.IntervalStyle(M::AffaDualBraidMonoid)=Garside.Interval()

"""
`DualBraidMonoid(W::Atilde)`

If  `W=Atilde(n)`, constructs  the dual  braid monoid  for `Ãₙ₋₁`  and the
Coxeter element `PPerm([1-n;3:n;2+n])`
"""
function Garside.DualBraidMonoid(W::Atilde;c=(n=ngens(W);PPerm([1-n;3:n;2+n])),
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

Base.show(io::IO, M::AffaDualBraidMonoid)=print(io,"DualBraidMonoid(",M.W,",c=",
                                                word(M.W,M.δ),")")

function (M::AffaDualBraidMonoid)(r::PPerm)
  if !isdualsimple(r) error("r is not a dual simple") end
  if r==one(M) GarsideElt(M,PPerm[];check=false)
  elseif r==M.δ GarsideElt(M,PPerm[],1;check=false)
  else GarsideElt(M,[r];check=false)
  end
end

# a,b are results of dualleftdescents
function firstintersectionleftdescents(a,b)
  for t in a[1] if t in b[1] || perm(t) in map(x->x.d,b[2]) return t end end
  for t in b[1] if perm(t) in map(x->x.d,a[2]) return t end end
  for t in a[2] if t in b[2] return t end end
end

function Garside.leftgcdc(M::AffaDualBraidMonoid,a::PPerm,b::PPerm)
  x=one(M)
  while true
    t=firstintersectionleftdescents(dualleftdescents(a),dualleftdescents(b))
    if isnothing(t) return x,(a,b) end
    x*=t
    a=t\a
    b=t\b
  end
end

function CoxGroups.firstleftdescent(M::AffaDualBraidMonoid,w::PPerm)
  l1,l2=dualleftdescents(w)
  if !isempty(l1) return whichatom(first(l1)) end
  if !isempty(l2) return whichatom(first(l2)) end
end
  
Garside.δad(M::AffaDualBraidMonoid,x::PPerm,i::Integer)=iszero(i) ? x : x^(M.δ^i)

function Garside.atom(M::AffaDualBraidMonoid,i::Integer)
  s=refls(M.W,i)
  if !isdualatom(s) error("$s=atom($i) is not a dual atom") end
  s
end

function isdualatom(a::PPerm)
  n=length(a.d)
  if count(x->x[1]!=x[2],pairs(a.d))!=2 return false end
  i=findfirst(i->i!=a.d[i],1:n)
  i==1 || abs(i-a.d[i])<n
end

end
