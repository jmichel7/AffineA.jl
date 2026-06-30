export coxsymPPerm, coxeter_sym_PPerm_group, σ
@GapObj struct Ctilde <: CoxeterGroup{PPerm}
  gens::Vector{PPerm}
end

"""
`coxeter_sym_PPerm_group(n::Integer)`

Returns the affine Coxeter group of type ``\\tilde C_{n/2}`` as the group of fixed
points under an automorphism `σ` in the affine Coxeter group of type `Ãₙ₋₁` itself
represented as a group of periodic permutations of period `n`.
The automorphism `σ` fixes the `n`-th generator and otherwise swaps the `i`-th and `n-i`-th generator.
"""
function coxeter_sym_PPerm_group(n::Integer)
  if n<4 || isodd(n) error(n," should be even >=4") end
  gens=[PPerm([0;2:n-1;n+1];check=false)]
  append!(gens,map(i->PPerm(n,(i,i+1),(n+1-i,n-i);check=false),1:div(n,2)-1))
  push!(gens,PPerm(n,(div(n,2),div(n,2)+1);check=false))
  Ctilde(gens,Dict{Symbol,Any}())
end

const coxsymPPerm=coxeter_sym_PPerm_group

PermRoot.refltype(W::Ctilde)=get!(W,:refltype)do
  n=ngens(W)
  [TypeIrred(Dict(:series=>Symbol("B̃"),:cartanType=>1,:indices=>1:n,:rank=>n-1))]
end

function Chevie.PermRoot.cartan(W::Ctilde)
  r=ngens(W)
  res=cartan(:A,r)*E(1)
  res[1,2]=res[2,1]=-root(2)
  res[end-1,end]=res[end,end-1]=-root(2)
  res
end

function PermRoot.reflrep(W::Ctilde)
  get!(W,:matgens)do
    c=cartan(W)
    reflectionMatrix.(eachrow(one(c)),eachrow(c))
  end
end

PermRoot.reflrep(W::Ctilde,w)=prod(reflrep(W)[word(W,w)];init=one(cartan(W)))
PermRoot.reflrep(W::Ctilde,i::Integer)=reflrep(W)[i]

Base.one(W::Ctilde)=one(first(gens(W)))
Base.isfinite(W::Ctilde)=false

Base.show(io::IO,G::Ctilde)=print(io,"coxeter_sym_PPerm_group(",2ngens(G)-2,")")

CoxGroups.isrightdescent(W::Ctilde,w::PPerm,i)= i==1 ? w.d[end]>w.d[1]+period(w) : w.d[i-1]>w.d[i]

CoxGroups.isleftdescent(W::Ctilde,w::PPerm,i)=isrightdescent(W,w^-1,i)

function σ(w::PPerm)
  m=period(w)
  res=similar(w.d)
  for i in 1:m
    res[i]=m+1-w.d[m+1-i]
  end
PPerm(res)
end

@GapObj struct AffcDualBraidMonoid{T,TW}<:Garside.GarsideMonoid{T}
  δ::T
  stringδ::String
  W::TW
end

Garside.IntervalStyle(M::AffcDualBraidMonoid)=Garside.Interval()

function
Garside.DualBraidMonoid(W::Ctilde;c=(n=ngens(W)-1;PPerm(2n,Tuple(1:2:2n-1)=>1,Tuple(2n:-2:2)=>-1;check=false)),
  revMonoid=nothing)
# If revMonoid is given, constructs the reversed monoid
# which allows to fill the field M.revMonoid after building the dual monoid
  M=AffcDualBraidMonoid(c,"c",W,Dict{Symbol,Any}())
  if revMonoid===nothing
    M.revMonoid=DualBraidMonoid(W;c=inv(c),revMonoid=M)
  else M.revMonoid=revMonoid
  end
  M
end

Base.show(io::IO, M::AffcDualBraidMonoid)=print(io,"DualBraidMonoid(",M.W,
                                                ",c=",word(M.W,M.δ),")")

function (M::AffcDualBraidMonoid)(r::PPerm)
  if reflength(r)+reflength(M.δ/r)==2ngens(M.W)-2 GarsideElt(M,[r]) end
end

function Garside.leftgcdc(M::AffcDualBraidMonoid,a::PPerm,b::PPerm)
  x=one(M)
  while true
    t=firstintersectionleftdescents(leftdescents(M,a),leftdescents(M,b))
    if isnothing(t) return x,(a,b) end
    x*=t
    a=t\a
    b=t\b
  end
end

"""
`refls(W::Ctilde,i::Integer)`

In `C̃ₙ` the reflections are `(a,b=>j)*σ((a,b=>j))`  if `a+b != 2n+1` and `(a,b=>j)` otherwise.
We obtain all reflections for `a` in `1:n` and `a< b≤ n+1-a` and `j>=0` if
for each `a,b` as above and `j>=0` we consider `(a,b=>j)` and `(a,b=>-j-1)` (there are
`2n²` reflections for each `j≥0`).
For each `j≥d` in 2 lexicographically ordered blocks (on `(a, b-a)`):
first the `(a,b=>j)`, then the `(a,b=>-j)`
for each `a`, `b-a` varies from `1` to `2(n-a)+1`: the number of ̀b̀ for ̀1≤ n-a ≤ x` is thus
`2(x+1)²`. 
"""
function PermRoot.refls(W::Ctilde,i::Integer)
  n=ngens(W)-1
  j,r=divrem(i-1,2n^2)
  r+=1
  if r <= n^2
    a=n-isqrt(r-1)
    b=a+r-isqrt(r-1)^2
    w=PPerm(2n,(a,b+j*2n))
  else
    r=r-n^2
    a=n-isqrt(r-1)
    b=a+r-isqrt(r-1)^2
    w=PPerm(2n,(a,b-(j+1)*2n))
  end  
  w==σ(w) ? w : w*σ(w)
end

# finds i such that pp=refls(W,i). It is assumed that pp=(a,b=>j)is a reflection.
function whichrefl(W::Ctilde,pp::PPerm)
  n=period(pp)
  a=findfirst(i->pp.d[i]!=i,1:n)
  m=div(n,2)
  j,b=divrem(pp.d[a],n)
   if b<=0 
     j-=1
     b+=n
   end
   if j>=0
     r=(m-a)^2+b-a
     return j*2m^2+r
   else
     r=(m-a)^2+b-a
     return -(2j+1)*m^2+r
   end  
end

function isdualsimple(M::AffcDualBraidMonoid,y::PPerm)
  if y!=σ(y) error(y, " is not in ", M.W) end
  reflength(y)+reflength(M.δ/y)==period(y)
end

"""
`leftdescents(M::AffcDualBraidMonoid,a::PPerm)`

Descent  sets are encoded as a pair: a list of atoms, and a list of atoms
`(a,b)`  representing all atoms `(a,bᵢ)`. This uses lemma 2.20 of [Digne] and
is valid only if a is a dual simple.
"""
function CoxGroups.leftdescents(M::AffcDualBraidMonoid,a::PPerm)
  n=period(a)
  if !isdualsimple(M,a) error(a," is not a dual simple") end
  res=[PPerm[],PPerm[]]
  for (x,d) in cycles(a)
    if d==0
      for j in 1:length(x)
        for k in j+1:length(x)
          push!(res[1],PPerm(n,(x[j],x[k]);check=false))
        end
      end
    end
    if d!=0
      for j in 1:length(x)
        for k in j+1:length(x)
           push!(res[1],PPerm(n,(x[j],x[k]);check=false)) 
           x[j]<x[k] ? push!(res[1],PPerm(n,(x[j]+n,x[k]);check=false)) : 
                       push!(res[1],PPerm(n,(x[k]+n,x[j]);check=false))
           push!(res[2],PPerm(n,(mod1(x[j],n),mod1(2n+1-x[k],n));check=false)) 
           push!(res[2],PPerm(n,(mod1(x[k],n),mod1(2n+1-x[j],n));check=false)) 
	end
        push!(res[2],PPerm(n,(mod1(x[j],n),mod1(2n+1-x[j],n));check=false)) 
      end
    end
  end
  res[1]=map(x->σ(x)==x ? x : x*σ(x),res[1])
  res[2]=map(x->σ(x)==x ? x : x*σ(x),res[2])
  [unique(res[1]),unique(res[2])]
end

function CoxGroups.firstleftdescent(M::AffcDualBraidMonoid,w::PPerm)
  l1,l2=leftdescents(M,w)
  if !isempty(l1) return whichrefl(M.W,first(l1)) end
  if !isempty(l2) return whichrefl(M.W,first(l2)) end
end

function Garside.atom(M::AffcDualBraidMonoid,i::Integer)
  s=refls(M.W,i)
  # if a reflection is simple it is an atom
  if !isdualsimple(M,s) error("$s=atom($i) is not a dual atom") end
  s
end

#"translation": on augmente de 1 les images paires des entiers impairs
# T^2 est la conjugaison par δ^n dans coxeter_sym_PPerm_group(n)
function T(w::PPerm)
  v=w.d
  x=copy(v)
  n=length(v)
  for i in 1:length(v)
    if isodd(i) && iseven(v[i])
      x[i]=w.d[i]-n
    elseif iseven(i) && isodd(v[i])
      x[i]=w.d[i]+n
    end  
  end
  return PPerm(x)
end


