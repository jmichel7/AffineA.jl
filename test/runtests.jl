# auto-generated tests from julia-repl docstrings
using Test, Chevie, AffineA
function mytest(file::String,cmd::String,man::String)
  println(file," ",cmd)
  exec=repr(MIME("text/plain"),eval(Meta.parse(cmd)),context=:limit=>true)
  if endswith(cmd,";") return true end
  exec=replace(exec,r"\s*$"m=>""); exec=replace(exec,r"\s*$"s=>"")
  exec=replace(exec,r"^\s*"=>"")
  if exec==man return true end
  inds=collect(eachindex(exec))
  i=findfirst(i->i<=lastindex(man) && exec[i]!=man[i],inds)
  if i==nothing i=ncodeunits(man)+1 else i=inds[i] end
  print("exec=$(repr(exec[i:end]))\nmanl=$(repr(man[i:end]))\n")
  false
end
@testset "AffineA.jl" begin
@test mytest("AffineA.jl","W=coxPPerm(3)","coxeter_PPerm_group(3)")
@test mytest("AffineA.jl","l=elements(W,3)","9-element Vector{PPerm}:\n (1,3)\n (1)₁(2,3₋₁)₋₁\n (2,3₋₁)\n (1,3)₁(2)₋₁\n (1,2)₁(3)₋₁\n (1,3₋₁)₋₁(2)₁\n (1,2₋₁)\n (1,2₋₁)₋₁(3)₁\n (1)₋₁(2,3)₁")
@test mytest("AffineA.jl","mod1.([-1, 1, 6],3)","3-element Vector{Int64}:\n 2\n 1\n 3")
@test mytest("AffineA.jl","l[8]","PPerm(3): (1,2₋₁)₋₁(3)₁")
@test mytest("AffineA.jl","diagram(W)","3\n      / \\\nÃ₂   1———2")
@test mytest("AffineA.jl","cartan(W)","3×3 Matrix{Int64}:\n  2  -1  -1\n -1   2  -1\n -1  -1   2")
@test mytest("AffineA.jl","reflrep(W)","3-element Vector{Matrix{Int64}}:\n [-1 0 0; 1 1 0; 1 0 1]\n [1 1 0; 0 -1 0; 0 1 1]\n [1 0 1; 0 1 1; 0 0 -1]")
@test mytest("AffineA.jl","refls(W,1:5)","5-element Vector{PPerm}:\n (1,2)\n (2,3)\n (1,3₋₁)\n (1,3)\n (1,2₋₁)")
@test mytest("AffineA.jl","B=BraidMonoid(W)","BraidMonoid(coxeter_PPerm_group(3))")
@test mytest("AffineA.jl","B(1,2,3)^6","123123123123123123")
@test mytest("AffineA.jl","B=DualBraidMonoid(W)","DualBraidMonoid(coxeter_PPerm_group(3),c=[1, 2, 3])")
@test mytest("AffineA.jl","b=prod(B.(refls(W,1:7)))","c.4.5.46")
@test mytest("AffineA.jl","repr(b)","\"B(2,4,3,4,5,4,6)\"")
@test mytest("AffineA.jl","leftdescents(B,B.δ)","2-element Vector{Vector{PPerm}}:\n [(2,3), (2,3₋₁)]\n [(1,2), (1,3)]")
@test mytest("AffineA.jl","B(3,-1,2)","(1)⁻¹6.2")
@test mytest("AffineA.jl","σ.(gens(W))","3-element Vector{PPerm}:\n (2,3)\n (1,2)\n (1,3₋₁)")
@test mytest("AffineA.jl","W1=coxsymPPerm(6)","coxeter_sym_PPerm_group(6)")
@test mytest("AffineA.jl","diagram(W1)","C̃₃ 1 > 2 - 3 < 4")
@test mytest("AffineA.jl","gens(W1)","4-element Vector{PPerm}:\n (1,6₋₁)\n (1,2)(5,6)\n (2,3)(4,5)\n (3,4)")
@test mytest("AffineA.jl","DualBraidMonoid(W1)","DualBraidMonoid(coxeter_sym_PPerm_group(6),c=[1, 3, 2, 4])")
@test mytest("AffineA.jl","PPerm([-1,1,6])","PPerm(3): (1,2₋₁)₋₁(3)₁")
@test mytest("AffineA.jl","PPerm(3,(1,2=>-1)=>-1,3=>1)","PPerm(3): (1,2₋₁)₋₁(3)₁")
@test mytest("AffineA.jl","PPerm(3,(1,-1)=>-1,3=>1)","PPerm(3): (1,2₋₁)₋₁(3)₁")
@test mytest("AffineA.jl","PPerm(3,(3,2=>1,4))","PPerm(3): (1,3₋₁,2)")
@test mytest("AffineA.jl","cycles(PPerm([-1,1,6]))","2-element Vector{Pair{Vector{Int64}, Int64}}:\n [1, -1] => -1\n     [3] => 1")
@test mytest("AffineA.jl","W=coxPPerm(3);l=elements(W,3)","9-element Vector{PPerm}:\n (1,3)\n (1)₁(2,3₋₁)₋₁\n (2,3₋₁)\n (1,3)₁(2)₋₁\n (1,2)₁(3)₋₁\n (1,3₋₁)₋₁(2)₁\n (1,2₋₁)\n (1,2₋₁)₋₁(3)₁\n (1)₋₁(2,3)₁")
@test mytest("AffineA.jl","B=DualBraidMonoid(W);B.(l)","9-element Vector{Union{Nothing, GarsideElt{PPerm, AffineA.AffaDualBraidMonoid{PPerm, AffineA.Atilde}}}}:\n 4\n nothing\n 6\n nothing\n nothing\n nothing\n 5\n nothing\n c")
end
