# auto-generated tests from julia-repl docstrings
using Test, Gapjm, AffineA
function mytest(file::String,cmd::String,man::String)
  println(file," ",cmd)
  exec=repr(MIME("text/plain"),eval(Meta.parse(cmd)),context=:limit=>true)
  if endswith(cmd,";") exec="nothing" 
  else exec=replace(exec,r"\s*$"m=>"")
       exec=replace(exec,r"\s*$"s=>"")
  end
  if exec!=man 
    i=1
    while i<=lastindex(exec) && i<=lastindex(man) && exec[i]==man[i]
      i=nextind(exec,i)
    end
    print("exec=$(repr(exec[i:end]))\nmanl=$(repr(man[i:end]))\n")
  end
  exec==man
end
@testset verbose = true "Gapjm" begin
@testset "AffineA.jl" begin
@test mytest("AffineA.jl","W=Atilde(3)","Atilde(3)")
@test mytest("AffineA.jl","l=elements(W,3)","9-element Vector{PPerm}:\n (2,3₋₁)\n (1,2)₁(3)₋₁\n (1)₋₁(2,3)₁\n (1,3)\n (1,3)₁(2)₋₁\n (1,2₋₁)₋₁(3)₁\n (1)₁(2,3₋₁)₋₁\n (1,3₋₁)₋₁(2)₁\n (1,2₋₁)")
@test mytest("AffineA.jl","mod1.([-1, 1, 6],3)","3-element Vector{Int64}:\n 2\n 1\n 3")
@test mytest("AffineA.jl","l[6]","(1,2₋₁)₋₁(3)₁")
@test mytest("AffineA.jl","B=DualBraidMonoid(W)","DualBraidMonoid(Atilde(3),c=[1, 2, 3])")
@test mytest("AffineA.jl","b=prod(B.(refls(W,1:7)))","c.4.5.46")
@test mytest("AffineA.jl","refls(W,2:5)","4-element Vector{PPerm}:\n (2,3)\n (1,3₋₁)\n (1,3)\n (1,2₋₁)")
@test mytest("AffineA.jl","PPerm([-1,1,6])","(1,2₋₁)₋₁(3)₁")
@test mytest("AffineA.jl","PPerm(3,(1,2=>-1)=>-1,3=>1)","(1,2₋₁)₋₁(3)₁")
@test mytest("AffineA.jl","cycles(PPerm([-1,1,6]))","2-element Vector{Pair{Vector{Int64}, Int64}}:\n [1, -1] => -1\n     [3] => 1")
end
end
