# AffineA
Periodic permutations, Coxeter group `Ãₙ` and its dual braid monoid.

This  package  implements  the  periodic  permutations of the integers, the
Coxeter group of type `Ãₙ₋₁` as a group of periodic permutations of period
`n`, and the function dual braid monoid for such groups.

### Installing
To install this package, start Julia, enter package mode with ], and then do
```
(@v1.12) pkg> add https://github.com/jmichel7/AffineA.jl
```
exit package mode with backspace and then do
```
julia> using Chevie, AffineA
```
and you are set up. To update later to the latest version, do
```
(@v1.12) pkg> update AffineA
```
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jmichel7.github.io/AffineA.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jmichel7.github.io/AffineA.jl/dev/)
[![Build Status](https://github.com/jmichel7/AffineA.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jmichel7/AffineA.jl/actions/workflows/CI.yml?query=branch%3Amain)
