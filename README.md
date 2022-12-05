# BCHcoeff.jl

This Julia package provides functions for computing coefficients
<img src="https://render.githubusercontent.com/render/math?math=h_w=\mathrm{coeff}(w,H)">
in the Baker-Campbell-Hausdorff series


<p align="center">
<img src="https://render.githubusercontent.com/render/math?math=H=\sum_wh_ww=\log(e^Ae^B)">
</p>

for words <img src="https://render.githubusercontent.com/render/math?math=w"> over the alphabet <img src="https://render.githubusercontent.com/render/math?math=\{A,B\}">,
or, more generally, in

<p align="center">
<img src="https://render.githubusercontent.com/render/math?math=H=\sum_wh_ww=\log(e^{A_1}{\cdots}e^{A_K})">
</p>

for words <img src="https://render.githubusercontent.com/render/math?math=w"> over the alphabet 
<img src="https://render.githubusercontent.com/render/math?math=\{A_1,\dots,A_K\}">.
 
A detailed description of the basic algorithm can be found in

>[1] [H. Hofstätter](http://www.harald-hofstaetter.at), [A simple and efficient algorithm for computing the
Baker-Campbell-Hausdorff series](https://arxiv.org/pdf/2212.01290).  

## Installation
```julia
julia> using Pkg; Pkg.add(PackageSpec(url="https://github.com/HaraldHofstaetter/BCHcoeff.jl"))
```

## Example

```julia
julia> using BCHcoeff
julia> using Combinatorics

julia> Q = collect(partitions(7))
15-element Array{Array{Int64,1},1}:
 [7]                  
 [6, 1]               
 [5, 2]               
 [5, 1, 1]            
 [4, 3]               
 [4, 2, 1]            
 [4, 1, 1, 1]         
 [3, 3, 1]            
 [3, 2, 2]            
 [3, 2, 1, 1]         
 [3, 1, 1, 1, 1]      
 [2, 2, 2, 1]         
 [2, 2, 1, 1, 1]      
 [2, 1, 1, 1, 1, 1]   
 [1, 1, 1, 1, 1, 1, 1]
 
julia> c = [bchcoeff(q) for q in Q]
15-element Array{Rational{Int64},1}:
  0//1    
  1//30240
 -1//5040 
 -1//5040 
  1//3780 
  1//2016 
  1//2016 
 -1//1512 
 -1//5040 
 -1//5040 
 -1//630  
 -1//1120 
  1//840  
  1//840  
 -1//140  
 
 julia> bchcoeff([1,1,2,2,3,3,4,4,5,5,6,6,7,7], T=BigInt)
 1225927349791//277060715315634921630893140554547200000000
 
 julia> bchcoeff("ABAABBAAABBBAAAABBBBAAAAABBBBBAAAAAABBBBBBAAAAAAABBBBBBB")
 1225927349791//277060715315634921630893140554547200000000
 
 julia> bchcoeff([1,1,2,2,3,3,4,4,5,5,6,6,7,7], T=BigInt, Afirst=false)
 -1225927349791//277060715315634921630893140554547200000000
 
 julia> bchcoeff("BABBAABBBAAABBBBAAAABBBBBAAAAABBBBBBAAAAAABBBBBBBAAAAAAA")
 -1225927349791//277060715315634921630893140554547200000000
 
 julia> bchcoeff([1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3], [1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7], T=BigInt)
 -181365590859710868850573//3587775871359251968716509214115480195096319355948556615680000000000
 
 julia> bchcoeff("xyzxxyyzzxxxyyyzzzxxxxyyyyzzzzxxxxxyyyyyzzzzzxxxxxxyyyyyyzzzzzzxxxxxxxyyyyyyyzzzzzzz")
 -181365590859710868850573//3587775871359251968716509214115480195096319355948556615680000000000
```
