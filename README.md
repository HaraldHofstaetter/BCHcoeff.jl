# BCH_series.jl

Efficient algorithms for computing coefficients of the Baker-Campbell-Hausdorff series.

This Julia package provides functions for computing coefficients
<img src="https://render.githubusercontent.com/render/math?math=h_w=\mathrm{coeff}(w,H)">
in the Baker-Campbell-Hausdorff series


<p align="center">
<img src="https://render.githubusercontent.com/render/math?math=H=\log(e^Ae^B)">
</p>

for words <img src="https://render.githubusercontent.com/render/math?math=w"> over the alphabet <img src="https://render.githubusercontent.com/render/math?math=\{A,B\}">,
or, more generally, in

<p align="center">
<img src="https://render.githubusercontent.com/render/math?math=H=\log(e^{A_1}{\cdots}e^{A_K})">
</p>

for words <img src="https://render.githubusercontent.com/render/math?math=w"> over the alphabet 
<img src="https://render.githubusercontent.com/render/math?math=\{A_1,\dots,A_K\}">.
 
A detailed description of the basic algorithm can be found in the appendix of

>[1] [H. Hofstätter](http://www.harald-hofstaetter.at), [Smallest common denominators for the homogeneous components of the Baker-Campbell-Hausdorff series](https://arxiv.org/pdf/2012.03818), submitted.  

## Example

```julia
julia> using BCH_series
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
 
julia> c = [BCH_coefficient(q) for q in Q]
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
 
 julia> BCH_coefficient([1,1,2,2,3,3,4,4,5,5,6,6,7,7], T=BigInt)
 1225927349791//277060715315634921630893140554547200000000
 
 julia> BCH_coefficient("ABAABBAAABBBAAAABBBBAAAAABBBBBAAAAAABBBBBBAAAAAAABBBBBBB")
 1225927349791//277060715315634921630893140554547200000000
 
 julia> BCH_coefficient([1,1,2,2,3,3,4,4,5,5,6,6,7,7], T=BigInt, Afirst=false)
 -1225927349791//277060715315634921630893140554547200000000
 
 julia> BCH_coefficient("BABBAABBBAAABBBBAAAABBBBBAAAAABBBBBBAAAAAABBBBBBBAAAAAAA")
 -1225927349791//277060715315634921630893140554547200000000
 
 julia> BCH_coefficient([1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3], [1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7], T=BigInt)
 -181365590859710868850573//3587775871359251968716509214115480195096319355948556615680000000000
 
 julia> BCH_coefficient("xyzxxyyzzxxxyyyzzzxxxxyyyyzzzzxxxxxyyyyyzzzzzxxxxxxyyyyyyzzzzzzxxxxxxxyyyyyyyzzzzzzz")
 -181365590859710868850573//3587775871359251968716509214115480195096319355948556615680000000000
```
