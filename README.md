# KDtree.jl

Kd tree for Julia.

[![Build Status](https://travis-ci.org/KristofferC/KDtree.jl.svg)](https://travis-ci.org/KristofferC/KDtree.jl) [![Coverage Status](https://coveralls.io/repos/KristofferC/KDtree.jl/badge.svg)](https://coveralls.io/r/KristofferC/KDtree.jl)

Currently supports KNN-search and finding all points inside an hyper sphere.

## Author
Kristoffer Carlsson (@KristofferC)

## Examples

### Points inside hyper sphere
```julia
using KDtree
tree = KDTree(randn(3, 1000))
query_ball_point(tree, [0.0, 0.0, 0.0], 0.4) # All nodes closer than 0.4 of (0.0, 0.0, 0.0)
```
gives the indices:
```
12-element Array{Int64,1}:
  66
 938
 227
 213
 542
 823
 586
  15
 178
 464
 887
 555
```

### KNN
```julia
using KDtree
tree = KDTree(randn(3, 1000))
k_nearest_neighbour(tree, [0.0, 0.0, 0.0], 5)
```
gives both the indices and distances:
```
([823,213,542,960,586],[0.207008,0.390674,0.372823,0.636221,0.224177])
```

### TODOs
* Implement a leaf size argument where the sub tree stop splitting after
   only a certain number of nodes are left in the sub tree.
* Add proper benchmarks, compare with others implementations.
* Priority Queue for storing the K best points in KNN instead of a linear array (should only matter for large K).
* Proper investigation of memory allocations and where time is spent.
* Throw errors at dimension mismatch in the functions etc.





