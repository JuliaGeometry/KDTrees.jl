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
8-element Array{Int64,1}:
 206
 180
 119
 199
 300
 186
 845
 161
```

### KNN
```julia
using KDtree
tree = KDTree(randn(3, 1000))
k_nearest_neighbour(tree, [0.0, 0.0, 0.0], 5)
```
gives both the indices and distances:
```
([300,119,206,180,845],[0.052019,0.200885,0.220441,0.22447,0.235882])
```

### TODOs
* Implement a leaf size argument where the sub tree stop splitting after
   only a certain number of nodes are left in the sub tree.
* Add proper benchmarks, compare with others implementations.
* Priority Queue for storing the K best points in KNN instead of a linear array (should only matter for large K).
* Proper investigation of memory allocations and where time is spent.
* Throw errors at dimension mismatch in the functions etc.





