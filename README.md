# KDtree.jl

Kd tree for Julia.

[![Build Status](https://travis-ci.org/KristofferC/KDtree.jl.svg)](https://travis-ci.org/KristofferC/KDtree.jl) [![Coverage Status](https://coveralls.io/repos/KristofferC/KDtree.jl/badge.svg)](https://coveralls.io/r/KristofferC/KDtree.jl)

Currently supports KNN-search and finding all points inside an hyper sphere centered at a given point. Currently only
uses Euclidean distance.

Some care has been taken with regards to performance. For example the tree is not implemented as nodes pointing to other nodes but instead as an ensamble of densely packed arrays. This should give better cache locality. The negative aspect of this storage method is that the tree is immutable and new data can not be entered into the tree after it has been created.

There are some benchmarks for the creation of the tree and the different searches in the benchmark folder.

Since this is a new project there are still some obvious improvements which are listed in the TODO list.

## Author
Kristoffer Carlsson (@KristofferC)

## Requirements

Currently the "ArrayViews" package is required because `sub` in v0.3 of Julia
is quite bad.

## Examples

In the examples, notice that the module is called `KDtree` and the actual tree type is called `KDTree`. This is because modules and types can currently not have the same name in Julia.

The tree is created from a matrix of floats of dimension `(n_dim, n_points)`.

### Points inside hyper sphere

Finds all points inside an hyper sphere centered at a given point. Returns the indices of these points. 

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

### K-Nearest-Neighbours

Finds the *k* nearest neighbours to a given point. Returns a tuple of two lists with the indices and the distances
from the given points respectively. These are sorted in the order of smallest to largest distance.

The current implementation is a bit slower than it has to be for very large *k*.

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
* Add other measures than Euclidean distance.
* Use a bounded priority queue for storing the K best points in KNN instead of a linear array (should only matter for large K).
* Proper investigation of memory allocations and where time is spent.
* Throw errors at dimension mismatch in the functions etc.


