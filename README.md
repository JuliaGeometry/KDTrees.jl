# KDTrees

Kd trees for Julia.

[![Build Status](https://travis-ci.org/KristofferC/KDTrees.jl.svg?branch=master)](https://travis-ci.org/KristofferC/KDTrees.jl) [![Coverage Status](https://coveralls.io/repos/KristofferC/KDTrees.jl/badge.svg)](https://coveralls.io/r/KristofferC/KDTrees.jl)

Currently supports KNN-search and finding all points inside an hyper sphere centered at a given point

Care has been taken with regards to performance. The tree is for example not naively implemented as nodes pointing to other nodes but instead as a collection of densely packed arrays. This gives better cache locality. This
however means that the tree is immutable and new points can not be entered into the tree after it has been created.

There are some benchmarks for the creation of the tree and the different searches in the benchmark folder. 

## Author
Kristoffer Carlsson (@KristofferC)

## Examples

### Creating the tree

The tree is created with the command:
```julia
using KDTrees
data = rand(3,10^3)
tree = KDTree(data, leaf_size=5)
```
The `data` argument for the tree should be a matrix of floats of dimension `(n_dim, n_points)`. The `leaf_size` determines for what number of points the tree should stop splitting. 5 is a good number and is also the default value.

### Points inside hyper sphere

The exported `query_ball_point(tree, point, radius)` finds all points inside a hyper sphere centered at a given point with the given radius. The function
returns a sorted list of the indices of the points in the sphere.

```julia
using KDTrees
tree = KDTree(randn(3, 1000))
query_ball_point(tree, [0.0, 0.0, 0.0], 0.4)
```
gives the indices:
```
8-element Array{Int64,1}:
 184
 199
 307
 586
 646
 680
 849
 906
 926
```

### K-Nearest-Neighbours

The exported function `k_nearest_neighbour(tree, point, k)` finds the *k* nearest neighbours to a given point. The function returns a tuple of two lists with the indices and the distances from the given points respectively. These are sorted in the order of smallest to largest distance.

```julia
using KDTrees
tree = KDTree(randn(3, 1000))
k_nearest_neighbour(tree, [0.0, 0.0, 0.0], 5)
```
gives both the indices and distances:
```
([300,119,206,180,845],[0.052019,0.200885,0.220441,0.22447,0.235882])
```

## Benchmarks

The benchmarks have been made with computer with a 4 core Intel i5-2500K @ 3.3 GHz with Julia v0.4.0-dev+3034.

Clicking on a plot takes you to the Plotly site for the plot where the exact data can be seen.

### KNN benchmark

[![bench_knn](https://plot.ly/~kcarlsson89/346.png)](https://plot.ly/~kcarlsson89/346/)

### Build time benchmark

[![bench_knn](https://plot.ly/~kcarlsson89/168.png)](https://plot.ly/~kcarlsson89/168/)

## TODOs
* Add proper benchmarks, compare with others implementations. Update: Partly done
* Add other measures than Euclidean distance.
* Use a heap for storing the K best points in KNN instead of a linear array (should only matter for large K).

### Contribution

Contributions are more than welcome. If you have an idea that would make the tree have better 
performance or be more general please create a PR. Make sure you run the benchmarks before and
after your changes.
