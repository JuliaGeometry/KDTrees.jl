# KDTrees

Kd trees for Julia.

[![Build Status](https://travis-ci.org/KristofferC/KDTrees.jl.svg?branch=master)](https://travis-ci.org/KristofferC/KDTrees.jl) [![Coverage Status](https://coveralls.io/repos/KristofferC/KDTrees.jl/badge.svg)](https://coveralls.io/r/KristofferC/KDTrees.jl)

Currently supports KNN-search and finding all points inside an hyper sphere centered at a given point. Currently only
uses Euclidean distance.

Some care has been taken with regards to performance. For example the tree is not implemented as nodes pointing to other nodes but instead as a collection of densely packed arrays. This should give better cache locality. The negative aspect of this storage method is that the tree is immutable and new data can not be entered into the tree after it has been created.

There are some benchmarks for the creation of the tree and the different searches in the benchmark folder. 

Since this is a new project there are still some obvious improvements which are listed in the TODO list.

## Author
Kristoffer Carlsson (@KristofferC)

## Examples

### Creating the tree

The tree is created with the command:
```julia
using KDTrees
data = rand(3,10^3)
tree = KDTree(data, leaf_size=15)
```
The `data` argument for the tree should be a matrix of floats of dimension `(n_dim, n_points)`. The `leaf_size` determines for what number of points the tree should stop splitting. 15 is a good number and is also the default value.

### Points inside hyper sphere

Finds all points inside an hyper sphere centered at a given point. This is done with the exported function `query_ball_point(tree, point, radius)`. Returns the sorted indices of these points. 

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

Finds the *k* nearest neighbours to a given point. his is done with the exported function `k_nearest_neighbour((tree, point, k)`. Returns a tuple of two lists with the indices and the distances
from the given points respectively. These are sorted in the order of smallest to largest distance.

The current implementation is a bit slower than it has to be for large *k*.

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

Clicking on a plot takes you to the Plotly site for the plot where the exact data can be seen.

### KNN benchmark

[![bench_knn](https://plot.ly/~kcarlsson89/168.png)](https://plot.ly/~kcarlsson89/168/)

### Build time benchmark

[![bench_knn](https://plot.ly/~kcarlsson89/164.png)](https://plot.ly/~kcarlsson89/164/)


## TODOs
* Add proper benchmarks, compare with others implementations. Update: Partly done
* Add other measures than Euclidean distance.
* Use a bounded priority queue for storing the K best points in KNN instead of a linear array (should only matter for large K). Julias built in PQ is slower than a normal array. 

### Contribution

Contributions are more than welcome. If you have an idea that would make the tree have better 
performance or be more general please create a PR. Make sure you run the benchmarks before and
after your changes.
