# KDTrees

Kd trees for Julia.

[![Build Status](https://travis-ci.org/KristofferC/KDTrees.jl.svg?branch=master)](https://travis-ci.org/KristofferC/KDTrees.jl) [![Coverage Status](https://coveralls.io/repos/KristofferC/KDTrees.jl/badge.svg)](https://coveralls.io/r/KristofferC/KDTrees.jl)

This package contains a highly optimized kd tree to perform *k* nearest neighbour searches and range searches.

The readme contains some examples, different benchmarks and a comparison for kNN to scipy's cKDTree.

## Author
Kristoffer Carlsson (@KristofferC)

## Examples

### Creating the tree

The tree is created with the command:
```julia
using KDTrees
data = rand(3,10^3)
tree = KDTree(data, leafsize=10, reorder_data=true)
```
The `data` argument for the tree should be a matrix of floats of dimension `(n_dim, n_points)`. The `leafsize` determines for what number of points the tree should stop splitting. The default value is `leafsize = 10` which is a decent value. However, the optimal leafsize is dependent on the cost of the 
distance function used.

The `reorder_data` argument is a bool which determines if the input data should
be reordered to optimize for memory access. Points that are likely to be accessed close in time are also put close in memory. A copy is made of the data
so the original data given is untouched.

### Range searches

The exported `query_ball_point(tree, point, radius)` finds all points closer than the `radius` argument to the `point`. The function
returns a sorted list of the indices of the points in range.

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

The benchmarks have been made with computer with a 4 core Intel i5-2500K @ 4.2 GHz with Julia v0.4.0-dev+3034.

Clicking on a plot takes you to the Plotly site for the plot where the exact data can be seen.

### KNN benchmark

[![bench_knn](https://plot.ly/~kcarlsson89/397.png)](https://plot.ly/~kcarlsson89/397/)

### Build time benchmark

[![bench_knn](https://plot.ly/~kcarlsson89/413.png)](https://plot.ly/~kcarlsson89/413/)

### Short comparison vs scipy's cKDTree

One of the most popular packages for scientific computing in Python
is the scipy package. It can therefore be interesting to see how
KDTrees.jl compares against scipy's cKDTree.

A KNN search for a 100 000 point tree was performed for the five closest neighbours. The code and the resulting search speed are shown, first for
cKDTree and then for KDTrees.jl

**cKDTree:**

```python
>>> import numpy as np
>>> from scipy.spatial import cKDTree
>>> import timeit

>>> t = timeit.Timer("tree.query(queries, k=5)",
"""
import numpy as np
from scipy.spatial import cKDTree
data = np.random.rand(10**5, 3)
tree = cKDTree(data)
queries = np.random.rand(10**5, 3)
""")
>>> t = min(t.repeat(3, 10)) / 10

>>> print("knn / sec: ", 10**5 / t)
('knn / sec: ', 251394)
```

**KDTrees.jl:**
```julia
julia> tree = KDTree(rand(3,10^5));
julia> t = @elapsed for i = 1:10^5
       k_nearest_neighbour(tree, rand(3), 5)
       end;
julia> print("knn / sec: ", 10^5 / t)
knn / sec: 730922
```

### Contribution

Contributions are more than welcome. If you have an idea that would make the
tree have better performance or be more general please create a PR. Make 
sure you run the benchmarks before and after your changes.
