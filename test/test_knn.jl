facts("KDTrees.knn") do

context("KDTrees.knn") do

# 8 node rectangle
data = [0.0 0.0 0.0 0.5 0.5 1.0 1.0 1.0;
        0.0 0.5 1.0 0.0 1.0 0.0 0.5 1.0]
tree = KDTree(data)

idxs, dists = knn(tree, [0.8, 0.8], 1)
@fact idxs[1] --> 8 # Should be closest to top right corner
@fact sqrt(0.2^2 + 0.2^2) --> roughly(dists[1])

idxs, dists = knn(tree, [0.1, 0.8], 3)
@fact idxs --> [3, 2, 5]

idxs, dists = knn(tree, [1//10, 8//10], 3)
@fact idxs --> [3, 2, 5]

@fact_throws knn(tree, [0.1, 0.8], 10) # k > n_points

@fact_throws knn(tree, [0.1], 10) # n_dim != trees dim

end  #context

end # facts
