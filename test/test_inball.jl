facts("KDTrees.inball") do

context("KDTrees.inball.point") do

for rn in [true, false]
    data = [0.0 0.0 0.0 0.0 1.0 1.0 1.0 1.0;
            0.0 0.0 1.0 1.0 0.0 0.0 1.0 1.0;
            0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0] # 8 node cube

    tree = KDTree(data, leafsize = 2, reorder = rn)

    idxs = inball(tree, [1.1, 1.1, 1.1], 0.2)
    @fact idxs --> [8] # Only corner 8 at least 0.2 distance away from [1.1, 1.1, 1.1]

    idxs = inball(tree, [0.0, 0.0, 0.5], 0.6)
    @fact allin(idxs, [1, 2]) --> true# Corner 1 and 2 at least 0.6 distance away from [0.0, 0.0, 0.5]

    idxs = inball(tree, [0, 0, 0], 0.6)
    @fact idxs --> [1]

    idxs = inball(tree, [1//3, 1//3, 1//3], 1)
    @fact allin(idxs, [1, 2, 3, 5]) --> true

    idxs = inball(tree, [0.5, 0.5, 0.5], 0.2)
    @fact idxs --> [] #

    idxs = inball(tree, [0.5, 0.5, 0.5], 1.0)
    @fact allin(idxs, [1, 2, 3, 4, 5, 6, 7, 8]) --> true #

    @fact_throws inball(tree, [0.1], 1.0) # n_dim != trees dim
end
end # context

context("KDTrees.inball.tree") do

    data = [0.0 0.0 0.0 0.0 1.0 1.0 1.0 1.0;
            0.0 0.0 1.0 1.0 0.0 0.0 1.0 1.0;
            0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0] # 8 node cube

    tree1 = KDTree(data, leafsize = 2, reorder = false)
    tree2 = KDTree([[0.0, 0.0, 0.5] [0.0, 0.0, 0.0]], leafsize = 1, reorder = false)

    idxs = inball(tree1, tree2, 0.6)
    @fact allin(idxs[1], [1, 2]) --> true
    @fact allin(idxs[2], [1]) --> true
end # context

end # facts
