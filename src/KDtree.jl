module KDtree

    using ArrayViews

    export KDTree
    export k_nearest_neighbour,  query_ball_point

    include("kd_tree.jl")

end
