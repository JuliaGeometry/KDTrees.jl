module KDTrees

    import Base.show

    using Compat
    using Devectorize

    export KDTree
    export k_nearest_neighbour,  query_ball_point

    include("hyper_rec.jl")
    include("kd_tree.jl")

end
