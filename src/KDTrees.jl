module KDTrees

    import Base.show

    using Compat
    using Devectorize

    export KDTree
    export knn, inball

    # Deprecated!
    export k_nearest_neighbour, query_ball_point

    include("hyper_rec.jl")
    include("kd_tree.jl")
    include("deprecations.jl")

end
