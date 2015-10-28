VERSION >= v"0.4.0-dev+6521" && __precompile__()

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
