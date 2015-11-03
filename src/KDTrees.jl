VERSION >= v"0.4.0-dev+6521" && __precompile__()

module KDTrees

    import Base.show

    using Compat
    using Devectorize

    export KDTree
    export knn, inball

    # Deprecated!
    export k_nearest_neighbour, query_ball_point

    warn(string("KDTrees.jl is deprecated in favor of NearestNeighbors.jl.",
            " If you still want to use KDTrees.jl or rely on a package that",
            " uses it, please run Pkg.pin(\"KDTrees\", v\"0.5.13\") and this warning will go away"))

    include("hyper_rec.jl")
    include("kd_tree.jl")
    include("deprecations.jl")

end
