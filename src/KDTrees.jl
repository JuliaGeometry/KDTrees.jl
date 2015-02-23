module KDTrees

    import Base.show
    using Base.Collections

    using Compat
    using Devectorize

    export KDTree
    export k_nearest_neighbour,  query_ball_point

    include("kd_tree.jl")

end
