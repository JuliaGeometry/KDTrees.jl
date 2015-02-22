module TestKDTrees
    using FactCheck
    using Devectorize
    using Base.Test
    using Base.Collections

    using KDTrees

    FactCheck.onlystats(true)

    include("test_knn.jl")
    include("test_query_ball.jl")

    FactCheck.exitstatus()

end #module
