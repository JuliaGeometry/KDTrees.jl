module TestKDtree
    using FactCheck
    using Base.Test
    using Base.Collections

    using KDtree

    FactCheck.onlystats(true)

    include("test_knn.jl")
    include("test_query_ball.jl")

    FactCheck.exitstatus()

end #module
