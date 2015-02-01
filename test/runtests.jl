module TestKDtree
    using FactCheck
    using Base.Test
    using KDtree

    FactCheck.onlystats(true)

    include("test_kd_tree.jl")

    FactCheck.exitstatus()

end #module
