module TestKDTrees
    using FactCheck
    using Devectorize
    using Base.Test
    using Base.Collections

    using KDTrees

    FactCheck.onlystats(true)

    # Test helper to see that two arrays are equal not counting order.
    # This is probably done in a better way and having this function
    # here makes me look dumb.
    function allin(idxs, data)
        return length(symdiff(Set(idxs), Set(data))) == 0
    end


    include("test_knn.jl")
    include("test_build.jl")
    include("test_inball.jl")
    include("test_monkey.jl")

    FactCheck.exitstatus()

end #module
