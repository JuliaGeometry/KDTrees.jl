facts("KDTrees.build") do

# Issue # 15
context("KDTrees.build.15") do
    # Just check if it can be created
    @fact typeof(KDTree(ones(1, 11))) --> KDTree{Float64}
end

end
