function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(KDTrees._in_ball, (KDTrees.KDTree{Float64}, Int64, Array{Float64, 1}, Float64, Array{Int64, 1}, Float64, Float64,))
    precompile(KDTrees._knn, (KDTrees.KDTree{Float64}, Array{Float64, 1}, Int64, Array{Int64, 1}, Array{Float64, 1}, Int64, Float64,))
    precompile(KDTrees.knn, (KDTrees.KDTree{Float64}, Array{Float64, 1}, Int64, Int64,))
    precompile(KDTrees.inball, (KDTrees.KDTree{Float64}, Array{Float64, 1}, Float64, Bool,))
    precompile(KDTrees.call, (Array{Any, 1}, Type{KDTrees.KDTree}, Array{Float64, 2},))
    precompile(KDTrees.percolate_down!, (Array{Float64, 1}, Array{Int64, 1}, Float64, Int64, Int64,))
    precompile(KDTrees.node_indices, (KDTrees.KDTree{Float64}, Int64,))
    precompile(KDTrees._knn_small, (KDTrees.KDTree{Float64}, Array{Float64, 1}, Int64, Array{Int64, 1}, Array{Float64, 1}, Int64,))
    precompile(KDTrees.addall, (KDTrees.KDTree{Float64}, Int64, Array{Int64, 1},))
    precompile(KDTrees.heap_sort_inplace!, (Array{Float64, 1}, Array{Int64, 1},))
    precompile(KDTrees.knn, (KDTrees.KDTree{Float64}, Array{Float64, 1}, Int64,))
end
