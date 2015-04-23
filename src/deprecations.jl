@deprecate k_nearest_neighbour{T <: FloatingPoint}(tree::KDTree, point::Vector{T}, k::Int) knn(tree::KDTree, point::Vector{T}, k::Int)

@deprecate query_ball_point{T <: FloatingPoint}(tree::KDTree{T}, point::Vector{T}, radius::T) inball(tree::KDTree{T}, point::Vector{T}, radius::T)

@deprecate  KDTree{T <: FloatingPoint}(data::Matrix{T}, ls::Int, reord::Bool = true)  KDTree(data, leafsize = ls, reorder = reord)