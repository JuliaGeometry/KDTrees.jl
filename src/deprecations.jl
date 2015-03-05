@deprecate k_nearest_neighbour{T <: FloatingPoint}(tree::KDTree, point::Vector{T}, k::Int) knn(tree::KDTree, point::Vector{T}, k::Int)

@deprecate query_ball_point{T <: FloatingPoint}(tree::KDTree{T}, point::Vector{T}, radius::T) inball(tree::KDTree{T}, point::Vector{T}, radius::T)