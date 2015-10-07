@deprecate k_nearest_neighbour{T <: AbstractFloat}(tree::KDTree, point::Vector{T}, k::Int) knn(tree::KDTree, point::Vector{T}, k::Int)

@deprecate query_ball_point{T <: AbstractFloat}(tree::KDTree{T}, point::Vector{T}, radius::T) inball(tree::KDTree{T}, point::Vector{T}, radius::T)

@deprecate  KDTree{T <: AbstractFloat}(data::Matrix{T}, ls::Int, reord::Bool = true)  KDTree(data, leafsize = ls, reorder = reord)
