immutable HyperRectangle{T <: FloatingPoint}
  mins::Vector{T}
  maxes::Vector{T}
end


function split_hyper_rec{T <: FloatingPoint}(hyper_rec::HyperRectangle, dim::Int, value::T)
  new_max = copy(hyper_rec.maxes)
  new_max[dim] = value

  new_min = copy(hyper_rec.mins)
  new_min[dim] = value

  return HyperRectangle(hyper_rec.mins, new_max),
  HyperRectangle(new_min, hyper_rec.maxes)
end

function min_max_distance{T <: FloatingPoint}(rec::HyperRectangle, point::Vector{T})
  min_d = zero(T)
  max_d = zero(T)
  for dim in 1:size(point,1)
    d1 = (rec.maxes[dim] - point[dim]) * (rec.maxes[dim] - point[dim])
    d2 = (rec.mins[dim] - point[dim]) * (rec.mins[dim] - point[dim])
    if (rec.mins[dim] > point[dim]) || (rec.maxes[dim] < point[dim])
      min_d += min(d1, d2)
    end
    max_d += max(d1,d2)
  end
  return sqrt(min_d), sqrt(max_d)
end


immutable KDTree{T <: FloatingPoint}
  data::Matrix{T}
  split_vals::Vector{T}
  split_dims::Vector{Int8}
  hyper_recs::Vector{HyperRectangle} # Each hyper rectangle bounds a node
  n_internal_nodes::Int
  n_points::Int
  index2::Vector{Int}
end


function is_leaf_node(tree::KDTree, node_idx::Int)
  if node_idx > tree.n_internal_nodes
    return true
  end
  return false
end

function get_left_node(node_idx::Int)
  node_idx * 2
end

function get_right_node(node_idx::Int)
  node_idx * 2 + 1
end

function get_parent_node(node_idx::Int)
  div(node_idx, 2)
end

function get_point(tree::KDTree, node_idx::Int)
  #view(tree.data, : , tree.index2[get_point_index(tree, node_idx)])
  tree.data[: , tree.index2[get_point_index(tree, node_idx)]]
end

function get_point_index(tree::KDTree, node_idx::Int)
  node_idx - tree.n_internal_nodes
end

# Constructor
function KDTree{T <: FloatingPoint}(data::Matrix{T},
                                    leaf_size::Int=1)



  n_dim, n_points = size(data)
  n_internal_nodes = div(n_points, leaf_size) - 1
  if (n_internal_nodes % leaf_size != 0)
    n_internal_nodes += 1
  end
  n_total_nodes = n_internal_nodes + div(n_points,leaf_size) +
    ((n_internal_nodes % leaf_size != 0) ? 1 : 0)

  # 128 dimensions should be enough
  indices = collect(1:n_points)
  index2 = Array(Int, n_points)

  #println(T)
  split_vals = Array(T, n_internal_nodes)
  split_dims = Array(Int8, n_internal_nodes)

  hyper_recs = Array(HyperRectangle, n_total_nodes)
  #mins = Vector{T, n_dim}
  #maxes = Vector{T, n_dim}
  #nodes = Array(KDNode, n_total_nodes)

  # Create first bounding hyper rectangle
  maxes = Array(T, n_dim)
  mins = Array(T, n_dim)
  for j in 1:n_dim
    max_dim = typemin(T)
    min_dim = typemax(T)
    for k in 1:n_points
      max_dim = max(data[j, k], max_dim)
      min_dim = min(data[j, k], min_dim)
    end
    maxes[j] = max_dim
    mins[j] = min_dim
  end
  hyper_recs[1] = HyperRectangle(mins, maxes)


  build_KDTree!(1, data, sub(indices,1:length(indices)), split_vals,
                split_dims, hyper_recs, n_internal_nodes, index2)

  KDTree(data, split_vals, split_dims,
         hyper_recs, n_internal_nodes, n_points, index2)
end



function build_KDTree!{T <: FloatingPoint}(index::Int,
                                           data::Matrix{T},
                                           indices::SubArray{Int,1},
                                           split_vals::Vector{T},
                                           split_dims::Vector{Int8},
                                           hyper_recs::Vector{HyperRectangle},
                                           n_internal_nodes::Int,
                                           index2::Vector{Int})


  n_points = length(indices)
  n_dim = size(data, 1)

  if n_points <= 1
    index2[index - n_internal_nodes] = indices[1]
    return  # Only one index left
  end

  split_dim = 1
  max_spread = zero(T)
  for dim in 1:n_dim
    xmin = typemax(T)
    xmax = typemin(T)
    # Find max and min in this dim

    for coordinate in 1:n_points
      xmin = min(xmin, data[dim, indices[coordinate]])
      xmax = max(xmax, data[dim, indices[coordinate]])
    end

    if xmax - xmin > max_spread # Found new max_spread, update split dimension
      max_spread = xmax - xmin
      split_dim = dim
    end
  end

  split_dims[index] = split_dim

  # find the median and partition
  #mid_idx = div(length(indices), 2)
  k = ifloor(log2(length(indices))) # Yuck
  rest = length(indices) - 2^k

  mid_idx = 2^(k-1) +  ((rest>2^(k-1)) ?  2^(k-1) : rest)

  #mid_idx = max(mid_idx, length(indices) - mid_idx )
  select!(indices, mid_idx, by=i -> data[split_dim, i])
  split = data[split_dim, indices[mid_idx]]

  split_vals[index] = split

  hyper_rec_1, hyper_rec_2 = split_hyper_rec(hyper_recs[index], split_dim, split)

  hyper_recs[get_left_node(index)] = hyper_rec_1
  hyper_recs[get_right_node(index)] = hyper_rec_2

  build_KDTree!(get_left_node(index), data,  sub(indices, 1:mid_idx),
                split_vals, split_dims, hyper_recs, n_internal_nodes, index2 )

  build_KDTree!(get_right_node(index), data, sub(indices, mid_idx+1:length(indices)),
                split_vals, split_dims, hyper_recs, n_internal_nodes, index2 )
end






function nearest_neighbour{T <: FloatingPoint}(tree::KDTree, point::Array{T, 1})

  function _nearest_neighbour{T <: FloatingPoint}(tree::KDTree, point::Array{T, 1},
                                                  index::Int=1, best_dist::T=typemax(T),
                                                  best_point::Int=0)

    min_d, max_d = min_max_distance(tree.hyper_recs[index], point)
    if min_d > best_dist
      return best_point , best_dist
    end
    if is_leaf_node(tree, index)
      if best_point == 0
        best_point = index
        best_dist = euclidean_distance(get_point(tree, index), point)
      else
        dist_d = euclidean_distance(get_point(tree, index), point)
        if dist_d < best_dist
          best_point = index
          best_dist = dist_d
        end
      end
      return  best_point, best_dist
    end

    if point[tree.split_dims[index]] < tree.split_vals[index]
      best_point, best_dist = _nearest_neighbour(tree, point, get_left_node(index), best_dist, best_point)
      best_point, best_dist = _nearest_neighbour(tree, point, get_right_node(index), best_dist, best_point)
    else
      best_point, best_dist = _nearest_neighbour(tree, point, get_right_node(index), best_dist, best_point)
      best_point, best_dist = _nearest_neighbour(tree, point, get_left_node(index), best_dist, best_point)
    end
  end

  index, dist = _nearest_neighbour(tree, point)
  return tree.index2[get_point_index(tree, index)], dist
end




function euclidean_distance{T <: FloatingPoint}(point_1::Array{T, 1}, point_2::Array{T, 1})

  dist = 0.0
  for i in 1:size(point_1, 1)
    dist += (point_1[i] - point_2[i]) * (point_1[i] - point_2[i])
  end
  return sqrt(dist)
end

function euclidean_distance_red{T <: FloatingPoint}(point_1::Array{T, 1}, point_2::Array{T, 1})

  dist = 0.0
  for i in 1:size(point_1, 1)
    dist += (point_1[i] - point_2[i]) * (point_1[i] - point_2[i])
  end
  return dist
end

function  traverse_no_check(tree::KDTree, index::Int, indices::Vector{Int})
if is_leaf_node(tree, index)
  push!(indices, tree.index2[get_point_index(tree, index)])
  return
else
  traverse_no_check(tree, get_left_node(index), indices)
  traverse_no_check(tree, get_right_node(index), indices)
end
end

 function traverse_check{T <: FloatingPoint}(tree::KDTree, index::Int, point::Array{T, 1}, r::T, indices::Vector{Int})
    min_d, max_d = min_max_distance(tree.hyper_recs[index], point)
    if min_d > r
      return
    elseif (max_d < r) && (max_d > zero(T))
       traverse_no_check(tree, index, indices)
    elseif is_leaf_node(tree, index)
      if euclidean_distance(get_point(tree, index), point) < r
        push!(indices, tree.index2[get_point_index(tree, index)])
      end
      return
    else
      traverse_check(tree, get_left_node(index), point, r, indices)
      traverse_check(tree, get_right_node(index), point, r, indices)
    end
  end

function query_ball_point{T <: FloatingPoint}(tree::KDTree, point::Vector{T},
                                              radius::T)

  index = 1
  indices = Int[]
  traverse_check(tree, index, point, radius, indices)
  return indices
end








