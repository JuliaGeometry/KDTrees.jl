# Todo: update to mikowski, p = 1, inf
function euclidean_distance{T <: FloatingPoint}(point_1::Array{T, 1},
                                                point_2::Array{T, 1})
    dist = 0.0
    for i in 1:size(point_1, 1)
        dist += (point_1[i] - point_2[i]) * (point_1[i] - point_2[i])
    end
    return sqrt(dist)
end

# Hyper rectangles are used to bound points in space.
# For an inner node all it's children are bounded by the
# inner nodes hyper rectangle.
immutable HyperRectangle{T <: FloatingPoint}
    mins::Vector{T}
    maxes::Vector{T}
end


# When creating the tree we split the parents hyper rectangles
# so that each children gets its own.
function split_hyper_rec{T <: FloatingPoint}(hyper_rec::HyperRectangle,
                                             dim::Int,
                                             value::T)
    new_max = copy(hyper_rec.maxes)
    new_max[dim] = value

    new_min = copy(hyper_rec.mins)
    new_min[dim] = value

    return HyperRectangle(hyper_rec.mins, new_max),
               HyperRectangle(new_min, hyper_rec.maxes)
end


# From a hyper rectangle we can find the minimum and maximum distance to a point.
# If the point is inside the hyper cube the minimum dist is 0
function get_min_max_distance{T <: FloatingPoint}(rec::HyperRectangle, point::Vector{T})
    min_d = zero(T)
    max_d = zero(T)
    for dim in 1:size(point,1)
        d1 = (rec.maxes[dim] - point[dim]) * (rec.maxes[dim] - point[dim])
        d2 = (rec.mins[dim] - point[dim]) * (rec.mins[dim] - point[dim])
        if (rec.mins[dim] > point[dim]) || (rec.maxes[dim] < point[dim]) # Point is outside
            min_d += min(d1, d2)
        end
        max_d += max(d1,d2)
    end
    return sqrt(min_d), sqrt(max_d)
end


# The KDTree type
immutable KDTree{T <: FloatingPoint}
  data::Matrix{T} # dim x n_points array with floats
  split_vals::Vector{T} # what values we split the tree at for a given internal node
  split_dims::Vector{Int8} # what dimension we split the tree for a given internal node
  hyper_recs::Vector{HyperRectangle} # Each hyper rectangle bounds its children
  n_internal_nodes::Int
  indices::Vector{Int}
end


# Helper functions to get node numbers and points
get_left_node(idx::Int) = idx * 2
get_right_node(idx::Int) = idx * 2 + 1
get_parent_node(idx::Int) = div(idx, 2)
get_point_index(tree::KDTree, idx::Int) = idx - tree.n_internal_nodes

# From node index -> point in data
get_point(tree::KDTree, idx::Int) = tree.data[: , tree.indices[get_point_index(tree, idx)]]

is_leaf_node(tree::KDTree, idx::Int) = idx > tree.n_internal_nodes


# Constructor for KDTree
function KDTree{T <: FloatingPoint}(data::Matrix{T})

    n_dim, n_points = size(data)

    if n_dim > 20
        warn(string("You are sending in data with a large dimension, n_dim = ", n_dim, 
            ". K-d trees are not optimal for high dimensional data.", 
            " The data matrix should be given in dimensions (n_dim, n_points).",
            " Did you acidentally flip them?"))
    end

    n_internal_nodes = n_points - 1
    n_total_nodes = n_internal_nodes + n_points


    perm= collect(1:n_points)
    indices = Array(Int, n_points)
    split_vals = Array(T, n_internal_nodes)
    split_dims = Array(Int8, n_internal_nodes) # 128 dimensions should be enough
    hyper_recs = Array(HyperRectangle, n_total_nodes)

    # Create first bounding hyper rectangle
    maxes = Array(T, n_dim)
    mins = Array(T, n_dim)
    for j in 1:n_dim
        dim_max = typemin(T)
        dim_min = typemax(T)
        for k in 1:n_points
            dim_max = max(data[j, k], dim_max)
            dim_min = min(data[j, k], dim_min)
        end
        maxes[j] = dim_max
        mins[j] = dim_min
    end
    hyper_recs[1] = HyperRectangle(mins, maxes)

    # Call the recursive KDTree builder
    build_KDTree(1, data, sub(perm,1:length(perm)), split_vals,
                  split_dims, hyper_recs, n_internal_nodes, indices)

    KDTree(data, split_vals, split_dims,
           hyper_recs, n_internal_nodes,  indices)
  end


# Recursive function to build the tree.
# Calculates what dimension has the maximum spread,
# and how many points to send to each side.
# Splits the hyper cubes and calls recursively
# with the new cubes and node indices.
function build_KDTree{T <: FloatingPoint}(index::Int,
                                          data::Matrix{T},
                                          perm::SubArray{Int,1},
                                          split_vals::Vector{T},
                                          split_dims::Vector{Int8},
                                          hyper_recs::Vector{HyperRectangle},
                                          n_internal_nodes::Int,
                                          indices::Vector{Int})

    n_points = length(perm) # Points left
    n_dim = size(data, 1)

    if n_points <= 1
        indices[index - n_internal_nodes] = perm[1]
        return
    end

    # Find the dimension where we have the largest spread.
    split_dim = -1
    max_spread = zero(T)
    for dim in 1:n_dim
        xmin = typemax(T)
        xmax = typemin(T)
        # Find max and min in this dim

      for coordinate in 1:n_points
        xmin = min(xmin, data[dim, perm[coordinate]])
        xmax = max(xmax, data[dim, perm[coordinate]])
      end

      if xmax - xmin > max_spread # Found new max_spread, update split dimension
        max_spread = xmax - xmin
        split_dim = dim
      end
    end

    split_dims[index] = split_dim


    # Decide where to split
    # k = floor(Integer, log2(n_points)) for v 0.4
    k = ifloor(log2(n_points)) # <- deprecated in v 0.4
    rest = n_points - 2^k

    if rest > 2^(k-1)
        mid_idx = 2^k
    else
        mid_idx = 2^(k-1) + rest
    end

    # select! works like n_th element in c++
    # sorts the points in the maximum spread dimension s.t
    # data[split_dim, a]) < data[split_dim, b]) for all a > mid_idx, b > mid_idx
    select_spec!(perm, mid_idx, 1, length(perm), data, split_dim)

    split = data[split_dim, perm[mid_idx]]
    split_vals[index] = split

    # Create the hyper rectangles for the children
    hyper_rec_1, hyper_rec_2 = split_hyper_rec(hyper_recs[index], split_dim, split)
    hyper_recs[get_left_node(index)] = hyper_rec_1
    hyper_recs[get_right_node(index)] = hyper_rec_2

    build_KDTree(get_left_node(index), data,  sub(perm, 1:mid_idx),
                  split_vals, split_dims, hyper_recs, n_internal_nodes, indices )

    build_KDTree(get_right_node(index), data, sub(perm, mid_idx+1:length(perm)),
                  split_vals, split_dims, hyper_recs, n_internal_nodes, indices )
end


# Finds the k nearest neighbour to a given point in space.
function k_nearest_neighbour{T <: FloatingPoint}(tree::KDTree, point::Array{T, 1}, k::Int)

    if k > size(tree.data, 2) || k <= 0
        error("k > number of points in tree or <= 0")
    end

    if size(point,1) != size(tree.data, 1)
        error(string("Wrong dimension of input point, points in the tree",
                     " have dimension ", size(tree.data, 1), " you",
                     " gave a point with dimension ", size(point,1), "."))
    end

    best_idxs = [-1 for i in 1:k]
    best_dists = [typemax(T) for i in 1:k]

    _k_nearest_neighbour(tree, point, k, best_idxs, best_dists)

   # Convert from indices in tree to indices in data
    true_indices = [tree.indices[get_point_index(tree, x)] for x in best_idxs]

    return true_indices, best_dists
end


function _k_nearest_neighbour{T <: FloatingPoint}(tree::KDTree,
                                                  point::Array{T, 1},
                                                  k::Int,
                                                  best_idxs ::Array{Int, 1},
                                                  best_dists::Array{T, 1},
                                                  index::Int=1)

    min_d, max_d = get_min_max_distance(tree.hyper_recs[index], point)
    if min_d > best_dists[k]
       return
    end
    if is_leaf_node(tree, index)
        dist_d = euclidean_distance(get_point(tree, index), point)
        if dist_d <= best_dists[k] # Closer than the currently k closest.
            ins_point = 1
            while(ins_point < k && dist_d > best_dists[ins_point]) # Look through best_dists for insertion point
                ins_point +=1
            end
            for (i in k:-1:ins_point+1) # Shift down
                 best_idxs[i] = best_idxs[i-1]
                 best_dists[i]  = best_dists[i-1]
            end
            # Update with new values
            best_idxs[ins_point] = index
            best_dists[ins_point] = dist_d
        end
        return
    end

    if point[tree.split_dims[index]] < tree.split_vals[index]
        _k_nearest_neighbour(tree, point, k, best_idxs, best_dists, get_left_node(index))
        _k_nearest_neighbour(tree, point, k, best_idxs, best_dists, get_right_node(index))
    else
        _k_nearest_neighbour(tree, point, k, best_idxs, best_dists, get_right_node(index))
        _k_nearest_neighbour(tree, point, k,best_idxs, best_dists,  get_left_node(index))
    end
end
# Returns the indices for all points in the tree inside a
# hypersphere of a given point with a given radius
function query_ball_point{T <: FloatingPoint}(tree::KDTree,
                                              point::Vector{T},
                                              radius::T)

    if size(point,1) != size(tree.data, 1)
        error(string("Wrong dimension of input point, points in the tree",
                     " have dimension ", size(tree.data, 1), " you",
                     " gave a point with dimension ", size(point,1), "."))
    end

    index = 1
    idx_in_ball = Int[]
    traverse_check(tree, index, point, radius, idx_in_ball)
    return idx_in_ball
end


# Explicitly check the distance between leaf node and point
function traverse_check{T <: FloatingPoint}(tree::KDTree,
                                            index::Int,
                                            point::Array{T, 1},
                                            r::T,
                                            idx_in_ball::Vector{Int})

    min_d, max_d = get_min_max_distance(tree.hyper_recs[index], point)
    if min_d > r # Hyper shpere is outside hyper rectangle, skip the whole sub tree
        return
    elseif max_d < r
        traverse_no_check(tree, index, idx_in_ball)
    elseif is_leaf_node(tree, index)
        if euclidean_distance(get_point(tree, index), point) < r
            push!(idx_in_ball, tree.indices[get_point_index(tree, index)])
        end
        return
    else
        traverse_check(tree, get_left_node(index), point, r, idx_in_ball)
        traverse_check(tree, get_right_node(index), point, r, idx_in_ball)
    end
end


# Adds everything in this subtree since we have determined
# that the hyper rectangle completely encloses the hyper sphere
function traverse_no_check(tree::KDTree, index::Int, idx_in_ball::Vector{Int})
    if is_leaf_node(tree, index)
        push!(idx_in_ball, tree.indices[get_point_index(tree, index)])
        return
    else
        traverse_no_check(tree, get_left_node(index), idx_in_ball)
        traverse_no_check(tree, get_right_node(index), idx_in_ball)
      end
end


# Taken from https://github.com/JuliaLang/julia/blob/v0.3.5/base/sort.jl
# and modified because I couldn't figure out how to get rid of
# the memory consumption when I passed in a new anonymous function
# to the "by" argument in each node. I also removed the return value.
function select_spec!{T <: FloatingPoint}(v::AbstractVector, k::Int, lo::Int, hi::Int, data::Matrix{T}, dim::Int)
    lo <= k <= hi || error("select index $k is out of range $lo:$hi")
     @inbounds while lo < hi
        if hi-lo == 1
            if data[dim, v[hi]] < data[dim, v[lo]]
                v[lo], v[hi] = v[hi], v[lo]
            end
            return
        end
        pivot = v[(lo+hi)>>>1]
        i, j = lo, hi
        while true
            while data[dim, v[i]] < data[dim, pivot]; i += 1; end
            while  data[dim, pivot] <  data[dim, v[j]] ; j -= 1; end
            i <= j || break
            v[i], v[j] = v[j], v[i]
            i += 1; j -= 1
        end
        if k <= j
            hi = j
        elseif i <= k
            lo = i
        else
            return
        end
    end
    return
end


#=
# TODO: To avoid the sqrt maybe?
function euclidean_distance_red{T <: FloatingPoint}(point_1::Array{T, 1},
                                                    point_2::Array{T, 1})
  dist = 0.0
  for i in 1:size(point_1, 1)
    dist += (point_1[i] - point_2[i]) * (point_1[i] - point_2[i])
  end
  return dist
end
=#
