####################################################################
# Hyper rectangles
####################################################################
# Hyper rectangles are used to bound points in space.
# For an inner node all it's children are bounded by the
# inner nodes hyper rectangle.
immutable HyperRectangle{T <: FloatingPoint}
    mins::Vector{T}
    maxes::Vector{T}
end


# When creating the tree we split the parents hyper rectangles
# so that each children gets its own.
function split_hyper_rec{T <: FloatingPoint}(hyper_rec::HyperRectangle{T},
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
# We do not return the sqrt here.
@inline function get_min_max_distance{T <: FloatingPoint}(rec::HyperRectangle{T}, point::Vector{T})
    min_d = zero(T)
    max_d = zero(T)
    @inbounds for dim in 1:size(point,1)
        d1 = (rec.maxes[dim] - point[dim]) * (rec.maxes[dim] - point[dim])
        d2 = (rec.mins[dim] - point[dim]) * (rec.mins[dim] - point[dim])
        if (rec.mins[dim] > point[dim]) || (rec.maxes[dim] < point[dim]) # Point is outside
            min_d += min(d1, d2)
        end
        max_d += max(d1,d2)
    end
    return min_d, max_d
end

####################################################################
# KD Tree
####################################################################
# The KDTree type
immutable KDTree{T <: FloatingPoint}
  data::Matrix{T} # dim x n_points array with floats
  hyper_recs::Vector{HyperRectangle{T}} # Each hyper rectangle bounds its children
  split_vals::Vector{T} # what values we split the tree at for a given internal node
  split_dims::Vector{Int} # what dimension we split the tree for a given internal node
  indices::Vector{Int} # Translates from a point index to the actual point in the data
  n_dim::Int
  last_node_size::Int
  leaf_size::Int
  n_leafs::Int
  n_internal_nodes::Int
  cross_node::Int # Index of first node in last row
  first_leaf_row::Int
  offset::Int
end


function show(io::IO, tree::KDTree)
    print(io, string("KDTree from ", size(tree.data, 2),
                 " point in ", size(tree.data, 1), " dimensions."))
end


# Helper functions to get node numbers and points
get_left_node(idx::Int) = idx * 2
get_right_node(idx::Int) = idx * 2 + 1
get_parent_node(idx::Int) = div(idx, 2)
is_leaf_node(tree::KDTree, idx::Int) = idx > tree.n_internal_nodes


# Gets number of points in a node, leaf_size for every node
# except the last one
function get_n_points(tree::KDTree, idx::Int)
    if idx != tree.n_leafs + tree.n_internal_nodes
        return tree.leaf_size
    else
        return tree.last_node_size
    end
end


# Gets the location of the node in the kd trees index vector
function get_point_index(tree::KDTree, idx::Int)
    if idx >= tree.cross_node
        return (idx - tree.cross_node) * tree.leaf_size + 1
    else
        return ((tree.offset - 1 + idx - 1 - tree.n_internal_nodes) * tree.leaf_size
                  + tree.last_node_size + 1)
    end
end


# Constructor for KDTree
function KDTree{T <: FloatingPoint}(data::Matrix{T},
                                    leaf_size::Int = 5)

    n_dim, n_points = size(data)

    n_leaf_nodes =  ceil(Integer, n_points / leaf_size)
    n_internal_nodes = n_leaf_nodes - 1

    l = floor(Integer, log2(n_leaf_nodes))

    offset = 2(n_leaf_nodes - 2^l)
    cross_node = int(2^(l+1))

    last_node_size = n_points % leaf_size
    if last_node_size == 0
        last_node_size = leaf_size
    end

    # This only happens when n_points / leaf_size is a power of 2?
    if cross_node >= n_internal_nodes + n_leaf_nodes
        cross_node = div(cross_node, 2)
    end

    indices = collect(1:n_points)
    split_vals = Array(T, n_internal_nodes)
    split_dims = Array(Int, n_internal_nodes) # 128 dimensions should be enough
    hyper_recs = Array(HyperRectangle{T}, n_internal_nodes + n_leaf_nodes)

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

    hyper_recs[1] = HyperRectangle{T}(mins, maxes)

    low = 1
    high = n_points

    # Call the recursive KDTree builder
    build_KDTree(1, data, split_vals, split_dims, hyper_recs,
                 indices, leaf_size, low, high)

    KDTree(data, hyper_recs, split_vals, split_dims, indices, n_dim, last_node_size,
          leaf_size, n_leaf_nodes, n_internal_nodes, cross_node, l ,offset)
end


# Recursive function to build the tree.
# Calculates what dimension has the maximum spread,
# and how many points to send to each side.
# Splits the hyper cubes and calls recursively
# with the new cubes and node indices.
function build_KDTree{T <: FloatingPoint}(index::Int,
                                          data::Matrix{T},
                                          split_vals::Vector{T},
                                          split_dims::Vector{Int},
                                          hyper_recs::Vector{HyperRectangle{T}},
                                          indices::Vector{Int},
                                          leaf_size::Int,
                                          low::Int,
                                          high::Int)

    n_points = high - low + 1 # Points left

    if n_points <= leaf_size
        return
    end

    # The number of leafs left, ceil to count a partially filled node as 1.
    n_leafs = ceil(Integer, n_points / leaf_size)

    # Rows left in the sub tree
    k = floor(Integer, log2(n_leafs))

    # Number of leftover nodes needed
    rest = n_leafs - 2^k

    # If the last leaf node will be on the right side of the tree we
    # send points so that left tree will be perfectly filled,
    # else we do the opposite.
    if k == 0
        mid_idx = low
    elseif n_points <= 2*leaf_size
        mid_idx = leaf_size + low
    elseif  rest > 2^(k-1)
        mid_idx = 2^k * leaf_size + low
    elseif rest == 0
        mid_idx = 2^(k-1)* leaf_size + low
    else
        mid_idx = n_points - 2^(k-1) * leaf_size + low
    end

    # Find the dimension where we have the largest spread.
    n_dim = size(data, 1)
    split_dim = -1
    max_spread = zero(T)
    for dim in 1:n_dim
        xmin = typemax(T)
        xmax = typemin(T)
        # Find max and min in this dim

        for coordinate in 1:n_points
            xmin = min(xmin, data[dim, indices[coordinate + low - 1]])
            xmax = max(xmax, data[dim, indices[coordinate + low - 1]])
        end

        if xmax - xmin > max_spread
            max_spread = xmax - xmin
            split_dim = dim
        end
    end

    split_dims[index] = split_dim

    # select! works like n_th element in c++
    # sorts the points in the maximum spread dimension s.t
    # data[split_dim, a]) < data[split_dim, b]) for all a > mid_idx, b > mid_idx
    select_spec!(indices, mid_idx, low, high, data, split_dim)

    split_val = data[split_dim, indices[mid_idx]]
    split_vals[index] = split_val

    # Create the hyper rectangles for the children
    hyper_rec_1, hyper_rec_2 = split_hyper_rec(hyper_recs[index], split_dim, split_val)
    hyper_recs[get_left_node(index)] = hyper_rec_1
    hyper_recs[get_right_node(index)] = hyper_rec_2

    build_KDTree(get_left_node(index), data,
                  split_vals, split_dims, hyper_recs,
                   indices, leaf_size, low, mid_idx - 1)

    build_KDTree(get_right_node(index), data,
                  split_vals, split_dims, hyper_recs,
                  indices, leaf_size, mid_idx, high)
end


####################################################################
# Distances
####################################################################
# Reduced euclidian distances
@inline function euclidean_distance_red{T <: FloatingPoint}(point_1::AbstractVector{T},
                                                        point_2::AbstractVector{T})
    dist = 0.0
    for i = 1:size(point_1, 1)
        @inbounds dist += abs2(point_1[i] - point_2[i])
    end
    return dist
end


@inline function euclidean_distance_red{T <: FloatingPoint}(tree::KDTree{T},
                                                        idx::Int,
                                                        point::AbstractVector{T})
    dist = 0.0
    for i = 1:tree.n_dim
        @inbounds dist += abs2(tree.data[i, tree.indices[idx]] - point[i])
    end
    return dist
end

####################################################################
# Query functions
####################################################################
# Finds the k nearest neighbour to a given point in space.
function k_nearest_neighbour{T <: FloatingPoint}(tree::KDTree, point::Vector{T}, k::Int)

    if k > size(tree.data, 2) || k <= 0
        error("k > number of points in tree or <= 0")
    end

    if size(point,1) != size(tree.data, 1)
        error(string("Wrong dimension of input point, points in the tree",
                     " have dimension ", size(tree.data, 1), " you",
                     " gave a point with dimension ", size(point,1), "."))
    end

    best_idxs = [-1 for i in 1:k+1]
    best_dists = [typemax(T) for i in 1:k+1]

    _k_nearest_neighbour(tree, point, k, best_idxs, best_dists, 1)

    # Sqrt here because storing in reduced distance format
    return best_idxs[1:end-1], sqrt(best_dists[1:end-1])
end

k_nearest_neighbour{P <: Real}(tree::KDTree, point::Vector{P}, k::Int) =
  k_nearest_neighbour(tree, float(point), k)


function _k_nearest_neighbour{T <: FloatingPoint}(tree::KDTree{T},
                                                  point::Vector{T},
                                                  k::Int,
                                                  best_idxs ::Vector{Int},
                                                  best_dists::Vector{T},
                                                  index::Int)

    if is_leaf_node(tree, index)
        point_index = get_point_index(tree, index)
        for z in point_index:point_index + get_n_points(tree, index) - 1
            dist_d = euclidean_distance_red(tree, z, point)
            if dist_d <= best_dists[1]
                best_dists[end] = dist_d
                idx = tree.indices[z]
                best_idxs[end] = idx
                percolate_down!(best_dists, best_idxs, dist_d, idx)
            end
        end
        return
    end

    if point[tree.split_dims[index]] < tree.split_vals[index]
        close_node = get_left_node(index)
        far_node = get_right_node(index)
    else
        far_node = get_left_node(index)
        close_node = get_right_node(index)
    end

    _k_nearest_neighbour(tree, point, k, best_idxs, best_dists, close_node)

    # Only go far node if the distance from the k-th best node crosses hyperplane
    if abs2(point[tree.split_dims[index]] - tree.split_vals[index]) < best_dists[1]
         _k_nearest_neighbour(tree, point, k, best_idxs, best_dists, far_node)
    end
    return
end


# Returns the sorted list of indices for all points in the tree inside a
# hypersphere of a given point with a given radius.
function query_ball_point{T <: FloatingPoint}(tree::KDTree{T},
                                              point::Vector{T},
                                              radius::T)

    if size(point,1) != size(tree.data, 1)
        error(string("Wrong dimension of input point, points in the tree",
                     " have dimension ", size(tree.data, 1), " you",
                     " gave a point with dimension ", size(point,1), "."))
    end

    index = 1
    idx_in_ball = Int[]
    traverse_check(tree, index, point, radius^2 , idx_in_ball)
    sort!(idx_in_ball)
    return idx_in_ball
end

query_ball_point{P <: Real, R <: Real}(tree::KDTree, point::Vector{P}, radius::R) =
  query_ball_point(tree, float(point), float(radius))


# Explicitly check the distance between leaf node and point while traversing
function traverse_check{T <: FloatingPoint}(tree::KDTree{T},
                                            index::Int,
                                            point::Vector{T},
                                            r::T,
                                            idx_in_ball::Vector{Int})

    min_d, max_d = get_min_max_distance(tree.hyper_recs[index], point)
    if min_d > r # Hyper sphere is outside hyper rectangle, skip the whole sub tree
        return
    end

    if is_leaf_node(tree, index)
        point_index = get_point_index(tree, index)
        for z in point_index:point_index + get_n_points(tree, index) - 1
            dist_d = euclidean_distance_red(tree, z, point)
            if dist_d < r
                push!(idx_in_ball, tree.indices[z])
            end
        end
        return
    end

    if max_d < r
        traverse_no_check(tree, index, idx_in_ball)
    else
        traverse_check(tree, get_left_node(index), point, r, idx_in_ball)
        traverse_check(tree, get_right_node(index), point, r, idx_in_ball)
    end
end


# Adds everything in this subtree since we have determined
# that the hyper rectangle completely encloses the hyper sphere
function traverse_no_check(tree::KDTree, index::Int, idx_in_ball::Vector{Int})
   if is_leaf_node(tree, index)
        point_index = get_point_index(tree, index)
        for z in point_index:point_index + get_n_points(tree, index) - 1
            push!(idx_in_ball, tree.indices[z])
        end
        return
    else
        traverse_no_check(tree, get_left_node(index), idx_in_ball)
        traverse_no_check(tree, get_right_node(index), idx_in_ball)
    end
end

####################################################################

# Taken from https://github.com/JuliaLang/julia/blob/v0.3.5/base/sort.jl
# and modified because I couldn't figure out how to get rid of
# the overhead when I passed in a new anonymous function
# to the "by" argument in each node. I also removed the return value.
function select_spec!{T <: FloatingPoint}(v::AbstractVector, k::Int, lo::Int,
                                          hi::Int, data::Matrix{T}, dim::Int)
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


# Binary min-heap percolate down.
function percolate_down!{T <: FloatingPoint}(xs::AbstractArray{T},
                                              xis::AbstractArray{Int},
                                              dist::T,
                                              index::Int)
    i = 1
    len = length(xs)
    @inbounds while (l = get_left_node(i)) <= len
        r = get_right_node(i)
        j = r >= len || (xs[l] > xs[r]) ? l : r
        if xs[j] > dist
            xs[i] = xs[j]
            xis[i] = xis[j]
            i = j
        else
            break
        end
    end
    xs[i] = dist
    xis[i] = index
end
