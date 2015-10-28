####################################################################
# KD Tree
####################################################################
# The KDTree type
immutable KDTree{T <: AbstractFloat}
    data::Matrix{T} # dim x n_p array with floats
    hyper_recs::Vector{HyperRectangle{T}} # Each hyper rectangle bounds its children
    indices::Vector{Int} # Translation between tree index to point indes or adsa
    split_vals::Vector{T} # what values we split the tree at for a given internal node
    split_dims::Vector{Int} # what dimension we split the tree for a given internal node
    n_d::Int # Number of dimension in tree
    last_size::Int # Number of points in last node
    leafsize::Int # Number of points in each node except the last node
    n_leafs::Int # Number of leaf nodes
    n_internal::Int # Number of internal nodes
    cross_node::Int # Index of first node in last row
    first_leaf_row::Int # The row which contains the first leaf node
    offset::Int # How many nodes in first_leaf_row before the first leaf node
    data_reordered::Bool # If the data has been reordered for optimizing memory access.
end


function show(io::IO, tree::KDTree)
    print(io, string("KDTree from ", size(tree.data, 2),
                     " points in ", size(tree.data, 1), " dimensions."))
end


# Helper functions to get node numbers and points
getleft(idx::Int) = 2idx
getright(idx::Int) = 2idx + 1
getparent(idx::Int) = div(idx, 2)
isleaf(tree::KDTree, idx::Int) = idx > tree.n_internal


# Returns the number of points in a node.
# All nodes have leafsize points in them except possibly the last node.
function n_ps(n_leafs::Int, n_internal::Int, leafsize::Int,
              last_size::Int, idx::Int)
    if idx != n_leafs + n_internal
        return leafsize
    else
        return last_size
    end
end


# Returns the index for the first point for a given leaf node.
function point_index(cross_node::Int, offset::Int, last_size:: Int,
                     leafsize::Int, n_internal::Int, idx::Int)
    if idx >= cross_node
        return (idx - cross_node) * leafsize + 1
    else
        return ((offset + idx - 1 - n_internal) * leafsize
                  + last_size + 1)
    end
end

# Returns the range of indices for any node, internal or terminal.
# For an internal node the range of indices are of all the points
# that that subtree contains.
# This range might be split into two subranges.
# The first two return values are the indices for the first range and the
# second two are the indices for the second range.
function node_indices(tree, index)
    if isleaf(tree, index)
        p_index = point_index(tree.cross_node, tree.offset, tree.last_size,
                              tree.leafsize, tree.n_internal, index)
        n_p =  n_ps(tree.n_leafs, tree.n_internal, tree.leafsize,
                    tree.last_size, index)
        l1 = p_index
        r1 = p_index + n_p -1
        return l1, r1, 0, -1
    else
        first_leaf = tree.n_internal + 1
        left, right = index, index
        # Find most left node
        while (left < first_leaf) left = getleft(left) end

        # Find most right node
        while (right < first_leaf) right = getright(right) end

        if left < right
            return node_indices(tree, left)[1], node_indices(tree, right)[2], 0, -1
        else # here we need to use the subranges
            l1 = node_indices(tree, first_leaf)[1]
            r1 = node_indices(tree, right)[2]
            l2 = node_indices(tree, left)[1]
            r2 = node_indices(tree, tree.n_leafs + tree.n_internal)[2]
            return l1, r1, l2, r2
        end
    end
end

# Constructor for KDTree
function KDTree{T <: AbstractFloat}(data::Matrix{T};
                                    leafsize::Int = 10,
                                    reorder::Bool = true)

    if size(data, 2) == 0
        error("Need at least 1 point")
    end

    if size(data, 1) == 0
        error("Need at least 1 dimensional points")
    end

    n_d, n_p = size(data)

    n_leaf =  ceil(Integer, n_p / leafsize)
    n_internal = n_leaf - 1

    # Row whith first leaf node
    l = floor(Integer, log2(n_leaf))

    # Node index of first leaf node in last row
    cross_node = @compat Int(2^(l+1))

    # This only happens when n_p / leafsize is a power of 2?
    if cross_node >= n_internal + n_leaf
        cross_node = div(cross_node, 2)
    end

    # Number of nodes into row l where leaf nodes start
    offset = 2(n_leaf - 2^l) - 1

    # Last node might not have leafsize points in it.
    last_size = n_p % leafsize
    if last_size == 0
        last_size = leafsize
    end

    indices = collect(1:n_p)

    if reorder
        indices_reordered = Array(Int, n_p)
        data_reordered = Array(T, n_d, n_p)
    else
        # Dummy variables here
        indices_reordered = Array(Int, 0)
        data_reordered = Array(T, n_d, 0)
    end

    split_vals = Array(T, n_internal)
    split_dims = Array(Int, n_internal)
    hyper_recs = Array(HyperRectangle{T}, n_internal + n_leaf)

    # Create first bounding hyper rectangle
    maxes = Array(T, n_d)
    mins = Array(T, n_d)
    for j in 1:n_d
        dim_max = typemin(T)
        dim_min = typemax(T)
        for k in 1:n_p
            dim_max = max(data[j, k], dim_max)
            dim_min = min(data[j, k], dim_min)
        end
        maxes[j] = dim_max
        mins[j] = dim_min
    end

    hyper_recs[1] = HyperRectangle{T}(mins, maxes)

    low = 1
    high = n_p

    # Call the recursive KDTree builder
    build_KDTree(1, data, data_reordered, split_vals, split_dims, hyper_recs,
                 indices, indices_reordered, leafsize, low, high, cross_node,
                 offset, last_size, n_internal, reorder)

    if reorder
        KDTree(data_reordered, hyper_recs, indices_reordered, split_vals, split_dims, n_d, last_size,
              leafsize, n_leaf, n_internal, cross_node, l ,offset, reorder)
    else
        KDTree(data, hyper_recs, indices, split_vals, split_dims, n_d, last_size,
              leafsize, n_leaf, n_internal, cross_node, l ,offset, reorder)
    end


end


# Recursive function to build the tree.
# Calculates what dimension has the maximum spread,
# and how many points to send to each side.
# Sorts the indices at the split_
# Splits the hyper cubes and calls recursively
# with the new cubes and node indices.
# TODO: The number of arguments are growing ridiculous.
function build_KDTree{T <: AbstractFloat}(index::Int,
                                          data::Matrix{T},
                                          data_reordered::Matrix{T},
                                          split_vals::Vector{T},
                                          split_dims::Vector{Int},
                                          hyper_recs::Vector{HyperRectangle{T}},
                                          indices::Vector{Int},
                                          indices_reorder::Vector{Int},
                                          leafsize::Int,
                                          low::Int,
                                          high::Int,
                                          cross_node::Int,
                                          offset::Int,
                                          last_size::Int,
                                          n_internal::Int,
                                          reorder::Bool)


    n_p = high - low + 1 # Points left in this subtree

    if n_p <= leafsize

        if reorder
            # Here we reorder the data points so that points contained
            # in nodes with an index close to each other are also themselves
            # close in memory.
            p_index = point_index(cross_node, offset, last_size,
                                  leafsize, n_internal, index)
            n_p_node = n_ps(n_internal + 1, n_internal, leafsize,
                            last_size, index)

            for i in p_index:(p_index + n_p_node - 1)
                idx = indices[i]
                @devec data_reordered[:, i] = data[:, idx]

                # Saves the inverse n
                indices_reorder[i] = idx
            end
        end
        return
    end

    # The number of leafs left, ceil to count a partially filled node as 1.
    n_leafs = ceil(Integer, n_p / leafsize)

    # Rows left in the sub tree
    k = floor(Integer, log2(n_leafs))

    # Number of leftover nodes needed
    rest = n_leafs - 2^k

    # If the last leaf node will be on the right side of the tree we
    # send points so that left tree will be perfectly filled,
    # else we do the opposite.
    if k == 0
        mid_idx = low
    elseif n_p <= 2 * leafsize
        mid_idx = leafsize + low
    elseif  rest > 2^(k-1) # Last node over the "half line" in the row
        mid_idx = 2^k * leafsize + low
    elseif rest == 0 # Perfectly filling the last row
        mid_idx = 2^(k-1)* leafsize + low
    else
        mid_idx = n_p - 2^(k-1) * leafsize + low
    end


    n_d = size(data, 1)
    split_dim = 1
    max_spread = zero(T)
    hyper_rec = hyper_recs[index]
    # Find dimension and and spread where the spread is maximal
    for d in 1:n_d
        spread = hyper_rec.maxes[d] - hyper_rec.mins[d]
        if spread > max_spread
            max_spread = spread
            split_dim = d
        end
    end

    split_dims[index] = split_dim

    # sorts the points in the maximum spread dimension s.t
    # data[split_dim, a]) < data[split_dim, b]) for all a <= mid_idx, b > mid_idx
    _select!(indices, mid_idx, low, high, data, split_dim)

    split_val = data[split_dim, indices[mid_idx]]
    split_vals[index] = split_val

    # Create the hyper rectangles for the children
    hyper_rec_1, hyper_rec_2 = split_hyper_rec(hyper_recs[index], split_dim, split_val)
    hyper_recs[getleft(index)] = hyper_rec_1
    hyper_recs[getright(index)] = hyper_rec_2

    # Recursive call left and right sub tree
    build_KDTree(getleft(index), data, data_reordered,
                 split_vals, split_dims, hyper_recs,
                 indices, indices_reorder, leafsize, low, mid_idx - 1,
                 cross_node, offset, last_size,
                 n_internal, reorder)

    build_KDTree(getright(index), data, data_reordered,
                 split_vals, split_dims, hyper_recs,
                 indices, indices_reorder, leafsize, mid_idx, high,
                 cross_node, offset, last_size,
                 n_internal, reorder)
end

####################################################################
# Distances
####################################################################
# Reduced euclidian distances
@inline function euclidean_distance_red{T <: AbstractFloat}(point_1::AbstractVector{T},
                                                            point_2::AbstractVector{T})
    dist = zero(T)
    @simd for i in eachindex(point_1)
        @inbounds dist += abs2(point_1[i] - point_2[i])
    end
    return dist
end


@inline function euclidean_distance_red{T <: AbstractFloat}(data::Matrix{T},
                                                            idx::Int,
                                                            point::AbstractVector{T})
    dist = zero(T)
    @simd for i in eachindex(point)
        @inbounds dist += abs2(data[i, idx] - point[i])
    end
    return dist
end

@inline function euclidean_distance_red{T <: AbstractFloat}(tree::KDTree{T},
                                                            idx::Int,
                                                            tree2::KDTree{T},
                                                            idx2::Int)
    dist = 0.0
    for i = 1:tree.n_d
        @inbounds dist += abs2(tree.data[i, idx] - tree2.data[i, idx2])
    end
    return dist
end


####################################################################
# Query functions
####################################################################
# Finds the k nearest neighbour to a given point in space..
# The full_rec_dim argument requires some explanation.
# There is a tradeoff between improving the culling in the tree and
# the extra time that improvement take compared to just checking some
# extra nodes. If full_rec_dim > the dimensions in the tree the actual
# minimal distance between the rectangle and the point will be used for
# culling. Else, only the distance between the point and the rectangle
# for the dimension the rectangle was split at is used. As the dimension
# is increasing it is more and more important to use the actual distance.
# The tipping point seems to be right now for around a dimension of 6
# which is the value that will be used default.
function knn{T <: AbstractFloat}(tree::KDTree, point::Vector{T}, k::Int, full_rec_dim::Int = 6)

    # Check that k is not greater than points in tree
    if k > size(tree.data, 2) || k <= 0
        error("k > number of points in tree or <= 0")
    end

    # Check consistent dimension between tree and point
    if size(point,1) != size(tree.data, 1)
        error(string("Wrong dimension of input point, points in the tree",
                     " have dimension ", size(tree.data, 1), " you",
                     " gave a point with dimension ", size(point,1), "."))
    end

    # Initiate heaps to store indices and distances
    best_idxs = [-1 for _ in 1:k]
    best_dists = [typemax(T) for _ in 1:k]

    if tree.n_d < full_rec_dim
        _knn_small(tree, point, k, best_idxs, best_dists, 1)
    else
        init_min = get_min_distance(tree.hyper_recs[1], point)
        _knn(tree, point, k, best_idxs, best_dists, 1, init_min)
    end

    # Sqrt here because distances are stored in reduced format.
    @inbounds @simd for i in eachindex(best_dists)
        best_dists[i] = sqrt(best_dists[i])
    end

    # Translate indices back to the original if we reoredered the data
    if tree.data_reordered
        for j in 1:k
            @inbounds best_idxs[j] = tree.indices[best_idxs[j]]
        end
    end

    # Sort both heaps according to the distance.
    heap_sort_inplace!(best_dists, best_idxs)

    return best_idxs, best_dists
end

# Convert
knn{P <: Real}(tree::KDTree, point::Vector{P}, k::Int) =
  knn(tree, float(point), k)


function _knn{T <: AbstractFloat}(tree::KDTree{T},
                                  point::Vector{T},
                                  k::Int,
                                  best_idxs ::Vector{Int},
                                  best_dists::Vector{T},
                                  index::Int,
                                  min_dist::T)

    data_reordered = tree.data_reordered
    indices = tree.indices


    # If leaf, brute force through the points in the node and
    # if the distance is smaller add both the distance and index
    # to their respective heaps.
    if isleaf(tree, index)
        p_index = point_index(tree.cross_node, tree.offset, tree.last_size,
                              tree.leafsize, tree.n_internal, index)
        n_p =  n_ps(tree.n_leafs, tree.n_internal, tree.leafsize,
                    tree.last_size, index)
        @inbounds for z in p_index:p_index + n_p - 1
            idx = data_reordered ? z : indices[z]
            dist_d = euclidean_distance_red(tree.data, idx, point)
            if dist_d <= best_dists[1]
                best_dists[1] = dist_d
                best_idxs[1] = idx
                percolate_down!(best_dists, best_idxs, dist_d, idx)
            end
        end
        return
    end

    # Find what subtree is closest to the point
    if point[tree.split_dims[index]] < tree.split_vals[index]
        close_node = getleft(index)
        far_node = getright(index)
    else
        far_node = getleft(index)
        close_node = getright(index)
    end

    # Call close subtree first to improve culling
    _knn(tree, point, k, best_idxs, best_dists, close_node, min_dist)

    rec = tree.hyper_recs[close_node]
    dim = tree.split_dims[index]
    min_d_dim = get_min_dim(rec, point, dim)

    far_min = min_dist - min_d_dim + abs2(point[tree.split_dims[index]] - tree.split_vals[index])

    # Only call the sub tree further away if it is close enough
    if far_min < best_dists[1]
            _knn(tree, point, k, best_idxs, best_dists, far_node, far_min)
    end
    return
end

function _knn_small{T <: AbstractFloat}(tree::KDTree{T},
                                        point::Vector{T},
                                        k::Int,
                                        best_idxs ::Vector{Int},
                                        best_dists::Vector{T},
                                        index::Int)

    data_reordered = tree.data_reordered
    indices = tree.indices

    # If leaf, brute force through the points in the node and
    # if the distance is smaller add both the distance and index
    # to their respective heaps.
    if isleaf(tree, index)
        p_index = point_index(tree.cross_node, tree.offset, tree.last_size,
                              tree.leafsize, tree.n_internal, index)
        n_p =  n_ps(tree.n_leafs, tree.n_internal, tree.leafsize,
                    tree.last_size, index)
        @inbounds for z in p_index:p_index + n_p - 1
            idx = data_reordered ? z : indices[z]
            dist_d = euclidean_distance_red(tree.data, idx, point)
            if dist_d <= best_dists[1]
                best_dists[1] = dist_d
                best_idxs[1] = idx
                percolate_down!(best_dists, best_idxs, dist_d, idx)
            end
        end
        return
    end

    # Find what subtree is closest to the point
    if point[tree.split_dims[index]] < tree.split_vals[index]
        close_node = getleft(index)
        far_node = getright(index)
    else
        far_node = getleft(index)
        close_node = getright(index)
    end

    # Call close subtree first to improve culling
    _knn_small(tree, point, k, best_idxs, best_dists, close_node)

    # Only call the sub tree further away if it is close enough
    if abs2(point[tree.split_dims[index]] - tree.split_vals[index]) < best_dists[1]
            _knn_small(tree, point, k, best_idxs, best_dists, far_node)
    end
    return
end


# Returns the list of indices closer than a given point at a given
# radius.
function inball{T <: AbstractFloat}(tree::KDTree{T},
                                              p::Vector{T},
                                              radius::T,
                                              sort::Bool = false)

    # Check dimension consistency
    if size(p, 1) != size(tree.data, 1)
        error(string("Wrong dimension of input point, points in the tree",
                     " have dimension ", size(tree.data, 1), " you",
                     " gave a point with dimension ", size(p,1), "."))
    end


    index = 1
    idx_in_ball = Int[]

    init_min, init_max = get_min_max_distance(tree.hyper_recs[1], p)
    _in_ball(tree, index, p, abs2(radius) , idx_in_ball, init_min, init_max)

    if tree.data_reordered
        for i in 1:length(idx_in_ball)
            idx_in_ball[i] = tree.indices[idx_in_ball[i]]
        end
    end

    if sort
        sort!(idx_in_ball)
    end

    return idx_in_ball
end

inball{P <: Real, R <: Real}(tree::KDTree, p::Vector{P}, r::R) =
  inball(tree, float(p), float(r))


# Explicitly check the distance between leaf node and point while traversing
function _in_ball{T <: AbstractFloat}(tree::KDTree{T},
                                            index::Int,
                                            p::Vector{T},
                                            r::T,
                                            idx_in_ball::Vector{Int},
                                            min_dist::T,
                                            max_dist::T)

    if min_dist > r # Hyper sphere is outside hyper rectangle, skip the whole sub tree
        return
    end

    if max_dist < r
        addall(tree, index, idx_in_ball)
        return
    end

    if isleaf(tree, index)
        p_index = point_index(tree.cross_node, tree.offset, tree.last_size,
                              tree.leafsize, tree.n_internal, index)
        n_p =  n_ps(tree.n_leafs, tree.n_internal, tree.leafsize,
                    tree.last_size, index)
        for z in p_index:p_index + n_p - 1
            idx = tree.data_reordered ? z : tree.indices[z]
            dist_d = euclidean_distance_red(tree.data, idx, p)
            if dist_d < r
                push!(idx_in_ball, idx)
            end
        end
        return
    end

    d =  tree.split_dims[index]

    # We are querying 3 times so makes no sense to do this optimization unless n_d > 3.
    if tree.n_d > 3
        # # Remove contribution from this rectangle
        min_d_dim, max_d_dim = get_min_max_dim(tree.hyper_recs[index], p, d)
        max_dist -= max_d_dim
        min_dist -= min_d_dim

        # Add contribution of left node
        min_d_dim, max_d_dim = get_min_max_dim(tree.hyper_recs[getleft(index)], p, d)
        _in_ball(tree, getleft(index), p, r, idx_in_ball, min_dist + min_d_dim, max_dist + max_d_dim)

        # Add contribution of right node
        min_d_dim, max_d_dim = get_min_max_dim(tree.hyper_recs[getright(index)], p, d)
        _in_ball(tree, getright(index), p, r, idx_in_ball, min_dist + min_d_dim, max_dist + max_d_dim)
    else
        min_dist, max_dist = get_min_max_distance(tree.hyper_recs[index], p)
        _in_ball(tree, getleft(index), p, r, idx_in_ball, min_dist, max_dist)
        _in_ball(tree, getright(index), p, r, idx_in_ball, min_dist, max_dist)
    end

end


# Adds everything in this subtree since we have determined
# that the hyper rectangle completely encloses the hyper sphere
function addall(tree::KDTree, index::Int, idx_in_ball::Vector{Int})
    l1, r1, l2, r2 =  node_indices(tree, index)
    for r in (l1:r1, l2:r2)
        for z in r
            idx = tree.data_reordered ? z : tree.indices[z]
            push!(idx_in_ball, idx)
        end
    end
    return
end



# Returns the sorted list of indices for all points in the tree inside a
# hypersphere of a given point with a given radius.
function inball{T <: AbstractFloat}(tree::KDTree{T},
                                    tree2::KDTree{T},
                                    r::T,
                                    sort::Bool = false)

    if tree.data_reordered || tree2.data_reordered
        error(string("Trees with reoredered data is not currently supported.",
              " Please create your trees using for example KDTree(data, 10, true)"))
    end

    # Check tree dimensions consistency
    if size(tree.data, 1) != size(tree2.data, 1)
        error(string("Wrong dimension of input tree, points in the tree",
                     " have dimension ", size(tree.data, 1), " you",
                     " gave a tree with dimension ", size(tree2.data, 1), "."))
    end

    index = 1
    idx_in_ball =  Array((Array{Int64, 1}), size(tree.data, 2))
    for i in 1:length(idx_in_ball)
        idx_in_ball[i] = Int[]
    end

    init_min, init_max = get_min_max_distance(tree.hyper_recs[1], tree2.hyper_recs[1])

    # Square r since we work in reduced distance
    _in_ball(tree, index, tree2, index, abs2(r) , idx_in_ball, init_min, init_max)

    if sort
        for i in 1:length(idx_in_ball)
            sort!(idx_in_ball[i])
        end
    end

    return idx_in_ball
end

inball{P <: Real, R <: Real}(tree::KDTree, p::Vector{P}, r::R) =
  inball(tree, float(p), float(r))


# Explicitly check the distance between leaf node and point while traversing
function _in_ball{T <: AbstractFloat}(tree::KDTree{T},
                                      index::Int,
                                      tree2::KDTree{T},
                                      index2::Int,
                                      r::T,
                                      idx_in_ball::Vector{Vector{Int}},
                                      min_dist::T,
                                      max_dist::T)

    # This is currently ineffective for large d since we recalculate the distance for all dimensions
    # We should instead just calculate the change in distance for the current split_dim
    min_dist, max_dist = get_min_max_distance(tree.hyper_recs[index], tree2.hyper_recs[index2])
    if min_dist > r # Hyper sphere is outside hyper rectangle, skip the whole sub tree
        return

    elseif max_dist < r
        addall(tree, index, tree2, index2, idx_in_ball)

    elseif isleaf(tree, index)
        if isleaf(tree2, index2)
            l1, r1, l2, r2 =  node_indices(tree, index)
            l1_2, r1_2, l2_2, r2_2 =  node_indices(tree2, index2)
            for i in l1:r1
                idx = tree.data_reordered ? i : tree.indices[i]
                for j in l1_2:r1_2
                    idx2 = tree2.data_reordered ? j : tree2.indices[j]
                    dist_d = euclidean_distance_red(tree, idx, tree2, idx2)
                    if dist_d < r
                        push!(idx_in_ball[idx], idx2)
                    end
                end
            end
        else # tree: leaf, tree2: internal
            _in_ball(tree, index, tree2, getleft(index2), r, idx_in_ball, min_dist, max_dist)
            _in_ball(tree, index, tree2, getright(index2), r, idx_in_ball, min_dist, max_dist)
        end
    else # tree: internal
        if isleaf(tree2, index2) # tree2: leaf
            _in_ball(tree, getleft(index), tree2, index2, r, idx_in_ball, min_dist, max_dist)
            _in_ball(tree, getright(index), tree2, index2, r, idx_in_ball, min_dist, max_dist)
        else # both internal
            _in_ball(tree, getleft(index), tree2, getleft(index2), r, idx_in_ball, min_dist, max_dist)
            _in_ball(tree, getright(index), tree2, getright(index2), r, idx_in_ball, min_dist, max_dist)
            _in_ball(tree, getleft(index), tree2, getright(index2), r, idx_in_ball, min_dist, max_dist)
            _in_ball(tree, getright(index), tree2, getleft(index2), r, idx_in_ball, min_dist, max_dist)
        end
    end
    return
end


# Adds everything in this subtree since we have determined
# that the hyper rectangle completely encloses the hyper sphere
function addall(tree::KDTree, index::Int, tree2::KDTree, index2::Int, idx_in_ball::Vector{Vector{Int}})
    l1, r1, l2, r2 = node_indices(tree, index)
    l1_2, r1_2, l2_2, r2_2 = node_indices(tree2, index2)
    for r in (l1:r1, l2:r2)
        for i in r
            idx = tree.data_reordered ? i : tree.indices[i]
            for r2 in (l1_2:r1_2, l2_2:r2_2)
                for j in r2
                    idx2 = tree2.data_reordered ? j : tree2.indices[j]
                    push!(idx_in_ball[idx], idx2)
                end
            end
        end
    end
    return
end


####################################################################
# Taken from https://github.com/JuliaLang/julia/blob/v0.3.5/base/sort.jl
# and modified because I couldn't figure out how to get rid of
# the overhead when I passed in a new anonymous function
# to the "by" argument in each node. I also removed the return value.
function _select!{T <: AbstractFloat}(v::AbstractVector, k::Int, lo::Int,
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


# In place heap sort
function heap_sort_inplace!{T <: AbstractFloat}(xs::AbstractArray{T}, xis::AbstractArray{Int})
    @inbounds for i in length(xs):-1:2
        xs[i], xs[1] = xs[1], xs[i]
        xis[i], xis[1] = xis[1], xis[i]
        percolate_down!(xs, xis, xs[1], xis[1], i-1)
    end
end

# Binary min-heap percolate down.
function percolate_down!{T <: AbstractFloat}(xs::AbstractArray{T},
                                             xis::AbstractArray{Int},
                                             dist::T,
                                             index::Int,
                                             len::Int=length(xs))
    i = 1
    @inbounds while (l = getleft(i)) <= len
        r = getright(i)
        j = ifelse(r > len || (xs[l] > xs[r]), l, r)
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
