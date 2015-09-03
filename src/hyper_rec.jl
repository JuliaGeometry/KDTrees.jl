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


# Splits a hyper rectangle into two rectangles by dividing the
# rectangle at a specific value in a given dimension.
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

############################################
# Rectangle - Point functions
############################################
@inline function get_min_dim{T <: FloatingPoint}(rec::HyperRectangle{T}, point::Vector{T}, dim::Int)
    @inbounds d = abs2(max(0, max(rec.mins[dim] - point[dim], point[dim] - rec.maxes[dim])))
    d
end

@inline function get_max_dim{T <: FloatingPoint}(rec::HyperRectangle{T}, point::Vector{T}, dim::Int)
    @inbounds d = abs2(max(rec.maxes[dim] - point[dim], point[dim] - rec.mins[dim]))
    d
end

# Max distance between rectangle and point
@inline function get_max_distance{T <: FloatingPoint}(rec::HyperRectangle{T}, point::Vector{T})
    max_dist = zero(T)
    @inbounds for dim in 1:size(point,1)
        max_dist += get_max_dim(rec, point, dim)
    end
    return max_dist
end

# Min distance between rectangle and point
@inline function get_min_distance{T <: FloatingPoint}(rec::HyperRectangle{T}, point::Vector{T})
    min_dist = zero(T)
    @inbounds for dim in 1:size(point,1)
        min_dist += get_min_dim(rec, point, dim)
    end
    return min_dist
end

# (Min, Max) distance between rectangle and point
@inline function get_min_max_distance{T <: FloatingPoint}(rec::HyperRectangle{T}, point::Vector{T})
    min_dist = get_min_distance(rec, point)
    max_dist = get_max_distance(rec, point)
    return min_dist, max_dist
end

# (Min, Max) distance between rectangle and point for a certain dim
@inline function get_min_max_dim{T <: FloatingPoint}(rec::HyperRectangle{T}, point::Vector{T}, dim::Int)
    min_dist_dim = get_min_dim(rec, point, dim)
    max_dist_dim = get_max_dim(rec, point, dim)
    return min_dist_dim, max_dist_dim
end


############################################
# Rectangle-rectangle functions
############################################
@inline function get_min_dim{T <: FloatingPoint}(rec1::HyperRectangle{T},  rec2::HyperRectangle{T}, dim::Int)
    return abs2(max(0, max(rec1.mins[dim] - rec2.maxes[dim], rec2.mins[dim] - rec1.maxes[dim])))
end

@inline function get_max_dim{T <: FloatingPoint}(rec1::HyperRectangle{T},  rec2::HyperRectangle{T}, dim::Int)
    return abs2(max(rec1.maxes[dim] - rec2.mins[dim],  rec2.maxes[dim] - rec1.mins[dim]))
end

@inline function get_min_max_distance{T <: FloatingPoint}(rec1::HyperRectangle{T}, rec2::HyperRectangle{T})
    min_dist = get_min_distance(rec1, rec2)
    max_dist = get_max_distance(rec1, rec2)
    return min_dist, max_dist
end

@inline function get_min_distance{T <: FloatingPoint}(rec1::HyperRectangle{T}, rec2::HyperRectangle{T})
    min_dist = zero(T)
    @inbounds for dim in 1:length(rec1.maxes)
        min_dist += get_min_dim(rec1, rec2, dim)
    end
    return min_dist
end

@inline function get_max_distance{T <: FloatingPoint}(rec1::HyperRectangle{T}, rec2::HyperRectangle{T})
    max_dist = zero(T)
    @inbounds for dim in 1:length(rec1.maxes)
        max_dist += get_max_dim(rec1, rec2, dim)
    end
    return max_dist
end

@inline function get_min_max_dim{T <: FloatingPoint}(rec1::HyperRectangle{T},  rec2::HyperRectangle{T}, dim::Int)
    min_dist_dim = get_min_distance(rec1, rec2, dim)
    max_dist_dim = get_max_distance(rec1, rec2, dim)
    return min_dist_dim, max_dist_dim
end
