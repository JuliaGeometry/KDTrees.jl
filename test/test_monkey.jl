# This contains a bunch of random tests that should hopefully detect if
# some edge case has been missed in the real tests

facts("KDTrees.monkey") do

context("KDTrees.monkey.knn") do

# Checks that we find existing point in the tree
# and that they are closest
for i in 1:10
dim_data = 3
size_data = 1000
data = rand(dim_data, size_data)

    for j = 1:15
        tree = KDTree(data, leafsize = rand(1:15), reorder = rand(Bool))
        n = rand(1:size_data)
        idx, dist = knn(tree, data[:,n], rand(1:30))
        @fact issorted(dist) --> true
        @fact n --> idx[1]

    end
end

# Compares vs Brute Force
for i in 1:10
    if Base.VERSION >= v"0.4.0-dev"
        pq = PriorityQueue(Int, Float64, Base.Order.Reverse)
    else
        pq = PriorityQueue{Int, Float64}(Base.Order.Reverse)
    end


    k = rand(1:12)
    for i in 1:k
        enqueue!(pq, -i, Inf)
    end

    dim_data = rand(1:5)
    size_data = rand(100:151)
    data = rand(dim_data, size_data)
    tree = KDTree(data, leafsize = rand(1:10))
    p = rand(dim_data)

    # Brute force
    for n in 1:size_data
        d = sqrt(KDTrees.euclidean_distance_red(data[:,n], p))
        if d <= peek(pq)[2] # Closer than the currently k closest.
            dequeue!(pq)
            enqueue!(pq, n, d)
        end
    end

    idx, dist = knn(tree, p, k)

    for i in 1:length(idx)
        @fact idx[i] in keys(pq) --> true
    end
end
end # context

context("KDTrees.monkey.inball") do

# Test against brute force
for rn in [true, false]
    for i in 1:10
        idx = Int[]
        dim_data = rand(1:6)
        size_data = rand(20:250)
        data = rand(dim_data, size_data)
        tree = KDTree(data, leafsize = 3, reorder = rn)
        p = zeros(dim_data)
        r = 0.3
        # Brute force
        for n in 1:size_data
            d = sqrt(KDTrees.euclidean_distance_red(data[:,n], p))
            if d <= r # Closer than the currently k closest.
                push!(idx, n)
            end
        end

        q_idxs = inball(tree, p, r)

        @fact allin(idx, q_idxs) --> true
    end
end

end # context


context("KDTrees.monkey.coupled") do
# Tests that the n-points in a random hyper sphere around
# a random point are all the n-closest points to that point.... yeah
for rn in [true, false]

    for i in 1:10
        dim_data = rand(1:5)
        size_data = rand(100:1000)
        data = randn(dim_data, size_data)
        tree = KDTree(data, leafsize = rand(1:10), reorder = rn)

        point = [randn() for x in 1:dim_data]
        idxs_ball = []
        r = 0.1
        while length(idxs_ball) < 10
            r *= 2.0
            idxs_ball = inball(tree,  point, r)
        end
        idxs_knn, dists = knn(tree, point, length(idxs_ball))

        @fact allin(idxs_knn, idxs_ball) --> true
    end
end
end # context

end # facts






