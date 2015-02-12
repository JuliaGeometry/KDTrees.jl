facts("KDTrees") do

    context("KDTrees.nearest_neighbour") do

        dim_data = 3
        size_data = 1000
        data = rand(dim_data, size_data)

        

        # Checking that we find existing points
        for i = 1:50
            tree = KDTree(data, rand(1:15))
            n = rand(1:size_data)
            idx, dist = k_nearest_neighbour(tree, data[:,n], 1)
            @fact n => idx[1]
            @fact KDTrees.euclidean_distance(data[:,idx[1]], data[:, n]) => roughly(0.0)
        end

        # Check results vs brute force
        if Base.VERSION >= v"0.4.0-dev"
            pq = PriorityQueue(Int, Float64, Base.Order.Reverse)
        else
            pq = PriorityQueue{Int, Float64}(Base.Order.Reverse)
        end

        k = 3
        for i in 1:k
            enqueue!(pq, -i, Inf)
        end

        dim_data = 3
        size_data = 500
        data = rand(dim_data, size_data)
        tree = KDTree(data)
        p = rand(dim_data)

        # Brute force
        for n in 1:size_data
            d = sqrt(KDTrees.euclidean_distance(data[:,n], p))
            if d <= peek(pq)[2] # Closer than the currently k closest.
                dequeue!(pq)
                enqueue!(pq, n, d)
            end
        end

        idx, dist = k_nearest_neighbour(tree, p, k)

        for i in 1:length(idx)
            @fact idx[i] in keys(pq) => true
        end

       # 8 node rectangle
        data = [0.0 0.0 0.0 0.5 0.5 1.0 1.0 1.0;
                0.0 0.5 1.0 0.0 1.0 0.0 0.5 1.0]
        tree = KDTree(data)

        idxs, dists = k_nearest_neighbour(tree, [0.8, 0.8], 1)
        @fact idxs[1] => 8 # Should be closest to top right corner
        @fact sqrt(0.2^2 + 0.2^2) => roughly(dists[1])

        idxs, dists = k_nearest_neighbour(tree, [0.1, 0.8], 3)
        @fact idxs => [3, 2, 5]

        @fact_throws k_nearest_neighbour(tree, [0.1, 0.8], 10) # k > n_points

        @fact_throws k_nearest_neighbour(tree, [0.1], 10) # n_dim != trees dim
    end  #context

end # facts
