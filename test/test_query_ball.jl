facts("KDTrees") do

context("KDTrees.ball_query") do

for rn in [true, false]
    data = [0.0 0.0 0.0 0.0 1.0 1.0 1.0 1.0;
            0.0 0.0 1.0 1.0 0.0 0.0 1.0 1.0;
            0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0] # 8 node cube

    tree = KDTree(data, 2, rn)

    idxs = query_ball_point(tree, [1.1, 1.1, 1.1], 0.2)
    @fact idxs => [8] # Only corner 8 at least 0.2 distance away from [1.1, 1.1, 1.1]

    idxs = query_ball_point(tree, [0.0, 0.0, 0.5], 0.6)
    @fact allin(idxs, [1, 2]) => true# Corner 1 and 2 at least 0.6 distance away from [0.0, 0.0, 0.5]

    idxs = query_ball_point(tree, [0, 0, 0], 0.6)
    @fact idxs => [1]

    idxs = query_ball_point(tree, [1//3, 1//3, 1//3], 1)
    @fact allin(idxs, [1, 2, 3, 5]) => true

    idxs = query_ball_point(tree, [0.5, 0.5, 0.5], 0.2)
    @fact idxs => [] #

    idxs = query_ball_point(tree, [0.5, 0.5, 0.5], 1.0)
    @fact allin(idxs, [1, 2, 3, 4, 5, 6, 7, 8]) => true #

    @fact_throws query_ball_point(tree, [0.1], 1.0) # n_dim != trees dim


    ####################################################################
    # Test against brute force
    for i in 1:10
        idx = Int[]
        dim_data = rand(1:6)
        size_data = rand(20:250)
        data = rand(dim_data, size_data)
        tree = KDTree(data, 3, rn)
        p = zeros(dim_data)
        r = 0.3
        # Brute force
        for n in 1:size_data
            d = sqrt(KDTrees.euclidean_distance_red([data[:,n]], p))
            if d <= r # Closer than the currently k closest.
                push!(idx, n)
            end
        end

        q_idxs = query_ball_point(tree, p, r)

        @fact allin(idx, q_idxs) => true
    end
    ####################################################################
end
end # context


context("KDTrees.yolo_testing") do

for rn in [true, false]
    # Tests that the n-points in a random hyper sphere around
    # a random point are all the n-closest points to that point.... yeah
    for i in 1:10
        dim_data = rand(1:5)
        size_data = rand(100:1000)
        data = randn(dim_data, size_data)
        tree = KDTree(data, 6, rn)

        point = [randn() for x in 1:dim_data]
        idxs_ball = []
        r = 0.1
        while length(idxs_ball) < 10
            r *= 2.0
            idxs_ball = query_ball_point(tree,  point, r)
        end
        idxs_knn, dists = k_nearest_neighbour(tree, point, length(idxs_ball))

        @fact allin(idxs_knn, idxs_ball) => true
    end
      
end

end # context

end # facts
