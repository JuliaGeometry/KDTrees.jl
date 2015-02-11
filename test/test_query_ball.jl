facts("KDTrees") do

    context("KDTrees.ball_query") do

        data = [0.0 0.0 0.0 0.0 1.0 1.0 1.0 1.0;
                0.0 0.0 1.0 1.0 0.0 0.0 1.0 1.0;
                0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0] # 8 node cube

        tree = KDTree(data)

        idxs = query_ball_point(tree, [1.1, 1.1, 1.1], 0.2)
        @fact idxs => [8] # Only corner 8 at least 0.2 distance away from [1.1, 1.1, 1.1]

        idxs = query_ball_point(tree, [0.0, 0.0, 0.5], 0.6)
        @fact idxs => [1, 2] # Corner 1 and 2 at least 0.6 distance away from [0.0, 0.0, 0.5]

        idxs = query_ball_point(tree, [0, 0, 0], 0.6)
        @fact idxs => [1]

        idxs = query_ball_point(tree, [1//3, 1//3, 1//3], 1)
        @fact idxs => [1, 2, 3, 5]

        idxs = query_ball_point(tree, [0.5, 0.5, 0.5], 0.2)
        @fact idxs => [] #

        idxs = query_ball_point(tree, [0.5, 0.5, 0.5], 1.0)
        @fact idxs => [1, 2, 3, 4, 5, 6, 7, 8] #

        @fact_throws query_ball_poin(tree, [0.1], 1.0) # n_dim != trees dim


        idx = Int[]
        dim_data = 3
        size_data = 100
        data = rand(dim_data, size_data)
        tree = KDTree(data)
        p = zeros(dim_data)
        r = 0.3
        # Brute force
        for n in 1:size_data
            d = sqrt(KDTrees.euclidean_distance([data[:,n]], p))
            if d <= r # Closer than the currently k closest.
                push!(idx, n)
            end
        end

        q_idxs = query_ball_point(tree, p, r)

        for i in 1:length(idx)
            @fact q_idxs[i] in idx => true
        end

    end #context


    context("KDTrees.yolo_testing") do

        # Tests that the n-points in a random hyper sphere around
        # a random point are all the n-closest points to that point.... yeah
        dim_data = rand(1:5)
        size_data = rand(100:10000)
        data = randn(dim_data, size_data)
        tree = KDTree(data)

        point = [randn() for x in 1:dim_data]
        idxs_ball = []
        r = 0.1
        while length(idxs_ball) < 10
            r *= 2.0
            idxs_ball = query_ball_point(tree,  point, r)
        end
        idxs_knn, dists = k_nearest_neighbour(tree, point, length(idxs_ball))



        for i in 1:length(idxs_ball)
            @fact idxs_ball[i] in idxs_knn => true
        end
    end #context

end # facts
