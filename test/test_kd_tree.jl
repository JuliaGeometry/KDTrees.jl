facts("KDtree") do

    context("KDtree.nearest_neighbour") do

        dim_data = 3
        size_data = 1000
        data = randn(dim_data, size_data )

        tree = KDTree(data)

        # Checking that we find existing points
        for i = 1:20
            n = rand(1:size_data)
            idx, dist = k_nearest_neighbour(tree, data[:,n], 1)
            @fact n => idx[1]
            @fact KDtree.euclidean_distance(data[:,idx[1]], data[:, n]) => roughly(0.0)
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

    context("KDtree.ball_query") do

        data = [0.0 0.0 0.0 0.0 1.0 1.0 1.0 1.0;
                0.0 0.0 1.0 1.0 0.0 0.0 1.0 1.0;
                0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0] # 8 node cube

        tree = KDTree(data)

        idxs = query_ball_point(tree, [1.1, 1.1, 1.1], 0.2)
        @fact idxs => [8] # Only corner 8 at least 0.2 distance away from [1.1, 1.1, 1.1]

        idxs = query_ball_point(tree, [0.0, 0.0, 0.5], 0.6)
        @fact idxs => [1, 2] # Corner 1 and 2 at least 0.6 distance away from [0.0, 0.0, 0.5]

        idxs = query_ball_point(tree, [0.5, 0.5, 0.5], 0.2)
        @fact idxs => [] #

        idxs = query_ball_point(tree, [0.5, 0.5, 0.5], 1.0)
        @fact idxs => [1, 2, 3, 4, 5, 6, 7, 8] #

        @fact_throws query_ball_poin(tree, [0.1], 1.0) # n_dim != trees dim
    end #context

    context("KDtree.yolo_testing") do

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
