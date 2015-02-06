using KDtree


function run_bench_knn()
    dims = 3
    n_points = [10^i for i in 3:6]
    ks = [1, 3, 10, 50, 100, 500]

    n_iters = 100

    times = fill(0.0, length(ks), length(n_points))

    # Compile it
    tree = KDTree(randn(2,2))
    k_nearest_neighbour(tree, zeros(2), 1)

    for (i, k) in enumerate(ks)
        for (j , n_point) in enumerate(n_points)
            data = rand(dims, n_point)
            tree = KDTree(data)
            for z in 1:n_iters
            	p = rand(dims)
                times[i,j]  += @elapsed k_nearest_neighbour(tree, p, k)
            end
            times[i,j] /= n_iters
        end
    end

    println(times)
    return
end

run_bench_knn()

#=
2015-02-06 (removed old bench since they used a too volatile method)
[2.54688e-6 7.427709999999996e-6 4.546329999999998e-6 4.439879999999999e-6
 3.1702599999999993e-6 3.8126800000000002e-6 5.36744e-6 1.047631e-5
 9.24091e-6 2.3145960000000002e-5 2.319162000000002e-5 2.1059090000000006e-5
 4.556972999999999e-5 5.366267e-5 6.215093000000004e-5 7.335329000000001e-5
 4.875141000000001e-5 5.9763709999999965e-5 0.00013036107000000003 0.00018395148999999993
 0.0007265081999999999 0.0011446717200000003 0.0012804762699999998 0.00121552762]
 =#
