using KDTrees

function run_bench_query_ball()
    dims = 3
    n_points = [10^i for i in 3:5]
    rs = [0.1 0.2 0.3 0.4]

    times = fill(0.0, length(rs), length(n_points))

    n_iters = 2000

    # Compile it
    tree = KDTree(randn(2,2))
    inball(tree, zeros(2), 0.1)

    for (i, r) in enumerate(rs)
        for (j , n_point) in enumerate(n_points)
            data = rand(dims, n_point)
            tree = KDTree(data, 10, true)
            t = time_ns()
            for z = 1:n_iters
                p = rand(dims)
                idxs = inball(tree, p, r)
            end
            t = (time_ns() - float(t)) / 10^9
            times[i,j] =  t / n_iters
        end
    end

    println(times)
end

run_bench_query_ball()

# Removed old benchmarks because apparently the sort in the
# end was taking most of the time making the results useless.