using KDTrees

function run_bench_query_ball()
    dims = 3
    n_points = [10^i for i in 3:5]
    rs = [0.1 0.2 0.3 0.4]

    times = fill(0.0, length(rs), length(n_points))

    # Compile it
    n_iters = 2000

    tree = KDTree(randn(2,2))
    query_ball_point(tree, zeros(2), 0.1)

    for (i, r) in enumerate(rs)
        for (j , n_point) in enumerate(n_points)
            data = rand(dims, n_point)
            tree = KDTree(data)
            for z = 1:n_iters
            	p = rand(dims)
            	times[i,j] += @elapsed query_ball_point(tree, p, r)
            end
            times[i,j] /= n_iters
        end
    end

    println(times)
end

run_bench_query_ball()


#=
2015-02-14
[1.512347500000003e-6 9.946277500000007e-6 5.1572081499999986e-5
 3.7438639999999997e-6 2.3865769500000025e-5 0.00023704194749999996
 7.439227999999991e-6 5.898463300000001e-5 0.0006188465289999998
 1.2508987000000021e-5 0.0001077092400000001 0.001196649754000001]

2015-02-06: (removed old inaccurate results)
[6.202539999999983e-6 3.595097900000001e-5 0.00021109113900000007
 1.3630627000000005e-5 7.9238053e-5 0.0006688952190000003
 5.310767400000003e-5 0.0001541270569999998 0.001355720711000002
 3.338056700000001e-5 0.00026903708300000007 0.0022353193569999976]
 =#
