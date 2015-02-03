using KDtree

dims = 3
n_points = [10^i for i in 3:5]
rs = [0.1 0.2 0.3 0.4]

times = Array(Float64, length(rs), length(n_points))

# Compile it
tree = KDTree(randn(2,2))
query_ball_point(tree, zeros(2), 0.1)


for (i, r) in enumerate(rs)
    for (j , n_point) in enumerate(n_points)
        data = randn(dims, n_point)
        tree = KDTree(data)
        times[i,j]  = @elapsed query_ball_point(tree, zeros(dims), r)
    end
end

println(times)

#=
2015-02-03: ArrayViews + no sqrt
[1.1196e-5 1.9593e-5 7.5885e-5
 1.5239e-5 4.1986e-5 0.000167011
 1.9905e-5 7.8374e-5 0.00041986
 2.6125e-5 0.000124713 0.000603665]


2015-02-03:
[2.1149e-5 3.2966e-5 8.3661e-5
 2.146e-5 5.1628e-5 0.000229523
 2.7369e-5 0.000116005 0.000453138
 4.012e-5 0.000204643 0.000789959]
 =#
