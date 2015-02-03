using KDtree

dims = 3
n_points = [10^i for i in 3:6]
ks = [1, 3, 10, 50, 100, 500]

times = Array(Float64, length(ks), length(n_points))

# Compile it
tree = KDTree(randn(2,2))
k_nearest_neighbour(tree, zeros(2), 1)

for (i, k) in enumerate(ks)
    for (j , n_point) in enumerate(n_points)
        data = randn(dims, n_point)
        tree = KDTree(data)
        times[i,j]  = @elapsed k_nearest_neighbour(tree, zeros(dims), k)
    end
end

println(times)

#=
2015-02-03: ArrayViews:
[1.3996e-5 2.1771e-5 5.1316e-5 4.4474e-5
 2.0837e-5 2.5502e-5 5.9402e-5 6.8732e-5
 3.2034e-5 4.1364e-5 8.7393e-5 8.2106e-5
 8.5838e-5 8.7394e-5 0.000183183 0.000170121
 0.000156437 0.000143063 0.000305409 0.000290792
 0.000695723 0.001162545 0.001638387 0.001731067]

2015-02-02:
[2.3015e-5 2.1771e-5 7.0288e-5 6.0958e-5
 2.7368e-5 3.5143e-5 8.957e-5 7.7752e-5
 6.0647e-5 5.0072e-5 0.000104499 0.000111029
 0.000194691 0.000216772 0.000293591 0.000318472
 0.000385027 0.000359835 0.000650316 0.000582517
 0.003119404 0.005313561 0.006453713 0.006805152]
 =#
