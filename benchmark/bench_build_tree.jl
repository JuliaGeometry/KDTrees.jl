using KDtree

dims = 1:6
n_points = [10^i for i in 1:5]

times = Array(Float64, length(dims), length(n_points))

# Compile
KDTree(randn(2,2))

i = 0
for dim in dims
    j = 0
    i+=1
    for n_point in n_points
        j+=1
        data = randn(dim, n_point)
        times[i,j]  = @elapsed KDTree(data)
    end
end

println(times)

#=
times 2015-02-01:

[4.6477e-5 0.000318537 0.004562246 0.054535408 0.796257552
 2.8026e-5 0.000350114 0.00443035 0.076212208 0.952349006
 2.895e-5 0.000329126 0.004856197 0.059021987 0.916525079
 3.1511e-5 0.000360546 0.004838845 0.075776018 0.921683179
 2.9422e-5 0.000347066 0.004848977 0.078650667 0.961300694
 3.1327e-5 0.00036603 0.005024072 0.080929123 0.992533382]
 =#
