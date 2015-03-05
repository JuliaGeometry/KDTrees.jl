using KDTrees
using StatsBase
using Plotly
function run_bench_knn_points(dim, knn, exps, rounds)
    println("Running KNN benchmark for 10^(", exps, ") points in ", dim, " dimensions with k = ", knn, ".")
    n_points = [10^i for i in exps]

    times = zeros(length(n_points), rounds)
    gc_disable()
    timer = 0.0
    for (i, n_point) in enumerate(n_points)
        for (j , round) in enumerate(1:rounds)
            n_iters = 25
            print("Round ", j, " out of ", rounds, " for ", dim, "x", int(n_point), "...\r")
            data = rand(dim, int(n_point))
            tree = KDTree(data, 10)
            while true
                gc_disable()
                timer = time_ns()
                for k in 1:n_iters
                	  p = rand(dim)
                    knn(tree, p, knn)
                end
                timer = (time_ns() - float(timer)) / 10^9 # To seconds
                gc_enable()
                if timer < 1.0
                    n_iters *= 3
                    continue
                end
                break # Ends this round
            end
        times[i, j] = timer / n_iters
        end
        print("\n")
    end
    println("\nDone!")
    return mean_and_std(1./times, 2)
end

data = Dict[]
for knn in [1, 5, 10]
  dim = 3
  exps = 3:0.5:6
  rounds = 5

  speeds, stds = run_bench_knn_points(dim, knn, exps, rounds)
  speeds = vec(speeds)
  stds = vec(stds)

  stderr = stds / sqrt(length(stds))
  ymins = speeds - 1.96*stderr
  ymaxs = speeds + 1.96*stderr
  sizes = [10^i for i in exps]
  n_points = [10^i for i in exps]
  trace = [
    "x" => sizes,
    "y" => speeds,
    "type" => "scatter",
    "mode" => "lines+markers",
    "name" => string("k = ", knn),
    "error_y" => [
        "type" => "data",
        "array" => stderr*1.96*2,
        "visible" => true
    ]
  ]
  push!(data, trace)
end

####################################################################
# Plotting
####################################################################

layout = [
  "title" => "KNN search speed (dim = 3)",
  "xaxis" => [
      "title" => "Number of points",
      "type" => "log"
  ],
  "yaxis" => [
      "title" => "KNN searches / second [1/s]",
      "type" => "log"
  ],
  "type" => "log",
  "autorange" => true
]

Plotly.signin("kcarlsson89", "lolololoololololo")

response = Plotly.plot(data, ["layout" => layout,
                              "filename" => "bench_linux_knn_3",
                              "fileopt" => "overwrite"])
plot_url = response["url"]


#=
p = plot(x = sizes, y = speeds,  Geom.point, Geom.line, Geom.errorbar, Scale.x_log10, Scale.y_log10,
         ymin = ymins,
         ymax = ymaxs,
         Guide.xticks(ticks=collect(exps)),
         Guide.yticks(ticks=collect(log10(speeds))),
         Guide.xlabel("Number of points"),
         Guide.ylabel("KNN searches / second [1/s]"),
         Guide.title("KNN search speed (dim = 3, k = 5)"))
g = render(p)
img = PNG(string("knn_plot_", knn, ".png"), 10inch, 8inch)
draw(img, g)
=#
