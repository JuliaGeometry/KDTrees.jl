using KDTrees
using StatsBase
using Plotly

function run_bench_build_tree(dim, knn, exps, rounds)
    println("Running build tree benchmark for 10^(", exps, ") points in ", dim, " dimensions.")
    n_points = [10^i for i in exps]

    times = zeros(length(n_points), rounds)

    timer = 0.0
    for (i, n_point) in enumerate(n_points)
        for (j , round) in enumerate(1:rounds)
            n_iters = 5
            println("Round ", j, " out of ", rounds, " for ", dim, "x", n_point, "...\r")
            while true
                timer = time_ns()
                for k in 1:n_iters
                  data = rand(dim, int(n_point))
                  tree = KDTree(data, 10, true)
                end
                timer = (time_ns() - float(timer)) / 10^9 # To seconds
                if timer < 1.0
                    n_iters *= 3
                    continue
                end
                break # Ends this round
            end
        times[i, j] = timer / n_iters
        end
        println("\n")
    end
    println("\nDone!")
    return mean_and_std(times, 2)
end

dim = 3
knn = 5
exps = 3:0.5:6
rounds = 3

times, stds = run_bench_build_tree(dim, knn, exps, rounds)
times = vec(times)
stds = vec(stds)


stderr = stds / sqrt(length(stds))
ymins = times - 1.96*stderr
ymaxs = times + 1.96*stderr
sizes = [10.0^i::Float64 for i in exps]

####################################################################
# Plotting
####################################################################
data = [
  [
    "x" => sizes,
    "y" => times,
    "type" => "scatter",
    "mode" => "lines+markers",
    "error_y" => [
        "type" => "data",
        "array" => stderr*1.96*2,
        "visible" => true
      ]
    ]
]

layout = [
  "title" => "Build time for KDTree (dim = 3)",
  "xaxis" => [
      "title" => "Number of points",
      "type" => "log"
  ],
  "yaxis" => [
      "title" => "Build time [s]",
      "type" => "log"
  ],
  "type" => "log",
  "autorange" => true
]


Plotly.signin("kcarlsson89", "lolololoololololo")


response = Plotly.plot(data, ["layout" => layout,
                              "filename" => "bench_build_x",
                              "fileopt" => "overwrite"])
plot_url = response["url"]



#=
using Gadfly
p = plot(x = sizes, y = times,  Geom.point, Geom.line, Geom.errorbar, Scale.x_log10, Scale.y_log10,
         ymin = ymins,
         ymax = ymaxs,
         Guide.xticks(ticks=collect(exps)),
         Guide.xlabel("Number of points"),
         Guide.ylabel("Build time [s]"),
         Guide.title("Build time for KDTree (dim = 3)"))
g = render(p)
img = PNG("build_plot.png", 5inch, 4inch)
draw(img, g)
=#


