include("../src/graphs/temporal_graph.jl")
include("../src/centralities/temporal_shortest_betweenness.jl")
include("../src/centralities/onbra.jl")
include("../src/centralities/pass_through_degree.jl")
include("../src/centralities/prefix_foremost_betweenness.jl")
include("../src/centralities/ego_shortest_betweenness.jl")
include("../src/centralities/utilities.jl")
include("../src/statistics/rankings.jl")

function one_centrality(args::Array{String})
    tg::temporal_graph = load_temporal_graph("graphs/" * args[1] * ".txt", " ")
    t::Float64 = -1.0
    if (args[2] == "ego")
        ego::Array{Float64}, t = temporal_ego_betweenness_centrality(tg, 100)
        mkpath("scores/" * args[1] * "/")
        save_centrality_values("scores/" * args[1] * "/ego.txt", ego)
        mkpath("times/" * args[1] * "/")
        f = open("times/" * args[1] * "/time_ego.txt", "w")
        write(f, string(t))
        close(f)
    elseif (args[2] == "tsb")
        tsb::Array{Float64}, t = temporal_shortest_betweenness(tg, 100)
        mkpath("scores/" * args[1] * "/")
        save_centrality_values("scores/" * args[1] * "/tsb.txt", tsb)
        mkpath("times/" * args[1] * "/")
        f = open("times/" * args[1] * "/time_tsb.txt", "w")
        write(f, string(t))
        close(f)
    elseif (args[2] == "prefix")
        prefix::Array{Float64}, t = temporal_prefix_foremost_betweenness(tg, 100)
        mkpath("scores/" * args[1] * "/")
        save_centrality_values("scores/" * args[1] * "/prefix.txt", prefix)
        mkpath("times/" * args[1] * "/")
        f = open("times/" * args[1] * "/time_prefix.txt", "w")
        write(f, string(t))
        close(f)
    elseif (args[2] == "ptd")
        ptd::Array{Float64}, t = pass_through_degree(tg)
        mkpath("scores/" * args[1] * "/")
        save_centrality_values("scores/" * args[1] * "/ptd.txt", ptd)
        mkpath("times/" * args[1] * "/")
        f = open("times/" * args[1] * "/time_ptd.txt", "w")
        write(f, string(t))
        close(f)
    end
end

# TO BE CONTINUED...
function main()
    network::Array{String} = ["01_venice"]
    centrality::Array{String} = ["tsb", "ego", "prefix", "ptd"]
    for fn in network
        println("Processing ", fn)
        one_centrality([fn, "tsb"])
    end
end

main()
