include("../src/graphs/temporal_graph.jl")
include("../src/centralities/temporal_shortest_betweenness.jl")
include("../src/centralities/onbra.jl")
include("../src/centralities/pass_through_degree.jl")
include("../src/centralities/prefix_foremost_betweenness.jl")
include("../src/centralities/ego_betweenness.jl")
include("../src/centralities/utilities.jl")
include("../src/statistics/rankings.jl")

function save_results(nn::String, cn::String, c::Array{Float64}, t::Float64)
    if (length(c) > 0)
        mkpath("scores/" * nn * "/")
        save_centrality_values("scores/" * nn * "/" * cn * ".txt", c)
        mkpath("times/" * nn * "/")
        f = open("times/" * nn * "/time_" * cn * ".txt", "w")
        write(f, string(t))
        close(f)
    else
        mkpath("times/" * nn * "/")
        f = open("times/" * nn * "/time_" * cn * ".txt", "w")
        write(f, uppercase(cn) * " not computed because too expensive or because of overflow")
        close(f)
    end
end

function tsb_time(nn::String, default_time::Float64)
    if (!isfile("times/" * nn * "/time_tsb.txt"))
        return default_time
    else
        f = open("times/" * nn * "/time_tsb.txt", "r")
        mt::Float64 = parse(Float64, readline(f))
        close(f)
        return mt
    end
end

function one_centrality(args::Array{String})
    tg::temporal_graph = load_temporal_graph("graphs/" * args[1] * ".txt", " ")
    aens::Float64 = average_ego_network_size(tg)
    t::Float64 = -1.0
    max_time::Float64 = tsb_time(args[1], 1000.0)
    centrality::Array{Float64} = zeros(tg.num_nodes)
    if (args[2] == "egotsb")
        if (aens > tg.num_nodes / 10)
            save_results(args[1], args[2], Array{Float64}([]), 0.0)
        else
            centrality, t = temporal_ego_betweenness_centrality(tg, max_time, 100)
            save_results(args[1], args[2], centrality, t)
        end
    elseif (args[2] == "egoprefix")
        if (aens > tg.num_nodes / 100)
            save_results(args[1], args[2], Array{Float64}([]), 0.0)
        else
            centrality, t = temporal_ego_prefix_foremost_betweenness(tg, max_time, 100)
            save_results(args[1], args[2], centrality, t)
        end
    elseif (args[2] == "prefix")
        centrality, t = temporal_prefix_foremost_betweenness(tg, 100)
        save_results(args[1], args[2], centrality, t)
    elseif (args[2] == "ptd")
        centrality, t = pass_through_degree(tg)
        save_results(args[1], args[2], centrality, t)
    elseif (args[2] == "tsb")
        centrality, t = temporal_shortest_betweenness(tg, 100)
        save_results(args[1], args[2], centrality, t)
    end
end

function main()
    network_name::Array{String} = ["00_hospital_ward", "01_venice", "02_college_msg", "03_email_eu", "04_bordeaux", "05_adelaide", "06_infectious", "07_SMS", "08_topology", "09_wiki_elections", "10_facebook_wall", "11_digg_reply", "12_mathoverflow"]
    centrality_name::Array{String} = ["tsb", "egotsb", "egoprefix", "prefix", "ptd"]
    for fn in network_name
        println("Processing ", fn)
        for cn in centrality_name
            println("Processing ", cn)
            one_centrality([fn, cn])
        end
    end
end

main()
