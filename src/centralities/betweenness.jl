using DataStructures
using Graphs

function betweenness(tg::temporal_graph, verbose::Bool)::Tuple{Array{Float64},Float64}
    if (verbose)
        println("Computing betweenness...")
    end
    start_time = time()
    if (verbose)
        println("    Creating underlying graph...")
    end
    ug::SimpleDiGraph = underlying_graph(tg)
    if (verbose)
        println("    Computing betweenness underlying graph...")
    end
    bc::Array{Float64} = betweenness_centrality(ug)
    finish_total::Float64 = round(time() - start_time; digits=4)
    return bc, finish_total
end
