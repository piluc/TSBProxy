function betweenness(tg::temporal_graph, verbose::Bool)::Tuple{Array{Float64},Float64}
    if (verbose)
        log("Computing betweenness...")
    end
    start_time = time()
    if (verbose)
        log("    Creating underlying graph...")
    end
    ug::SimpleDiGraph = underlying_graph(tg)
    if (verbose)
        log("    Computing betweenness underlying graph...")
    end
    bc::Array{Float64} = betweenness_centrality(ug)
    finish_total::Float64 = round(time() - start_time; digits=4)
    return bc, finish_total
end
