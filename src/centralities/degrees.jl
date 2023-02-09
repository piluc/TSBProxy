function degrees(tg::temporal_graph, verbose::Bool)::Tuple{Array{Float64},Array{Float64},Float64}
    if (verbose)
        log("Computing degrees...")
    end
    start_time = time()
    if (verbose)
        log("    Creating underlying graph...")
    end
    ug::SimpleDiGraph = underlying_graph(tg)
    if (verbose)
        log("    Computing indegree centrality underlying graph...")
    end
    idc::Array{Float64} = indegree_centrality(ug)
    if (verbose)
        log("    Computing outdegree centrality underlying graph...")
    end
    odc::Array{Float64} = outdegree_centrality(ug)
    finish_total::Float64 = round(time() - start_time; digits=4)
    return idc, odc, finish_total
end

function temporal_degrees(tg::temporal_graph, verbose::Bool)::Tuple{Array{Float64},Array{Float64},Float64}
    if (verbose)
        log("Computing temporal degrees...")
    end
    start_time = time()
    if (verbose)
        log("    Computing temporal degrees of temporal graph...")
    end
    tidc::Array{Float64} = zeros(tg.num_nodes)
    todc::Array{Float64} = zeros(tg.num_nodes)
    for te in tg.temporal_edges
        todc[te[1]] = todc[te[1]] + 1
        tidc[te[2]] = tidc[te[2]] + 1
    end
    finish_total::Float64 = round(time() - start_time; digits=4)
    return tidc, todc, finish_total
end
