using DataStructures

function temporal_ego_betweenness_centrality(tg::temporal_graph, verbose_step::Int64)::Tuple{Array{Float64},Float64}
    start_time = time()
    temporal_ego_betweenness_centrality::Array{Float64} = zeros(tg.num_nodes)
    processed_so_far::Int64 = 0
    for e in 1:tg.num_nodes
        temporal_ego_betweenness_centrality[e] = temporal_shortest_betweenness(ego_network(tg, e), 0)[1][1]
        processed_so_far += 1
        if (verbose_step > 0 && processed_so_far % verbose_step == 0)
            finish_partial::String = string(round(time() - start_time; digits=4))
            println("Processed " * string(processed_so_far) * "/" * string(tg.num_nodes) * " ego in " * finish_partial * " seconds")
        end
    end
    finish_total::Float64 = round(time() - start_time; digits=4)
    return temporal_ego_betweenness_centrality, finish_total
end
