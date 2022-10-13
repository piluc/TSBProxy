using DataStructures

function temporal_ego_betweenness_centrality(tg::temporal_graph, max_time::Float64, verbose_step::Int64)::Tuple{Array{Float64},Float64}
    start_time = time()
    temporal_ego_centrality::Array{Float64} = zeros(tg.num_nodes)
    processed_so_far::Int64 = 0
    for e in 1:tg.num_nodes
        temporal_ego_centrality[e] = temporal_shortest_betweenness(ego_network(tg, e), 0, false)[1][1]
        if ((time() - start_time) > max_time)
            return Array{Float64}[], 0.0
        end
        processed_so_far += 1
        if (verbose_step > 0 && processed_so_far % verbose_step == 0)
            finish_partial::String = string(round(time() - start_time; digits=4))
            println("EGOTSB. Processed " * string(processed_so_far) * "/" * string(tg.num_nodes) * " ego in " * finish_partial * " seconds")
        end
    end
    finish_total::Float64 = round(time() - start_time; digits=4)
    return temporal_ego_centrality, finish_total
end

function temporal_ego_prefix_foremost_betweenness(tg::temporal_graph, max_time::Float64, verbose_step::Int64)::Tuple{Array{Float64},Float64}
    start_time = time()
    temporal_ego_centrality::Array{Float64} = zeros(tg.num_nodes)
    processed_so_far::Int64 = 0
    for e in 1:tg.num_nodes
        temporal_ego_centrality[e] = temporal_prefix_foremost_betweenness(ego_network(tg, e), 0)[1][1]
        if ((time() - start_time) > max_time)
            return Array{Float64}[], 0.0
        end
        processed_so_far += 1
        if (verbose_step > 0 && processed_so_far % verbose_step == 0)
            finish_partial::String = string(round(time() - start_time; digits=4))
            println("EGOPREFIX. Processed " * string(processed_so_far) * "/" * string(tg.num_nodes) * " ego in " * finish_partial * " seconds")
        end
    end
    finish_total::Float64 = round(time() - start_time; digits=4)
    return temporal_ego_centrality, finish_total
end
