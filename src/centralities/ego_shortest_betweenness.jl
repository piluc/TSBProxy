using DataStructures

function ego_network(tg::temporal_graph, e::Int64)::temporal_graph
    temporal_ego_edges::Array{Tuple{Int64,Int64,Int64}} = Array{Tuple{Int64,Int64,Int64}}([])
    current_node_id::Int64 = 2
    graph_id::Vector{String} = Vector{String}([])
    graph_id_to_ego_id::Dict{Int64,Int64} = Dict{Int64,Int64}()
    graph_time::Vector{Int64} = Vector{Int64}([])
    graph_time_to_ego_time::Dict{Int64,Int64} = Dict{Int64,Int64}()
    i::Int64 = -1
    current_time::Int64 = 1
    graph_id_to_ego_id[e] = 1
    push!(graph_id, tg.file_id[e])
    for edge in tg.temporal_edges
        if (edge[1] == e || edge[2] == e)
            i = (edge[1] == e) ? 2 : 1
            if (!haskey(graph_id_to_ego_id, edge[i]))
                graph_id_to_ego_id[edge[i]] = current_node_id
                push!(graph_id, tg.file_id[edge[i]])
                current_node_id = current_node_id + 1
            end
            if (!haskey(graph_time_to_ego_time, edge[3]))
                graph_time_to_ego_time[edge[3]] = current_time
                push!(graph_time, edge[3])
                current_time = current_time + 1
            end
        end
    end
    for edge in tg.temporal_edges
        if (edge[1] != e && edge[2] != e && haskey(graph_id_to_ego_id, edge[1]) && haskey(graph_id_to_ego_id, edge[2]))
            if (!haskey(graph_time_to_ego_time, edge[3]))
                graph_time_to_ego_time[edge[3]] = current_time
                push!(graph_time, edge[3])
                current_time = current_time + 1
            end
        end
    end
    sort!(graph_time)
    for t in 1:lastindex(graph_time)
        graph_time_to_ego_time[graph_time[t]] = t
    end
    for edge in tg.temporal_edges
        if (haskey(graph_id_to_ego_id, edge[1]) && haskey(graph_id_to_ego_id, edge[2]))
            push!(temporal_ego_edges, (graph_id_to_ego_id[edge[1]], graph_id_to_ego_id[edge[2]], graph_time_to_ego_time[edge[3]]))
        end
    end
    return temporal_graph(length(graph_id_to_ego_id), temporal_ego_edges, graph_id, graph_time)
end

function temporal_ego_betweenness_centrality(tg::temporal_graph, verbose_step::Int64)::Tuple{Array{Float64},Float64}
    start_time = time()
    temporal_ego_betweenness_centrality::Array{Float64} = zeros(tg.num_nodes)
    processed_so_far::Int64 = 0
    for e in 1:tg.num_nodes
        temporal_ego_betweenness_centrality[e], _ = temporal_shortest_betweenness(ego_network(tg, e), 0)[1]
        processed_so_far += 1
        if (verbose_step > 0 && processed_so_far % verbose_step == 0)
            finish_partial::String = string(round(time() - start_time; digits=4))
            println("Processed " * string(processed_so_far) * "/" * string(tg.num_nodes) * " ego in " * finish_partial * " seconds")
        end
    end
    finish_total::Float64 = round(time() - start_time; digits=4)
    return temporal_ego_betweenness_centrality, finish_total
end
