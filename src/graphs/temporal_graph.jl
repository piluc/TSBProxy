struct temporal_graph
    num_nodes::Int64
    temporal_edges::Array{Tuple{Int64,Int64,Int64}}
    file_id::Array{String}
    file_time::Array{Int64}
end

function load_temporal_graph(file_name::String, sep::String)
    @assert isfile(file_name) "The temporal edge list file does not exist"
    current_node_id::Int64 = 1
    file_id::Vector{String} = []
    file_id_to_graph_id::Dict{String,Int64} = Dict{String,Int64}()
    current_time::Int64 = 1
    file_time::Vector{Int64} = []
    file_time_to_graph_time::Dict{Int64,Int64} = Dict{Int64,Int64}()
    temporal_edges::Array{Tuple{Int64,Int64,Int64}} = []
    t::Int64 = 0
    f::IOStream = open(file_name, "r")
    for line in eachline(f)
        split_line::Vector{String} = split(line, sep)
        @assert length(split_line) == 3 "Bad line format: " * line
        t = parse(Int64, split_line[3])
        if (!haskey(file_id_to_graph_id, split_line[1]))
            file_id_to_graph_id[split_line[1]] = current_node_id
            push!(file_id, split_line[1])
            current_node_id = current_node_id + 1
        end
        if (!haskey(file_id_to_graph_id, split_line[2]))
            file_id_to_graph_id[split_line[2]] = current_node_id
            push!(file_id, split_line[2])
            current_node_id = current_node_id + 1
        end
        if (!haskey(file_time_to_graph_time, t))
            file_time_to_graph_time[t] = current_time
            push!(file_time, t)
            current_time = current_time + 1
        end
    end
    close(f)
    sort!(file_time)
    for t in 1:lastindex(file_time)
        file_time_to_graph_time[file_time[t]] = t
    end
    f = open(file_name, "r")
    for line in eachline(f)
        split_line::Vector{String} = split(line, sep)
        t = parse(Int64, split_line[3])
        push!(temporal_edges, (file_id_to_graph_id[split_line[1]], file_id_to_graph_id[split_line[2]], file_time_to_graph_time[t]))
    end
    sort!(temporal_edges, by=te -> te[3])
    return temporal_graph(length(file_id_to_graph_id), temporal_edges, file_id, file_time)
end

function print_stats(tg::temporal_graph; graph_name="anonymous")
    println("====================================================")
    println("Temporal network: " * graph_name)
    println("====================================================")
    println("Temporal adj list loaded")
    println("Number of nodes " * string(tg.num_nodes))
    println("Number temporal of edges " * string(length(tg.temporal_edges)))
    println("Number of unique time stamps " * string(length(tg.file_time)))
    println("====================================================")
end

function temporal_adjacency_list(tg::temporal_graph)::Array{Array{Tuple{Int64,Int64}}}
    tal::Array{Array{Tuple{Int64,Int64}}} = Array{Array{Tuple{Int64,Int64}}}(undef, tg.num_nodes)
    for u in 1:tg.num_nodes
        tal[u] = Tuple{Int64,Int64,Int64}[]
    end
    te::Tuple{Int64,Int64,Int64} = (0, 0, 0)
    for i in 1:lastindex(tg.temporal_edges)
        te = tg.temporal_edges[i]
        push!(tal[te[1]], (te[2], te[3]))
    end
    return tal
end

function temporal_incidency_list(tg::temporal_graph)::Array{Array{Tuple{Int64,Int64}}}
    tal::Array{Array{Tuple{Int64,Int64}}} = Array{Array{Tuple{Int64,Int64}}}(undef, tg.num_nodes)
    for u in 1:tg.num_nodes
        tal[u] = Tuple{Int64,Int64,Int64}[]
    end
    te::Tuple{Int64,Int64,Int64} = (0, 0, 0)
    for i in 1:lastindex(tg.temporal_edges)
        te = tg.temporal_edges[i]
        push!(tal[te[2]], (te[1], te[3]))
    end
    return tal
end

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
