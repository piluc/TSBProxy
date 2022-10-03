function edge_time_range(tg::temporal_graph)::Tuple{Vector{Pair{Tuple{Int64,Int64},Int64}},Vector{Pair{Tuple{Int64,Int64},Int64}}}
    temporal_edge_min::Dict{Tuple{Int64,Int64},Int64} = Dict()
    temporal_edge_max::Dict{Tuple{Int64,Int64},Int64} = Dict()
    current_min::Int64 = 0
    current_max::Int64 = 0
    te::Tuple{Int64,Int64,Int64} = (0, 0, 0)
    for i in 1:lastindex(tg.temporal_edges)
        te = tg.temporal_edges[i]
        current_min = get(temporal_edge_min, (te[1], te[2]), 0)
        if (current_min == 0 || te[3] < current_min)
            temporal_edge_min[(te[1], te[2])] = te[3]
        end
        current_max = get(temporal_edge_max, (te[1], te[2]), 0)
        if (current_max == 0 || te[3] > current_max)
            temporal_edge_max[(te[1], te[2])] = te[3]
        end
    end
    return sort!(collect(temporal_edge_min), by=x -> x[2]), sort!(collect(temporal_edge_max), by=x -> x[2])
end

function count_less_than(a::Array{Int64}, t::Int64)::Int64
    left::Int64 = 1
    right::Int64 = length(a)
    pos::Int64 = 0
    mid::Int64 = -1
    while (left <= right)
        mid = (left + right) ÷ 2
        if (a[mid] < t)
            left = mid + 1
            pos = mid
        else
            right = mid - 1
        end
    end
    return pos
end

function pass_through_degree(tg::temporal_graph)::Tuple{Array{Float64},Float64}
    start::Float64 = time()
    te_min::Vector{Pair{Tuple{Int64,Int64},Int64}}, te_max::Vector{Pair{Tuple{Int64,Int64},Int64}} = edge_time_range(tg)
    min_arrival_times::Array{Array{Int64}} = [[] for j = 1:tg.num_nodes]
    head::Int64 = 0
    t::Int64 = 0
    for i in 1:lastindex(te_min)
        head = te_min[i][1][2]
        t = te_min[i][2]
        push!(min_arrival_times[head], t)
    end
    ptd::Array{Int64} = zeros(Int64, tg.num_nodes)
    for i in 1:lastindex(te_max)
        tail = te_max[i][1][1]
        t = te_max[i][2]
        ptd[tail] += count_less_than(min_arrival_times[tail], t)
    end
    finish_total::Float64 = round(time() - start; digits=4)
    return ptd, finish_total
end
