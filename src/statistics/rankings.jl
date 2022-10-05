function min_h_k(a::Array{Float64}, b::Array{Float64}, k::Int64)::Int64
    @assert length(a) == length(b) "The two rankings have different number of elements"
    ai::Vector{Int64} = sortperm(a, rev=true)
    bi::Vector{Int64} = sortperm(b, rev=true)
    min_h_b_k_a::Int64 = 0
    for a_j in 1:k
        b_j::Int64 = 1
        while (bi[b_j] != ai[a_j])
            b_j = b_j + 1
        end
        if (b_j > min_h_b_k_a)
            min_h_b_k_a = b_j
        end
    end
    return min_h_b_k_a
end

# fn = ["01_venice", "02_college_msg", "03_email_eu", "04_bordeaux", "06_infectious", "09_wiki_elections", "10_facebook_wall", "11_digg_reply"]
# is_ego = [true, true, false, true, true, true, true, true]
# fp = "graphs/www/"
# sp = "scores/"
# op = "results/"

function analyse(fp::String, sp::String, op::String, fn::Array{String}, is_ego::Array{Bool})
    mkpath(op)
    f::IOStream = open(op * "results.txt", "w")
    write(f, "Network cor_ego cor_prefix cor_ptd h_10_ego h_10_prefix h_10_ptd h_25_ego h_25_prefix h_25_ptd h_50_ego h_50_prefix h_50_ptd h_10%_ego h_10%_prefix h_10%_ptd\n")
    for g in 1:lastindex(fn)
        tg = load_temporal_graph(fp * fn[g] * ".txt", " ")
        tsb = read_centrality_values(sp * fn[g] * "/tsb.txt")
        prefix = read_centrality_values(sp * fn[g] * "/prefix.txt")
        ptd = read_centrality_values(sp * fn[g] * "/ptd.txt")
        ten_percent = Int64(round(0.01 * tg.num_nodes))
        write(f, fn[g] * " ")
        if (is_ego[g])
            ego = read_centrality_values(sp * fn[g] * "/ego.txt")
            write(f, string(round(corspearman(tsb, ego); digits=4)) * " ")
        else
            write(f, "0.0 ")
        end
        write(f, string(round(corspearman(tsb, prefix); digits=4)) * " ")
        write(f, string(round(corspearman(tsb, ptd); digits=4)) * " ")
        for k in [10, 25, 50, ten_percent]
            if (is_ego[g])
                ego = read_centrality_values(sp * fn[g] * "/ego.txt")
                write(f, string(min_h_k(tsb, ego, k)) * " ")
            else
                write(f, "0 ")
            end
            write(f, string(min_h_k(tsb, prefix, k)) * " ")
            write(f, string(min_h_k(tsb, ptd, k)) * " ")
        end
        write(f, "\n")
    end
    close(f)
end
