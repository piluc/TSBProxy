function jaccard(a::Array{Float64}, b::Array{Float64}, max_k::Int64)::Array{Int64}
    @assert length(a) == length(b) "The two rankings have different number of elements"
    if (max_k > length(a))
        max_k = length(a)
    end
    ai::Vector{Int64} = sortperm(a, rev=true)
    bi::Vector{Int64} = sortperm(b, rev=true)
    min_h::Array{Int64} = zeros(Int64, max_k)
    min_h_b_k_a::Int64 = 0
    for a_j in 1:max_k
        b_j::Int64 = 1
        while (bi[b_j] != ai[a_j])
            b_j = b_j + 1
        end
        if (b_j > min_h_b_k_a)
            min_h_b_k_a = b_j
        end
        min_h[a_j] = min_h_b_k_a
    end
    return min_h
end

function jaccard(a::Array{Float64}, b::Array{Float64}, max_k::Int64)::Array{Float64}
    @assert length(a) == length(b) "The two rankings have different number of elements"
    if (max_k > length(a))
        max_k = length(a)
    end
    ai::Vector{Int64} = sortperm(a, rev=true)
    bi::Vector{Int64} = sortperm(b, rev=true)
    jac::Array{Float64} = zeros(max_k)
    for k in 1:max_k
        jac[k] = length(intersect(Set(ai[1:k]), Set(bi[1:k]))) / length(union(Set(ai[1:k]), Set(bi[1:k])))
    end
    return jac
end

function jaccard(nn::String, cn1::String, cn2::String, max_k::Int64)::Array{Int64}
    @assert cn1 != "onbra" && cn2 != "onbra" "The centrality measure cannot be ONBRA"
    @assert isfile("scores/" * nn * "/" * cn1 * ".txt") "The values of " * cn1 * " are not available for the network " * nn
    @assert isfile("scores/" * nn * "/" * cn2 * ".txt") "The values of " * cn2 * " are not available for the network " * nn
    a::Array{Float64} = read_centrality_values("scores/" * nn * "/" * cn1 * ".txt")
    b::Array{Float64} = read_centrality_values("scores/" * nn * "/" * cn2 * ".txt")
    return jaccard(a, b, max_k)
end

function jaccard(nn::String, cn::String, ne::Int64, ss::Int64, max_k::Int64)::Array{Float64}
    @assert cn != "onbra" "The other centrality measure cannot be ONBRA"
    @assert isfile("scores/" * nn * "/" * cn * ".txt") "The values of " * cn1 * " are not available"
    @assert isfile("scores/" * nn * "/onbra/onbra_twice_1.txt") "The values of " * cn2 * " are not available"
    @assert isfile("scores/" * nn * "/onbra/onbra_twice_" * string(ne) * ".txt") "The values of " * cn2 * " are not available"
    all_min_h::Array{Array{Float64}} = []
    a::Array{Float64} = read_centrality_values("scores/" * nn * "/" * cn * ".txt")
    for e in 1:ne
        b::Array{Float64} = read_onbra_centrality_values("scores/" * nn * "/onbra/onbra_twice_" * string(e) * ".txt", ss, length(a))
        push!(all_min_h, jaccard(a, b, max_k))
    end
    min_h::Array{Float64} = zeros(max_k)
    for k in 1:max_k
        for e in 1:ne
            min_h[k] = min_h[k] + all_min_h[e][k]
        end
        min_h[k] = min_h[k] / ne
    end
    return min_h
end

function jaccard(nn::String, cn1::String, cn2::String, max_k::Int64)::Array{Float64}
    @assert cn1 != "onbra" && cn2 != "onbra" "The centrality measure cannot be ONBRA"
    @assert isfile("scores/" * nn * "/" * cn1 * ".txt") "The values of " * cn1 * " are not available for the network " * nn
    @assert isfile("scores/" * nn * "/" * cn2 * ".txt") "The values of " * cn2 * " are not available for the network " * nn
    a::Array{Float64} = read_centrality_values("scores/" * nn * "/" * cn1 * ".txt")
    b::Array{Float64} = read_centrality_values("scores/" * nn * "/" * cn2 * ".txt")
    return jaccard(a, b, max_k)
end

function jaccard(nn::String, cn::String, ne::Int64, ss::Int64, max_k::Int64)::Array{Float64}
    @assert cn != "onbra" "The other centrality measure cannot be ONBRA"
    @assert isfile("scores/" * nn * "/" * cn * ".txt") "The values of " * cn * " are not available"
    @assert isfile("scores/" * nn * "/onbra/onbra_twice_1.txt") "The values of " * cn2 * " are not available"
    @assert isfile("scores/" * nn * "/onbra/onbra_twice_" * string(ne) * ".txt") "The values of " * cn2 * " are not available"
    all_jaccard::Array{Array{Float64}} = []
    a::Array{Float64} = read_centrality_values("scores/" * nn * "/" * cn * ".txt")
    for e in 1:ne
        b::Array{Float64} = read_onbra_centrality_values("scores/" * nn * "/onbra/onbra_twice_" * string(e) * ".txt", ss, length(a))
        push!(all_jaccard, jaccard(a, b, max_k))
    end
    jac::Array{Float64} = zeros(max_k)
    for k in 1:max_k
        for e in 1:ne
            jac[k] = jac[k] + all_jaccard[e][k]
        end
        jac[k] = jac[k] / ne
    end
    return jac
end
