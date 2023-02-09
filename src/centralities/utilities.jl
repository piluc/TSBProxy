function read_centrality_values(file_name::String)::Array{Float64}
    @assert isfile(file_name) "The centrality value file does not exist"
    f::IOStream = open(file_name, "r")
    centrality::Array{Float64} = []
    value::Float64 = 0.0
    for l in eachline(f)
        value = parse(Float64, l)
        if (value < -0.1)
            log("ERROR. There are negative values with absolute big values")
            return Array{Float64}([])
        end
        if (value < 0)
            value = 0
        end
        push!(centrality, value)
    end
    close(f)
    return centrality
end

function read_onbra_centrality_values(file_name::String, ss::Int64, nn::Int64)::Array{Float64}
    @assert isfile(file_name) "The centrality value file does not exist"
    f::IOStream = open(file_name, "r")
    centrality::Array{Float64} = zeros(nn)
    value::Float64 = 0.0
    for _ in 1:ss
        for n in 1:nn
            value = parse(Float64, readline(f))
            if (value < -0.1)
                log("ERROR. There are negative values with absolute big values")
                return Array{Float64}([])
            end
            if (value < 0)
                value = 0
            end
            centrality[n] += value / ss
        end
    end
    close(f)
    return centrality
end

function read_time(nn::String, cn::String, default_time::Float64)::Float64
    if (!isfile("times/" * nn * "/time_" * cn * ".txt"))
        return default_time
    else
        f = open("times/" * nn * "/time_" * cn * ".txt", "r")
        t::Float64 = parse(Float64, readline(f))
        close(f)
        return t
    end
end

function save_centrality_values(file_name::String, centrality::Array{Float64})::Nothing
    f::IOStream = open(file_name, "w")
    for u in 1:lastindex(centrality)
        write(f, string(centrality[u]) * "\n")
    end
    close(f)
end

function save_results(nn::String, cn::String, c::Array{Float64}, t::Float64)
    if (length(c) > 0)
        mkpath("scores/" * nn * "/")
        save_centrality_values("scores/" * nn * "/" * cn * ".txt", c)
        mkpath("times/" * nn * "/")
        f = open("times/" * nn * "/time_" * cn * ".txt", "w")
        write(f, string(t))
        close(f)
    else
        mkpath("times/" * nn * "/")
        f = open("times/" * nn * "/time_" * cn * ".txt", "w")
        write(f, "-1.0\n")
        close(f)
    end
end

function save_onbra(nn::String, c::Array{Float64}, ss::Int64, t::Tuple{Float64,Float64,Float64})
    if (length(c) > 0)
        mkpath("scores/" * nn * "/")
        save_centrality_values("scores/" * nn * "/onbra.txt", c)
        mkpath("times/" * nn * "/")
        f = open("times/" * nn * "/time_onbra.txt", "w")
        write(f, string(round(t[3]; digits=4)) * " " * string(round(t[1]; digits=4)) * " " * string(round(t[2]; digits=4)) * " " * string(ss) * "\n")
        close(f)
    else
        mkpath("times/" * nn * "/")
        f = open("times/" * nn * "/time_onbra.txt", "w")
        write(f, "-1.0 -1.0 -1.0 -1.0")
        close(f)
    end
end
