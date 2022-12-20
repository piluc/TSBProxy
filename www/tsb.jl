using Graphs

include("../src/graphs/temporal_graph.jl")
include("../src/centralities/temporal_shortest_betweenness.jl")
include("../src/centralities/onbra.jl")
include("../src/centralities/pass_through_degree.jl")
include("../src/centralities/prefix_foremost_betweenness.jl")
include("../src/centralities/ego_betweenness.jl")
include("../src/centralities/utilities.jl")
include("../src/statistics/rankings.jl")
include("python_correlation.jl")

"""
 Save the values of a (not ONBRA) centrality measure on one file (one line per line) and the corresponding execution time on another file (if there are no values write -1.0)
"""
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

"""
Save the ONBRA values corresponding to an experiment on one file (for each sample pair, one line per node, so that the file contains a number of lines which is the product of the sample size times the number of nodes) and append the corresponding average execution time (along with the standard deviation, the total execution time, and the sample size) on another file (if there are no values write all -1.0)
"""
function save_onbra_results(nn::String, c::Array{Float64}, type::String, ss::Int64, e::Int64, t::Tuple{Float64,Float64,Float64})
    if (length(c) > 0)
        mkpath("scores/" * nn * "/onbra/")
        save_centrality_values("scores/" * nn * "/onbra/onbra_" * type * "_" * string(e) * ".txt", c)
        mkpath("times/" * nn * "/onbra/")
        f = open("times/" * nn * "/onbra/time_" * type * ".txt", "a")
        write(f, string(round(t[3]; digits=4)) * " " * string(round(t[1]; digits=4)) * " " * string(round(t[2]; digits=4)) * " " * string(ss) * "\n")
        close(f)
    else
        mkpath("times/" * nn * "/onbra/")
        f = open("times/" * nn * "/onbra/time_" * type * ".txt", "a")
        write(f, "-1.0 -1.0 -1.0 -1.0")
        close(f)
    end
end

"""
Read the execution time for computing one kind of centrality values
"""
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

"""
Compute and return the average execution time (with standard deviation) of ONBRA with a specific type of sample size
"""
function read_onbra_time(nn::String, te::String, nt::Int64, default_time::Float64)::Tuple{Float64,Float64}
    if (!isfile("times/" * nn * "/onbra/time_" * te * ".txt"))
        return default_time, 0.0
    else
        t::Array{Float64} = zeros(nt)
        f = open("times/" * nn * "/onbra/time_" * te * ".txt", "r")
        i::Int64 = 1
        for l in eachline(f)
            t[i] = parse(Float64, split(l, " ")[1])
            if (t[1] < 0)
                return default_time, 0.0
            end
            i += 1
        end
        close(f)
        return mean(t), std(t)
    end
end

"""
Compute and return the average Spearman correlation (with standard deviation) of ONBRA (with a specific type of sample size) and the TSB values (another measure could be used as well) 
"""
function compute_onbra_correlations(nn::String, te::String, nt::Int64, tsb::Array{Float64})
    if (!isfile("times/" * nn * "/onbra/time_" * te * ".txt"))
        return 0.0, 0.0
    end
    f = open("times/" * nn * "/onbra/time_" * te * ".txt", "r")
    sl = split(readline(f), " ")
    close(f)
    if (parse(Float64, sl[1]) < 0.0)
        return 0.0, 0.0
    end
    ss::Int64 = parse(Int64, sl[4])
    onbra_spearman = zeros(nt)
    onbra_ktau = zeros(nt)
    onbra_wktau = zeros(nt)
    for e in 1:nt
        if (!isfile("scores/" * nn * "/onbra/onbra_" * te * "_" * string(e) * ".txt"))
            return (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)
        end
        onbra = read_onbra_centrality_values("scores/" * nn * "/onbra/onbra_" * te * "_" * string(e) * ".txt", ss, length(tsb))
        correlations = compute_correlations(tsb, onbra, true)
        onbra_spearman[e] = correlations[1][1]
        onbra_ktau[e] = correlations[2][1]
        onbra_wktau[e] = correlations[3][1]
    end
    return (mean(onbra_spearman), mean(onbra_ktau), mean(onbra_wktau)), (std(onbra_spearman), std(onbra_ktau), std(onbra_wktau))
end

"""
Compute and return the average minimum value h (with standard deviation) such that the k-th element in the TSB ranking (another measure could be used as well) is in position h in the ONBRA (with a specific type of sample size) ranking 
"""
function compute_onbra_min_h_k(nn::String, te::String, nt::Int64, tsb::Array{Float64}, h::Int64)::Tuple{Float64,Float64}
    if (!isfile("times/" * nn * "/onbra/time_" * te * ".txt"))
        return 0.0, 0.0
    end
    f = open("times/" * nn * "/onbra/time_" * te * ".txt", "r")
    sl = split(readline(f), " ")
    close(f)
    if (parse(Float64, sl[1]) < 0.0)
        return 0.0, 0.0
    end
    ss::Int64 = parse(Int64, sl[4])
    onbra_min_h_k = zeros(nt)
    for e in 1:nt
        if (!isfile("scores/" * nn * "/onbra/onbra_" * te * "_" * string(e) * ".txt"))
            return 0.0, 0.0
        end
        onbra = read_onbra_centrality_values("scores/" * nn * "/onbra/onbra_" * te * "_" * string(e) * ".txt", ss, length(tsb))
        onbra_min_h_k[e] = min_h_k(tsb, onbra, h)
    end
    return mean(onbra_min_h_k), std(onbra_min_h_k)
end

"""
Compute and save the values of the specified centrality with respect to the specified network file 
"""
function one_centrality(nn::String, cn::String, bigint::Bool)
    tg::temporal_graph = load_temporal_graph("graphs/" * nn * ".txt", " ")
    aens::Float64 = 0.0
    if (cn == "egotsb" || cn == "egoprefix")
        aens = average_ego_network_size(tg)
    end
    t::Float64 = -1.0
    max_time::Float64 = read_time(nn, "tsb", 1000.0)
    centrality::Array{Float64} = zeros(tg.num_nodes)
    if (cn == "egotsb")
        if (aens > tg.num_nodes / 10)
            centrality, t = Array{Float64}([]), 0.0
        else
            centrality, t = temporal_ego_betweenness_centrality(tg, max_time, 100)
        end
    elseif (cn == "egoprefix")
        if (aens > tg.num_nodes / 100)
            centrality, t = Array{Float64}([]), 0.0
        else
            centrality, t = temporal_ego_prefix_foremost_betweenness(tg, max_time, 100)
        end
    elseif (cn == "prefix")
        centrality, t = temporal_prefix_foremost_betweenness(tg, 100)
    elseif (cn == "ptd")
        centrality, t = pass_through_degree(tg)
    elseif (cn == "tsb")
        centrality, t = temporal_shortest_betweenness(tg, 100, bigint)
    end
    save_results(nn, cn, centrality, t)
end

"""
Compute and save the values of all centralities apart from ONBRA 
"""
function execute_all_but_onbra(network_name::Array{String}, bigint::Bool)
    centrality_name::Array{String} = ["tsb", "egotsb", "egoprefix", "prefix", "ptd"]
    for fn in network_name
        println("Processing ", fn)
        for cn in centrality_name
            println("Processing ", cn)
            one_centrality(fn, cn, bigint)
        end
    end
end

"""
Compute and save the values of all versions of ONBRA (based on the execution time of PREFIX) 
"""
function execute_all_onbras(network_name::Array{String}, bigint::Bool)
    type::Array{String} = ["equal", "twice", "half"]
    factor::Array{Float64} = [1, 2, 0.5]
    time_estimate_trials::Int64 = 100
    num_trials::Int64 = 10
    for fn in network_name
        println("Processing ", fn)
        tg::temporal_graph = load_temporal_graph("graphs/" * fn * ".txt", " ")
        _, t = onbra(tg, time_estimate_trials, 0, bigint)
        println("ONBRA average one sample execution time: ", t[1], " (", t[2], ")")
        prefix_time::Float64 = read_time(fn, "prefix", 1.0)
        for te in 1:lastindex(type)
            sample_size::Int64 = Int64(round(factor[te] * prefix_time / t[1]))
            println("ONBRA sample size (" * type[te] * "): ", sample_size)
            execute_one_onbra(fn, tg, sample_size, type[te], num_trials, 0, bigint)
        end
    end
end

"""
Compute and save the values of one version of ONBRA (based on the execution time of PREFIX) 
"""
function execute_one_onbra(nn::String, tg::temporal_graph, ss::Int64, type::String, nt::Int64, verb::Int64, bigint::Bool)
    onbra_value::Array{Float64} = zeros(tg.num_nodes)
    onbra_time::Tuple{Float64,Float64,Float64} = (0.0, 0.0, 0.0)
    if (isfile("times/" * nn * "/onbra/time_" * type * ".txt"))
        rm("times/" * nn * "/onbra/time_" * type * ".txt")
    end
    for e in 1:nt
        print("ONBRA experiment " * string(e) * "/" * string(nt) * ". ")
        onbra_value, onbra_time = onbra(tg, ss, verb, bigint)
        save_onbra_results(nn, onbra_value, type, ss, e, onbra_time)
        println("ONBRA with " * string(ss) * " seeds computed in " * string(round(onbra_time[3]; digits=4)) * " seconds")
    end
end

"""
Analyse results on specified networks for all centrality measures apart from ONBRA
"""
function analyse_all_but_onbra_with_min_h_k(network_name::Array{String})
    centrality_name::Array{String} = ["egotsb", "egoprefix", "prefix", "ptd"]
    for fn in network_name
        mkpath("evaluation/")
        println("Analysing ", fn)
        t::Float64 = read_time(fn, "tsb", -1.0)
        if (t >= 0.0)
            f::IOStream = open("evaluation/" * fn * "_not_onbra.txt", "w")
            write(f, "Network:time_exact:time_egotsb:time_egoprefix:time_prefix:time_ptd:spearman_egotsb:spearman_egoprefix:spearman_prefix:spearman_ptd:h_10_egotsb:h_50_egotsb:h_100_egotsb:h_tenpercent_egotsb:h_10_egoprefix:h_50_egoprefix:h_100_egoprefix:h_tenpercent_egoprefix:h_10_prefix:h_50_prefix:h_100_prefix:h_tenpercent_prefix:h_10_ptd:h_50_ptd:h_100_ptd:h_tenpercent_ptd\n")
            write(f, fn[4:end] * ":" * string(round(t; digits=4)))
            for cn in centrality_name
                t = read_time(fn, cn, -1.0)
                if (t >= 0.0)
                    write(f, ":" * string(round(t; digits=4)))
                else
                    write(f, ":0.0")
                end
            end
            tsb::Array{Float64} = read_centrality_values("scores/" * fn * "/tsb.txt")
            for cn in centrality_name
                if (isfile("scores/" * fn * "/" * cn * ".txt"))
                    c::Array{Float64} = read_centrality_values("scores/" * fn * "/" * cn * ".txt")
                    write(f, ":" * string(round(corspearman(tsb, c); digits=2)))
                else
                    write(f, ":0.0")
                end
            end
            for cn in centrality_name
                if (isfile("scores/" * fn * "/" * cn * ".txt"))
                    c::Array{Float64} = read_centrality_values("scores/" * fn * "/" * cn * ".txt")
                    write(f, ":" * string(min_h_k(tsb, c, 10)))
                    write(f, ":" * string(min_h_k(tsb, c, 50)))
                    write(f, ":" * string(min_h_k(tsb, c, 100)))
                    write(f, ":" * string(min_h_k(tsb, c, Int64(round(0.01 * length(c))))))
                else
                    write(f, ":0:0:0:0")
                end
            end
            close(f)
        else
            println("ERROR. No TSB execution time.")
        end
    end
end

"""
Analyse results on specified networks for all centrality measures apart from ONBRA
"""
function analyse_all_but_onbra(network_name::Array{String})
    centrality_name::Array{String} = ["egotsb", "egoprefix", "prefix", "ptd"]
    for fn in network_name
        mkpath("evaluation/tmp/")
        println("Analysing ", fn)
        t::Float64 = read_time(fn, "tsb", -1.0)
        if (t >= 0.0)
            f::IOStream = open("evaluation/tmp/" * fn * "_not_onbra.txt", "w")
            write(f, "Network:time_exact:time_egotsb:spearman_egotsb:ktau_egotsb:wktau_egotsb:time_egoprefix:spearman_egoprefix:ktau_egoprefix:wktau_egoprefix:time_prefix:spearman_prefix:ktau_prefix:wktau_prefix:time_ptd:spearman_ptd:ktau_ptd:wktau_ptd\n")
            write(f, fn[4:end] * ":" * string(round(t; digits=4)))
            tsb::Array{Float64} = read_centrality_values("scores/" * fn * "/tsb.txt")
            for cn in centrality_name
                t = read_time(fn, cn, -1.0)
                if (t >= 0.0)
                    write(f, ":" * string(round(t; digits=4)))
                else
                    write(f, ":0.0000")
                end
                if (isfile("scores/" * fn * "/" * cn * ".txt"))
                    c::Array{Float64} = read_centrality_values("scores/" * fn * "/" * cn * ".txt")
                    correlations = compute_correlations(tsb, c, true)
                    for cor in 1:3
                        write(f, ":" * string(round(correlations[cor][1]; digits=2)))
                    end
                else
                    write(f, ":0.00:0.00:0.00")
                end
            end
            close(f)
        else
            println("ERROR. No TSB execution time.")
        end
    end
end

"""
Analyse results on specified networks for ONBRA
"""
function analyse_all_onbras(network_name::Array{String}, type::Array{String})
    num_trials::Int64 = 10
    for fn in network_name
        mkpath("evaluation/tmp/")
        println("Analysing ", fn)
        t::Float64 = read_time(fn, "tsb", -1.0)
        if (t >= 0.0)
            tsb::Array{Float64} = read_centrality_values("scores/" * fn * "/tsb.txt")
            f::IOStream = open("evaluation/tmp/" * fn * "_onbra.txt", "w")
            write(f, "Network:time_exact")
            for te in 1:lastindex(type)
                write(f, ":avg_time_onbra_" * type[te] * ":std_time_onbra_" * type[te] * ":avg_spearman_onbra_" * type[te] * ":std_spearman_onbra_" * type[te] * ":avg_ktau_onbra_" * type[te] * ":std_ktau_onbra_" * type[te] * ":avg_wktau_onbra_" * type[te] * ":std_wktau_onbra_" * type[te])
            end
            write(f, "\n" * fn[4:end] * ":" * string(round(t; digits=4)))
            for te in 1:lastindex(type)
                println("ONBRA " * type[te])
                avg_t::Float64, std_t::Float64 = read_onbra_time(fn, type[te], num_trials, -1.0)
                if (avg_t >= 0.0)
                    write(f, ":" * string(round(avg_t; digits=4)) * ":" * string(round(std_t; digits=4)))
                else
                    write(f, ":0.0000:0.0000")
                end
                println("    Times read")
                avg_cors, std_cors = compute_onbra_correlations(fn, type[te], num_trials, tsb)
                for cor in 1:3
                    write(f, ":" * string(round(avg_cors[cor]; digits=2)) * ":" * string(round(std_cors[cor]; digits=2)))
                end
            end
            close(f)
        else
            println("ERROR. No TSB execution time.")
        end
    end
end

"""
For all specified conferences, merge the two files of the analysis of non-ONBRA and ONBRA results into a unique file (one for all conferences) to be used by the R script
"""
function merge_analysis(network_name::Array{String})
    f::IOStream = open("evaluation/results.csv", "w")
    for fn in network_name
        if (isfile("evaluation/" * fn * "_not_onbra.txt") && isfile("evaluation/" * fn * "_onbra.txt"))
            lno::Array{Array{String}} = Array{String}(undef, 2, 26)
            lo::Array{Array{String}} = Array{String}(undef, 2, 38)
            fno::IOStream = open("evaluation/" * fn * "_not_onbra.txt")
            lno[1] = split(readline(fno), ":")
            lno[2] = split(readline(fno), ":")
            close(fno)
            fo::IOStream = open("evaluation/" * fn * "_onbra.txt")
            lo[1] = split(readline(fo), ":")
            lo[2] = split(readline(fo), ":")
            close(fo)
            bi::Int64 = 2
            if (fn == network_name[1])
                bi = 1
            end
            for l in bi:2
                for i in 1:6
                    write(f, lno[l][i] * ":")
                end
                for i in 3:8
                    write(f, lo[l][i] * ":")
                end
                for i in 7:10
                    write(f, lno[l][i] * ":")
                end
                for i in 9:14
                    write(f, lo[l][i] * ":")
                end
                for i in 11:(lastindex(lno[l])-1)
                    write(f, lno[l][i] * ":")
                end
                write(f, lno[l][lastindex(lno[l])])
                for i in 15:(lastindex(lo[l])-1)
                    write(f, lo[l][i] * ":")
                end
                write(f, lo[l][lastindex(lo[l])])
                write(f, "\n")
            end
        else
            println("ERROR. Results missing for " * fn)
        end
    end
    close(f)
end

function onbra_evolution(network_name::Array{String}, nt::Int64, as::Int64)
    mkpath("evaluation/tmp/")
    for fn in network_name
        println("Processing ", fn)
        if (isfile("evaluation/tmp/" * fn * "_onbra.txt"))
            mkpath("evaluation/onbra_evolution/" * fn * "/")
            of::IOStream = open("evaluation/onbra_evolution/" * fn * "/evolution.txt", "w")
            write(of, "num_samples:time:spearman:tau:wtau\n")
            tsb::Array{Float64} = read_centrality_values("scores/" * fn * "/tsb.txt")
            nn::Int64 = length(tsb)
            cum_onbra::Array{Float64} = zeros(nn)
            onbra::Array{Float64} = zeros(nn)
            ss_time, _ = read_onbra_time(fn, "twice", nt, 0.0)
            f::IOStream = open("times/" * fn * "/onbra/time_twice.txt", "r")
            l::String = readline(f)
            close(f)
            ss::Int64 = parse(Int64, split(l, " ")[4])
            sample_t::Float64 = ss_time / ss
            rem::Int64 = 0
            ns::Int64 = div((ss - rem), as)
            total_cum_ss::Int64 = 0
            cum_ss::Int64 = 0
            for e in 1:nt
                cum_ss = 0
                f = open("scores/" * fn * "/onbra/onbra_twice_" * string(e) * ".txt")
                if (rem > 0)
                    for _ in 1:rem
                        for u in 1:nn
                            cum_onbra[u] += parse(Float64, readline(f))
                        end
                        total_cum_ss += 1
                        cum_ss += 1
                    end
                    for u in 1:nn
                        onbra[u] = cum_onbra[u] / total_cum_ss
                    end
                    correlations = compute_correlations(tsb, onbra, false)
                    write(of, string(total_cum_ss) * ":" * string(round(total_cum_ss * sample_t; digits=4)) * ":" * string(correlations[1][1]) * ":" * string(correlations[2][1]) * ":" * string(correlations[3][1]) * "\n")
                end
                for _ in 1:ns
                    for _ in 1:as
                        for u in 1:nn
                            cum_onbra[u] += parse(Float64, readline(f))
                        end
                        total_cum_ss += 1
                        cum_ss += 1
                    end
                    for u in 1:nn
                        onbra[u] = cum_onbra[u] / total_cum_ss
                    end
                    correlations = compute_correlations(tsb, onbra, false)
                    write(of, string(total_cum_ss) * ":" * string(round(total_cum_ss * sample_t; digits=4)) * ":" * string(correlations[1][1]) * ":" * string(correlations[2][1]) * ":" * string(correlations[3][1]) * "\n")
                end
                rem = ss - cum_ss
                for _ in 1:rem
                    for u in 1:nn
                        cum_onbra[u] += parse(Float64, readline(f))
                    end
                    total_cum_ss += 1
                    cum_ss += 1
                end
                if (cum_ss != ss)
                    println("ERROR. The number of samples in a file does not correspond to the number of read samples.")
                end
                rem = as - rem
                ns = div((ss - rem), as)
                close(f)
            end
            close(of)
        end
    end
end

function ne_underlying_graph(network_name::String)::Int64
    tg::temporal_graph = load_temporal_graph("graphs/" * network_name * ".txt", " ")
    g::SimpleDiGraph = SimpleDiGraph(tg.num_nodes)
    for te in tg.temporal_edges
        add_edge!(g, te[1], te[2])
    end
    return ne(g)
end
# network_name::Array{String} = ["00_hospital_ward", "01_venice", "02_college_msg", "03_email_eu", "04_bordeaux", "05_adelaide", "06_infectious", "07_SMS", "08_topology", "09_wiki_elections", "10_facebook_wall", "11_digg_reply", "12_mathoverflow"]
# type::Array{String} = ["half", "equal", "twice"]

# execute_all_but_onbra(["04_bordeaux"], true)
# execute_all_onbras(network_name)
# analyse_all_but_onbra(network_name)
# analyse_all_onbras(network_name, type)
# merge_analysis(network_name)

function table_1(network_name::String)
    tg::temporal_graph = load_temporal_graph("graphs/" * network_name * ".txt", " ")
    print_stats(tg)
end
