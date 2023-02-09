using DataStructures
using Dates
using Graphs
using PlotlyJS
using StatsBase

include("graphs/temporal_graph.jl")
include("centralities/betweenness.jl")
include("centralities/degrees.jl")
include("centralities/temporal_shortest_betweenness.jl")
include("centralities/onbra.jl")
include("centralities/pass_through_degree.jl")
include("centralities/prefix_foremost_betweenness.jl")
include("centralities/ego_betweenness.jl")
include("centralities/utilities.jl")
include("statistics/rankings.jl")
include("statistics/python_correlation.jl")

function log(s::String)
    # println(now(), " ", s)
    println(s)
    flush(stdout)
end
