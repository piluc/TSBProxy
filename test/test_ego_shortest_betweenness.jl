using Test

include("../src/graphs/temporal_graph.jl")
include("../src/centralities/temporal_shortest_betweenness.jl")
include("../src/centralities/ego_betweenness.jl")

@testset "ego_shortest_betweenness" begin
    @testset "from star file" begin
        tg::temporal_graph = load_temporal_graph("test/graphs/star.txt", " ")
        tsb::Array{Float64}, _ = temporal_ego_betweenness_centrality(tg, 0)
        @test isapprox(collect.(tsb), collect.([0.0, 20.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))
    end
    @testset "from path file" begin
        tg::temporal_graph = load_temporal_graph("test/graphs/path.txt", " ")
        tsb::Array{Float64}, _ = temporal_ego_betweenness_centrality(tg, 0)
        @test isapprox(collect.(tsb), collect.([0.0, 1.0, 1.0, 1.0, 1.0, 0.0]))
    end
    @testset "from lollipop file" begin
        tg::temporal_graph = load_temporal_graph("test/graphs/lollipop.txt", " ")
        tsb::Array{Float64}, _ = temporal_ego_betweenness_centrality(tg, 0)
        @test isapprox(collect.(tsb), collect.([0.0, 0.0, 0.0, 3.0, 1.0, 1.0, 1.0, 0.0]))
    end
    @testset "from multipath file" begin
        tg::temporal_graph = load_temporal_graph("test/graphs/multipath.txt", " ")
        tsb::Array{Float64}, _ = temporal_ego_betweenness_centrality(tg, 0)
        @test isapprox(collect.(tsb), collect.([0.0, 1.0, 1.0, 2.0, 0.0]))
    end
end;
