using Test

include("../src/graphs/temporal_graph.jl")
include("../src/centralities/temporal_shortest_betweenness.jl")
include("../src/centralities/prefix_foremost_betweenness.jl")

@testset "prefix_foremost_betweenness" begin
    @testset "from star file" begin
        tg::temporal_graph = load_temporal_graph("test/graphs/star.txt", " ")
        tsb::Array{Float64}, _ = temporal_prefix_foremost_betweenness(tg, 0)
        @test isapprox(collect.(tsb), collect.([0.0, 20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))
    end
    @testset "from path file" begin
        tg::temporal_graph = load_temporal_graph("test/graphs/path.txt", " ")
        tsb::Array{Float64}, _ = temporal_prefix_foremost_betweenness(tg, 0)
        @test isapprox(collect.(tsb), collect.([0.0, 4.0, 6.0, 6.0, 4.0, 0.0]))
    end
    @testset "from lollipop file" begin
        tg::temporal_graph = load_temporal_graph("test/graphs/lollipop.txt", " ")
        tsb::Array{Float64}, _ = temporal_prefix_foremost_betweenness(tg, 0)
        @test isapprox(collect.(tsb), collect.([0.0, 5, 0, 12, 12, 10, 6.0, 0]))
    end
    @testset "from multipath file" begin
        tg::temporal_graph = load_temporal_graph("test/graphs/multipath.txt", " ")
        tsb::Array{Float64}, _ = temporal_prefix_foremost_betweenness(tg, 0)
        @test isapprox(collect.(tsb), collect.([0.0, 2, 0, 3.0, 0]))
    end
end;
