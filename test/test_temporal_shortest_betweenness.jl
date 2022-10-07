using Test

include("../src/graphs/temporal_graph.jl")
include("../src/centralities/temporal_shortest_betweenness.jl")

@testset "next_temporal_neighbours" begin
    @testset "from small file" begin
        tg::temporal_graph = load_temporal_graph("graphs/test/small.txt", " ")
        tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
        tn::Array{Tuple{Int64,Int64}} = next_temporal_neighbors(tal, 1, 3)
        @test tn == [(2, 4), (2, 5), (2, 6)]
        tn = next_temporal_neighbors(tal, 2, 3)
        @test tn == [(3, 7), (4, 8)]
        tn = next_temporal_neighbors(tal, 3, 12)
        @test tn == []
        tn = next_temporal_neighbors(tal, 4, 10)
        @test tn == []
    end
    @testset "one edge file" begin
        tg::temporal_graph = load_temporal_graph("graphs/test/one_edge.txt", " ")
        tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
        tn::Array{Tuple{Int64,Int64}} = next_temporal_neighbors(tal, 2, 1)
        @test tn == []
    end
end;

@testset "temporal_node_index" begin
    @testset "from small file" begin
        tg::temporal_graph = load_temporal_graph("graphs/test/small.txt", " ")
        tni::Dict{Tuple{Int64,Int64},Int64} = temporal_node_index(tg)
        @test tni[(2, 1)] == 1 && tni[(2, 6)] == 6
        @test tni[(3, 7)] == 7 && tni[(3, 10)] == 10
        @test tni[(4, 8)] == 8 && tni[(4, 9)] == 9
        @test tni[(5, 11)] == 11
        @test get(tni, (2, 10), 0) == 0
    end
    @testset "from self-loop file" begin
        tg::temporal_graph = load_temporal_graph("graphs/test/self_loop.txt", " ")
        tni::Dict{Tuple{Int64,Int64},Int64} = temporal_node_index(tg)
        @test tni[(1, 1)] == 1
        @test get(tni, (1, 2), 0) == 0
        @test get(tni, (2, 1), 0) == 0
    end
    @testset "from one edge file" begin
        tg::temporal_graph = load_temporal_graph("graphs/test/one_edge.txt", " ")
        tni::Dict{Tuple{Int64,Int64},Int64} = temporal_node_index(tg)
        @test tni[(2, 1)] == 1
        @test get(tni, (1, 1), 0) == 0
        @test get(tni, (2, 2), 0) == 0
    end
end;

@testset "temporal_shortest_betweenness" begin
    @testset "from small file" begin
        tg::temporal_graph = load_temporal_graph("graphs/test/small.txt", " ")
        tsb::Array{Float64}, _ = temporal_shortest_betweenness(tg, 0)
        @test isapprox(collect.(tsb), collect.([0.0, 3.0, 3.0, 0.0, 0.0]))
    end
    @testset "from one edge file" begin
        tg::temporal_graph = load_temporal_graph("graphs/test/one_edge.txt", " ")
        tsb::Array{Float64}, _ = temporal_shortest_betweenness(tg, 0)
        @test tsb == [0.0, 0.0]
    end
    @testset "from self-loop file" begin
        tg::temporal_graph = load_temporal_graph("graphs/test/self_loop.txt", " ")
        tsb::Array{Float64}, _ = temporal_shortest_betweenness(tg, 0)
        @test tsb == [0.0]
    end
    @testset "from star file" begin
        tg::temporal_graph = load_temporal_graph("graphs/test/star.txt", " ")
        tsb::Array{Float64}, _ = temporal_shortest_betweenness(tg, 0)
        @test isapprox(collect.(tsb), collect.([0.0, 20.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))
    end
    @testset "from path file" begin
        tg::temporal_graph = load_temporal_graph("graphs/test/path.txt", " ")
        tsb::Array{Float64}, _ = temporal_shortest_betweenness(tg, 0)
        @test isapprox(collect.(tsb), collect.([0.0, 4.0, 6.0, 6.0, 4.0, 0.0]))
    end
    @testset "from lollipop file" begin
        tg::temporal_graph = load_temporal_graph("graphs/test/lollipop.txt", " ")
        tsb::Array{Float64}, _ = temporal_shortest_betweenness(tg, 0)
        @test isapprox(collect.(tsb), collect.([0.0, 0.0, 0.0, 12.0, 12.0, 10.0, 6.0, 0.0]))
    end
    @testset "from multipath file" begin
        tg::temporal_graph = load_temporal_graph("graphs/test/multipath.txt", " ")
        tsb::Array{Float64}, _ = temporal_shortest_betweenness(tg, 0)
        @test isapprox(collect.(tsb), collect.([0.0, 1.0, 1.0, 3.0, 0.0]))
    end
end;
