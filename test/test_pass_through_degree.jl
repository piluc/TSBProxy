using Test

include("../src/graphs/temporal_graph.jl")
include("../src/centralities/pass_through_degree.jl")

@testset "edge_time_range" begin
    @testset "from small file" begin
        tg::temporal_graph = load_temporal_graph("test/graphs/small.txt", " ")
        te_min::Vector{Pair{Tuple{Int64,Int64},Int64}}, te_max::Vector{Pair{Tuple{Int64,Int64},Int64}} = edge_time_range(tg)
        @test te_min[1] == ((1, 2) => 1) && te_max[1] == ((1, 2) => 6)
        @test te_min[2] == ((2, 3) => 7) && te_max[2] == ((2, 3) => 7)
        @test te_min[3] == ((2, 4) => 8) && te_max[3] == ((2, 4) => 8)
        @test te_min[4] == ((3, 4) => 9) && te_max[4] == ((3, 4) => 9)
        @test te_min[5] == ((4, 3) => 10) && te_max[5] == ((4, 3) => 10)
        @test te_min[6] == ((3, 5) => 11) && te_max[6] == ((3, 5) => 11)
        @test length(te_min) == 6 && length(te_max) == 6
    end
    @testset "from one edge file" begin
        tg::temporal_graph = load_temporal_graph("test/graphs/one_edge.txt", " ")
        te_min::Vector{Pair{Tuple{Int64,Int64},Int64}}, te_max::Vector{Pair{Tuple{Int64,Int64},Int64}} = edge_time_range(tg)
        @test te_min[1] == ((1, 2) => 1) && te_max[1] == ((1, 2) => 1)
        @test length(te_min) == 1 && length(te_max) == 1
    end
    @testset "from self_loop file" begin
        tg::temporal_graph = load_temporal_graph("test/graphs/self_loop.txt", " ")
        te_min::Vector{Pair{Tuple{Int64,Int64},Int64}}, te_max::Vector{Pair{Tuple{Int64,Int64},Int64}} = edge_time_range(tg)
        @test te_min[1] == ((1, 1) => 1) && te_max[1] == ((1, 1) => 1)
        @test length(te_min) == 1 && length(te_max) == 1
    end
end;

@testset "count_less_than" begin
    a::Array{Int64} = [3, 5, 7, 9, 11]
    @test count_less_than(a, 2) == 0
    @test count_less_than(a, 12) == 5
    @test count_less_than(a, 3) == 0
    @test count_less_than(a, 11) == 4
    @test count_less_than(a, 6) == 2
    @test count_less_than(a, 7) == 2
    @test count_less_than(a, 8) == 3
end;

@testset "pass_through_degree" begin
    @testset "from small file" begin
        tg::temporal_graph = load_temporal_graph("test/graphs/small.txt", " ")
        ptd::Array{Float64}, _ = pass_through_degree(tg)
        @test isapprox(collect.(ptd), collect.([0.0, 2.0, 3.0, 2.0, 0.0]))
    end
    @testset "from one edge file" begin
        tg::temporal_graph = load_temporal_graph("test/graphs/one_edge.txt", " ")
        ptd::Array{Float64}, _ = pass_through_degree(tg)
        @test isapprox(collect.(ptd), collect.([0.0, 0.0]))
    end
    @testset "from self_loop file" begin
        tg::temporal_graph = load_temporal_graph("test/graphs/self_loop.txt", " ")
        ptd::Array{Float64}, _ = pass_through_degree(tg)
        @test isapprox(collect.(ptd), collect.([0.0]))
    end
end;
