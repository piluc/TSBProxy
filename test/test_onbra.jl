using Test

include("../src/graphs/temporal_graph.jl")
include("../src/centralities/temporal_shortest_betweenness.jl")
include("../src/centralities/onbra.jl")

@testset "onbra_sample" begin
    @testset "from small file" begin
        tg::temporal_graph = load_temporal_graph("graphs/test/small.txt", " ")
        sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 3)
        @test length(sample) == 3
        for i in 1:3
            @test sample[i][1] != sample[i][2]
        end
        sample = onbra_sample(tg, 20)
        @test length(sample) == 20
        for i in 1:20
            @test sample[i][1] != sample[i][2]
        end
    end
end;

@testset "onbra_sample" begin
    @testset "from small file" begin
        tg::temporal_graph = load_temporal_graph("graphs/test/small.txt", " ")
        o::Array{Float64}, _ = onbra_sample(tg, 10, 0, test_sample=[(3, 5), (2, 1), (5, 3), (5, 1), (2, 5), (1, 3), (1, 3), (5, 1), (2, 4), (1, 2)])
        @test isapprox(collect.(o), collect.([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))
    end
end;
