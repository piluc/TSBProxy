using Test

include("../src/statistics/rankings.jl")

@testset "min_h_k" begin
    @testset "from small file with exact and ptd" begin
        a::Array{Float64} = [0.0, 3.0, 3.0, 0.0, 0.0]
        b::Array{Float64} = [0.0, 2.0, 3.0, 2.0, 0.0]
        @test min_h_k(a, b, 1) == 2
        @test min_h_k(a, b, 2) == 2
        @test min_h_k(a, b, 3) == 4
        @test min_h_k(a, b, 4) == 4
        @test min_h_k(a, b, 5) == 5
    end

    @testset "from small file with exact and exact" begin
        a::Array{Float64} = [0.0, 3.0, 3.0, 0.0, 0.0]
        b::Array{Float64} = [0.0, 3.0, 3.0, 0.0, 0.0]
        @test min_h_k(a, b, 1) == 1
        @test min_h_k(a, b, 2) == 2
        @test min_h_k(a, b, 3) == 3
        @test min_h_k(a, b, 4) == 4
        @test min_h_k(a, b, 5) == 5
    end

    @testset "sorted versus all equal" begin
        a::Array{Float64} = [5.0, 4.0, 3.0, 2.0, 1.0]
        b::Array{Float64} = [1.0, 1.0, 1.0, 1.0, 1.0]
        @test min_h_k(a, b, 1) == 1
        @test min_h_k(a, b, 2) == 2
        @test min_h_k(a, b, 3) == 3
        @test min_h_k(a, b, 4) == 4
        @test min_h_k(a, b, 5) == 5
    end

    @testset "inverse sorted versus all equal" begin
        a::Array{Float64} = [1.0, 2.0, 3.0, 4.0, 5.0]
        b::Array{Float64} = [1.0, 1.0, 1.0, 1.0, 1.0]
        @test min_h_k(a, b, 1) == 5
        @test min_h_k(a, b, 2) == 5
        @test min_h_k(a, b, 3) == 5
        @test min_h_k(a, b, 4) == 5
        @test min_h_k(a, b, 5) == 5
    end
end;
