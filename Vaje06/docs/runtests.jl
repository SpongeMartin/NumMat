using Test, Vaje06

@testset "3x3 prehodna matrika" begin
    #lastni vektor je [1,1,1] za lastno vrednost 1
    P = [0.1 0.4 0.5; 0 0.5 0.5; 0.2 0 0.8]
    v,lambda = potencna(P,[-1,0,0])
    @test isapprox(lambda,1)
    @test isapprox(v,[1,1,1])
end