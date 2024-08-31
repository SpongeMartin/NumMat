using Main.DN2
using Test
using QuadGK

@testset "Računanje integrala" begin
    p(x) = 2*x^3-5*x^2+3*x-7
    result1,_ = QuadGK.quadgk(p, 0.0, 3.0)
    rezultat1 = DN2.GaussLegendre2(p,3.0,100)
    @test result1 ≈ rezultat1
    g(x) = 3*x^4 - 2*x^3 + π*x^2 - sin(x) + log(x + 1)
    result2,_ = QuadGK.quadgk(g, 0.0, 13.0)
    rezultat2 = DN2.GaussLegendre2(g,13.0,100)
    @test result2 ≈ rezultat2
    h(x) = sin(x) * (x^4 - 5*x^3 + 6*x^2) + (exp(-x) / (x^2 + 1)) * cos(2*x) - (1 / (1 + x^2))
    result3,_ = QuadGK.quadgk(h, 0.0, 2.7)
    rezultat3 = DN2.GaussLegendre2(h,2.7,120)
    @test result3 ≈ rezultat3 
end