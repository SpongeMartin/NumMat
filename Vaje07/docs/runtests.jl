using Test, Vaje07

@testset "Newtonova iteracija" begin
    f(x) = x^2 - 2
    df(x) = 2x
    tol = 1e-5
    x, it = newton(f,df,2.0;tol=tol)
    @test f(x).abs < tol
end