using Main.DN1
using Test
using LinearAlgebra

@testset "QR dekompozicija" begin
    H = [4.0 3.0 2.0 1.0;
    1.0 4.0 3.0 2.0;
    0.0 1.0 4.0 3.0;
    0.0 0.0 1.0 4.0]
    He = DN1.Hessenberg(H)
    _,R = DN1.qr(He)
    LA = LinearAlgebra.qr(H)
    @test abs.(R) ≈ abs.(LA.R) #abs ker se predznaki razlikujejo
end
@testset "Lastne vrednosti" begin
    H = [4.0 3.0 2.0 1.0;
    1.0 4.0 3.0 2.0;
    0.0 1.0 4.0 3.0;
    0.0 0.0 1.0 4.0]
    He = DN1.Hessenberg(H)
    lastneVrednosti,_ = DN1.eigen(He,100)
    LA_lV = LinearAlgebra.eigen(H)
    @test sort(LA_lV.values,rev = true) ≈ sort(lastneVrednosti,rev = true) # Potrebno sortirati, ker lastne vrednosti niso enako razvrščene.
end