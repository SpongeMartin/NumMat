#' # Robni problem za Laplaceovo enačbo
using Vaje04

rp = Vaje04.RobniProblemPravokotnikLaplace(
    [0,pi,-1,1],
    [sin, sin ,x->0.0 ,x->0.0]
)

#' funkcija resi ne deluje kot pričakovano

Z = resi(rp,0.1)


using Plots
surface(Z)

#' Narisana ploskev

#' Napolnitev Laplaceove matrike

L = laplaceova_matrika(10,10)
spy(L)

using LinearAlgebra

#' Obzervacija - na nediagonalnih mestih neničelni elementi

F = lu(L)
spy(F.L)

spy(F.U)

#' # Iterativne metode

Z = resi_iter(rp,0.1,1) #Gauss-Seidlova iteracija
surface(Z)

#' Naslednji graf

resi_iter(rp,0.1,1.5)