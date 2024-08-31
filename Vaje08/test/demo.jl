#' ## Spektralno razvrščanje v gruče
#' Spektralno razvrščanje v gruče uporabi Laplaceovo matriko (`L`) podobnostnega grafa podatkov, da podatke preslika v prostor kjer jih je lažje razvrstiti.
#' Algoritem deluje tako, da najprej poiščemo k najmanjših lastnih vrednosti za `L`.
#' `Q = [v_1,v_2,...,v_k]` je matrika lastnih vektorjev.
#' Za stolpce `Q^T` izvedemo algoritem k-povprečij.

using Vaje08 # To so pravzaprav vaje06, napačno sem jih že prej generiral.

#' Sprva generiramo oblak točk v ravnini
using Plots
using Random
m = 100;
Random.seed!(12)
x = [1 .+ randn(m, 1); -3 .+ randn(m,1); randn(m,1)];
y = [-2 .+ randn(m, 1); -1 .+ randn(m,1); 1 .+ randn(m,1)];
scatter(x, y, title="Oblak točk v ravnini")
savefig("1.png")
#' Tedaj ustvarimo Laplaceovo matriko na podlagi inverzne potenčne metode, ki se uporablja zato da dobimo le nekaj najnižjih lastnih vrednosti.
using SparseArrays
tocke = hcat(x,y)'
r = 1.4
G = graf_eps(tocke, r)
L = laplace(G)
spy(sparse(Matrix(L)), title="Porazdelitev neničelnih elementov v laplaceovi matriki")
savefig("2.png")
#' Izračunamo prvih 20 lastnih vrednosti
import LinearAlgebra.eigen
razcep = eigen(Matrix(L))
scatter(razcep.values[1:20], title="Prvih 20 lastnih vrednosti laplaceove matrike")
savefig("3.png")
scatter(razcep.vectors[:,4], razcep.vectors[:,5], title="Vložitev s komponentami 4. in 5. lastnega vektorja")
savefig("4.png")
#' Zgleda, kot da je matrika L napačno generirana saj graf vložitvje 4. in 5. lastnega vektorja ne izgleda pravilno

#' Sedaj pa uporabimo še clustering funkcijo kmeans, da izračunamo clusterje na podlagi podobnostnega grafa.
using Clustering
nove_tocke = hcat(razcep.vectors[:,4],razcep.vectors[:,5])'
gruce = kmeans(nove_tocke, 3).assignments

#' Originalen graf
p1 = scatter(tocke[findall(gruce .== 1)], color=:blue, title="Originalne točke")
scatter!(p1, tocke[findall(gruce .== 2)], color=:red)
scatter!(p1, tocke[findall(gruce .== 3)], color=:green)

#' Graf podobnosti
p2 = scatter(nove_tocke[findall(gruce .== 1)], color=:blue, title="Preslikane točke")
scatter!(p2, nove_tocke[findall(gruce .== 2)], color=:red)
scatter!(p2, nove_tocke[findall(gruce .== 3)], color=:green)

plot(p1,p2)
savefig("5.png")