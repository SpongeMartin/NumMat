

#' Potenčna metoda

#' Iščemo lastni vektor matrike matrike A za največjo lastno vrednost.

#' ## Invariantna mera za Markovsko verigo

P = [0.1 0.4 0.5;
     0 0.5 0.5;
     0.2 0 0.8]

using Vaje06

#' Invariantna mera je lastni vektor za $P^T$ za lastno vrednost 1.

p, lambda = potencna(P', ones(3))
p = p/sum(p)

#' Invariantna mera za Markovsko verigo na 6 vozliščih. (Bipartitni graf)

P = [0 0.3 0 0.4 0 0.3;
     0.1 0 0.2 0 0.7 0;
     0 0.5 0 0.2 0 0.3;
     0.4 0 0.2 0 0.4 0;
     0 0.5 0 0.2 0 0.3;
     0.4 0 0.2 0 0.4 0]

using LinearAlgebra
eigen(P')

p, lambda = potencna(P', rand(6))
p = p/sum(p)

#' Premik $A - δI$
#' Če so $λ_1,...,λ_n$ lastne vrednosti matrike A, potem so $λ_1-δ,...,λ_n-δ$ lastne vrednosti matrike $A - δI$

#' Če ima matrika A več različnih lastnih vrednosti, ki so po absolutni vrednosti največje, potem potenčna metoda ne konvergira za vsak
#' Začetni približek. V prejšnem primeru je imela matrika lastne vrednosti 1 in -1. Problem rešimo s premikom.

B = P + I

p,lambda = potencna(B',rand(6))

p = p/sum(p)

using Plots

scatter(p,title="Invariantna porazdelitev")