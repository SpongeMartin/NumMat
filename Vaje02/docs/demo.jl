#' # Vaja 02 - linearni sistemi

#' Reši sistem linearnih enačb
#' ```
#' x + 2y - 7 = 1
#' 2x - y - 3z = 2
#' x - z = 3
#' ```
#' Če želimo sistem rešit z računalnikom moramo sistem spremeniti v matrično obliko.
#' $Ax = b$.
#' Formalno lahko zapišemo rešitev $x = A^{-1} * b$, toda imamo boljše alternative.
#' alternative: reši sistem z LU razcepom. $A = L*U$, kjer L je spodnje trikotna matrika, U pa zgornje trikotna.
#' \\ v juliji je operator za deljenje z desne. Torej $A / b = A^{-1} * b$
#' $b^T / A = b^T * A^{-1}$


A = [1 2 -1; 2 -1 -3; 1 0 -1.0] # matrika sistema

b = [1,2,3] #desne strani

x = A\b # rešitev sistema $A x = b (A / b = A^{-1} b)$

#' Preizkus naredimo tako, da pomnožimo $Ax - b$ enak 0.

A * x - b

#' Sistem lahko rešimo tudi z LU razcepom, tako da sistem $LUx = b$ prevedemo
#' na dva sistema $Ly=b$ in $Ux = y$.

using LinearAlgebra

#' p je permutacijski vektor, ki premeša sisteme, tako da optimalno razporedi sisteme za minimiziranje napak.

L,U,p = lu(A)

x = U \ (L \ b[p])
norm(A*x-b, Inf)

#' Julia faktorje razcepa zapakira v en objekt.

F = lu(A)
#' rezultat je LU, ki je faktorizacija.
x = F \ b

#' ## Tridiagonalne matrike
using Vaje02
T = Tridiag([1,2],[3,4,5],[6,7])

#' Produkt matrike  z vektorjem 
T * [1,2,3]

#' Deljenje z leve oziroma operator '\'


#' # Slučajni sprehod
#' `X_{n} = \sum{n,i=1}Bin(p)`
#' Bin(p) ˜ [[1,-1],[p,1-p]]
"Generator naključnih števil porazdeljenih kot Bin(p)"
rand_bin(p) = (rand() < p) ? 1 : -1
[rand_bin(0.6) for i = 1:10]

sprehod(p,n) = cumsum([rand_bin(p) for i=1:n])

sprehod(0.6,10)

using Plots
scatter(sprehod(0.5,100))
#' # Markovske verige $\{X_{n}\}$
#' $P(X_{n}|X_{n-1}=x_{n-1},X_{n-2}=x_{n-2},...,X_{1}) = x_{1} = P(X_{n-1}) = x_{n-1}$

n = 100
p = 0.495
T = Tridiag(-p*ones(2*n-2),ones(2*n-1),-(1-p)*ones(2*n-2))
k = T\ones(2*n-1)
scatter(-n:n,k)