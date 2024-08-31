module Vaje08

using LinearAlgebra, SparseArrays
export inverse_iter, inverse_iterqr, graf_eps, laplace

"""
    v, lambda = inverse_iter(f, n)

Poišči lastno vrednost `lambda` in lastni vektor `v` za matriko z inverzno iteracijo.
Argument `f` je funkcija, ki za dani vektor `b` poišče rešitev sistema `Ax=b`.
Argument `n` je dimenzija matrike.

# Primer

```jl
A = [3 1 1; 1 1 3; 1 3 4]
F = lu(A)
f(b) = F \\ b
v, lambda = inverse_iter(f, 3)
```
"""
function inverse_iter(f, n, maxit=100, tol=1e-10)
    x = rand(n)
    ls, index = findmax(abs, x)
    ls = x[index]
    x = x / ls # normirano, da je maksimalni element enak 1
    for i = 1:maxit
        x = f(x) # poišči rešitev sistema Ay = x
        ln, index = findmax(abs, x)
        ln = x[index]
        x = x / ln # x normiramo, da je maksimalni element enak 1
        if abs(ln - ls) < tol
           println("Inverzna potenčna metoda se je končala po $i korakih.")
            return x, 1 / ln
        end
        ls = ln
    end
    throw("Potenčna metoda ne konvergira v $maxit korakih.")
end

"""
    x, lambda = inverse_iterqr(f, n, k)

Poišči najmanjših m lastnih vektorjev in lastnih vrednosti z inverzno iteracijo s QR 
razcepom.
"""
function inverse_iterqr(f, n, k, maxit=100, tol=1e-10)
    x0 = rand(n, k)
    x = x0
    _, m = size(x0)
    ls = Inf * ones(m)
    for i = 1:maxit
        x = f(x) # rešitev Ay = x
        Q, R = qr(x)
        x = Matrix(Q)
        ln = diag(R) # lastne vrednosti A^(-1) na diagonali R
        if norm(ln - ls, Inf) < tol
            println("Inverzna potenčna metoda se je končala po $i korakih.")
            return x, 1 ./ ln
        end
        ls = ln
    end
    throw("Inverzna potenčna metoda ne konvergira v $maxit korakih.")
end

"""
    A = graf_eps(oblak, epsilon)

Poišči podobnostni graf ε-okolic za dani oblak točk.
Argument `oblak` je `k x n` matrika, katere stolpci so koordinate točk v oblaku.
Funkcija vrne matriko sosednosti A za podobnostni graf.
"""
function graf_eps(oblak, epsilon)
    _,n = size(oblak)
    A = SparseArrays.spzeros(n, n)
    for i = 1:n
        for j = i+1:n-1
            if LinearAlgebra.norm(oblak[:, i] - oblak[:, j],3) < epsilon
                A[i, j] = 1
                A[j, i] = 1
            end
        end
    end
    return A
end

"""
  L = laplace(A)

Poišči Laplaceovo matriko za dani graf. Laplaceova matrika grafa ima na diagonali stopnje točk,
izven diagonale pa vrednosti -1 ali 0, odvisno ali sta indeksa povezana v grafu. Graf je podan z
matriko sosednosti `A`.
"""
laplace(A) = SparseArrays.spdiagm(sum(A, dims=2)[:, 1]) - A

end