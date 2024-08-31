module DN1
using LinearAlgebra
export ZgornjiHessenberg,qr,Givens,Hessenberg,eigen

import Base: *

struct ZgornjiHessenberg
  H::Matrix{Float64}
end

struct Givens
    G::Matrix{Float64}
end


function *(H::ZgornjiHessenberg,G::Matrix)
    return G.G * H.H 
end

function *(R::Matrix,G::Givens)
    return R * G.G
end
"""
    G = Givens(M)

Ustvari Givens matriko, ki ima 2 vrstici n stolpcev, vsak stolpec predstavlja rotacijo [cos(a),sin(a)]^T
"""
function Givens(n::Int)
    return Givens(zeros(2,n))
end

"""
    H = Hessenberg(M)

Ustvari Zgornjo Hessenberg matriko
"""
function Hessenberg(M::Matrix)
    return ZgornjiHessenberg(M)
end
"""
    R = GivensRotationMatrix(i,G,n)

Ustvari Givensovo rotacijsko matriko G(i,j,a)
"""
function GivensRotationMatrix(i::Int,G::Givens,n::Int)
    H = Matrix{Float64}(I,n,n)
    H[i,i] = G.G[1,i]
    H[i+1,i+1] = G.G[1,i]
    H[i+1,i] = -G.G[2,i]
    H[i,i+1] = G.G[2,i]
    return H
end


"""
    Q,R = qr(Zh)

QR razcep Zgornje Hessebergove matrike z Givensovimi rotacijami.
"""
function qr(Zh::ZgornjiHessenberg)
    velikost = size(Zh.H, 1)
    Gi = Givens(velikost-1)
    i = 1
    R = copy(Zh.H)
    while i < velikost
        Gi.G[1,i] = 1.0 #c = 1, s = 0 če želimo da je Givensova rotacija G(i,j,a) = I (indentiteta)
        if R[i+1,i] != 0
            a = R[i,i]
            b = R[i+1,i]
            c = sqrt(a^2 + b^2)
            Gi.G[1,i] = a/c #c
            Gi.G[2,i] = -b/c #s
            R = GivensRotationMatrix(i,Gi,velikost)' * R
        end
        i = i + 1
    end
    return Gi,R
end

"""
    evalues,evectors = eigen(H,iter)

Funkcija uporabi metodo QR iteracije za računanje lastnih vrednosti in lastnih vektorjev matrike H v iter iteracijah.
"""
function eigen(H::ZgornjiHessenberg,iter = 10)
    n = size(H.H,1)
    R = Matrix{Float64}(I, n, n)
    lastnivektorji = Matrix{Float64}(I, n, n)
    for _ in 1:iter
        G,R = qr(H)
        Q = GivensRotationMatrix(1,G,n)
        for i in 2:n-1
            Q = Q * GivensRotationMatrix(i,G,n)
        end
        H = Hessenberg(R*Q)
        lastnivektorji *= Q
    end
    lastnevrednosti = LinearAlgebra.diag(H.H)
    return lastnevrednosti, lastnivektorji
end

end # module DN1
