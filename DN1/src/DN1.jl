module DN1
using LinearAlgebra
export ZgornjiHessenberg,qr,Givens,Hessenberg,qrIteration

import Base: *

"""
Podatkovna struktura za 3diagonalno matriko.
"""
#Uporabi ] activate Vaje02 in pol using backspace Vaje02

struct ZgornjiHessenberg
  H::Matrix{Float64}
end

struct Givens
    G::Matrix{Float64}
end


function *(H::ZgornjiHessenberg,G::Matrix)
    return G.G * H.H 
end

"""
Ustvari Givens matriko, ki ima 2 vrstici n stolpcev, vsak stolpec predstavlja rotacijo [cos(a),sin(a)]^T
"""
function Givens(n::Int)
    return Givens(zeros(2,n))
end

"""
Ustvari Zgornjo Hessenberg matriko
"""
function Hessenberg(M::Matrix)
    return ZgornjiHessenberg(M)
end
"""
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
QR razcep Zgornje Hessebergove matrike z Givensovimi rotacijami.
"""
function qr(Zh::ZgornjiHessenberg)
    # TODO: QR razcep Hessebergove matrike z Givensovimi rotacijami.
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

function qrIteration(H::ZgornjiHessenberg,iter = 10)
    n = size(H.H,1)
    R = Matrix{Float64}(I,n,n)
    for _ in 1:iter
        G,R = qr(H)
        Q = GivensRotationMatrix(1,G,n)
        for i in 2:n-1
            Q = Q * GivensRotationMatrix(i,G,n)
        end
        H = Hessenberg(R*Q)
    end
    return R
end

end # module DN1
