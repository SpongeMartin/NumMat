module Vaje07

export newton,presecisce

"""
    x,it = newton(F,DF,x0;maxit=maxit,tol=tol)

Poišči funkcijo približek za rešitev enačbe F(x)=0 z Newtonovo metodo za dano 
funkcijo `F`, odvodom funkcije `DF` in začetnim približkom x0
"""

function newton(F,DF,x0,maxit=100,tol=1e-10)
    for i=1:maxit
        z = F(x0)
        x = x0 - F(x0)/DF(x0)
        if abs(z) < tol
            return x,i
        end
        x0 = x
    end
    throw("Metoda ne konvergira")
end

"""
    T = presecisce(poltrak,F,t0)

Poišči presečišče poltraka in implicitno podane ploskve z enačbo F(x,y,z) = 0
z Newtonovo metodo z začetnim približkom za parameter na poltraku `t0`
DF je funkcija odvoda(gradient) funkcije F
"""

function presecisce(poltrak, F, DF, t0)
    x0, e = poltrak
    # funkcija desnih strani
    x(t) = x0 + t*e
    f(t) = F(x(t))
    df(t) = DF(x(t))' * e # skalarni produkt
    t, it = newton(f,df,t0)
    if t < 0
        throw("Presečišče ni na poltraku.")
    end
    return x(t) # presečišče
end

end # module Vaje07
