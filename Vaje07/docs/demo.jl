#' # Nelinearne enačbe
#' Iščemo rešitev nelinearne enačbe $f(t) = F(x(t),y(t),z(t)) = 0$, kjer je t spremenljivka.
#' Uporabimo Newtonovo metodo.

#' Poiščemo presečišče $y = sin(x)$ in $y = cos(x)$
using Vaje07
using Plots
using ForwardDiff
f(x) = sin(x) - cos(x)
df(x) = cos(x) + sin(x)
x, it = newton(f,df,0.0)

plot(sin,0,2*π)
plot!(cos,0,2*π)
scatter!([x],[sin(x)])

#' Poiščemo presečišče med poltraka in ploskve
x0 = [0,0,0] # Začetna točka poltraka
e = [1,1,0.1] # Smerni vektor poltraka
r = 2
R = 3
F(x) = (R - sqrt(x[1]^2 + x[2]^2))^2 + x[3]^2 - r # Enačba torusa
DF(x) = 2 * vcat((-R/sqrt(x[1]^2 + x[2]^2) - 1) .* x[1:2],x[3])
autoDF(x) = ForwardDiff.gradient(F,x)
p = presecisce((x0,e), F,DF, 1.0)

#' spremenljivka
function z(x,y) 
    z2 = r^2 - (R - sqrt(x^2 + y^2))^2
    if z2 < 0
        return 0
    else
        return sqrt(z2)
    end
end
x = LinRange(-5,5,100)
y = LinRange(0,5,100)
surface(x,y,z)

plot!([x0[1],p[1]],[x0[2],p[2]],[x0[3],p[3]])
