module DN3

export nihalo

"""
Računanje kotnega odmika s pomočjo metode Runge-Kutta 4. stopnje.
"""
function nihalo(l,t,theta0,dtheta0,n)
    g = 9.80665
    f(y) = -g/l * sin(y)
    y0 = [theta0,dtheta0]
    h = t/n
    for i in 0:n-1
        y0 = rk4(t,i*h,y0,f)
    end
    return y0
end

"""
Runge-Kutta 4. stopnje.
"""
function rk4(t,h,y,f)
    k1 = f(t,y)
    k2 = f(t + h/2,y + k1/2)
    k3 = f(t + h/2,y + k2/2)
    k4 = f(t + h,y + k3)
    return y + (k1 + 2*k2 + 2*k3 + k4)/6
end

end # module DN3
