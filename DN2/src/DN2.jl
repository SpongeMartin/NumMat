module DN2

export GaussLegendre2

"""
Gauss-Legendrove kvadraturno pravilo na dveh točkah
"""
function GaussLegendre2(f,h::Float64,subintegrals::Int)
    x0 = -1/sqrt(3)
    x1 = 1/sqrt(3)
    value = 0.0
    h = h/subintegrals
    b = h
    a = 0.0
    for _ in 1:subintegrals
        t0 = (b-a)/2 * x0 + (b+a)/2
        t1 = (b-a)/2 * x1 + (b+a)/2
        value += (b-a)/2 * (f(t0) + f(t1))
        a += h
        b += h
    end
    
    return value
end

end # module DN2
