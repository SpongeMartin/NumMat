module DN2

export GaussLegendre2, GaussLegendre2error

using Roots,ForwardDiff

"""
    GaussLegendre2(f,h,subintegrals)

Sestavljeno Gauss-Legendrovo kvadraturno pravilo na dveh točkah za integrale od 0 do h.
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

"""
    maximumf(f4,a,h)

Pomožna funkcija za izračunanje maksimuma petega odvoda neke funkcije s podanim 4. odvodom.
"""
function maximumf(f4,a,h)
    critical_points = Roots.find_zeros(f4, 0, 5)

    points = [a,a+0.01,h-0.01, h, critical_points...]

    values = map(x -> ForwardDiff.derivative(f4, x), points)

    max_value = maximum(values)

    return max_value
end
"""
    GaussLegendre2error(f,f4,h,tol = 1^-10)

Sestavljeno Gauss-Legadrevo kvadraturno pravilo, kjer povečujemo število korakov glede na to ali je zadoščeno določeni toleranci napake.
"""
function GaussLegendre2error(f,f4,h::Float64,tol = 1^-10)
    x0 = -1/sqrt(3)
    x1 = 1/sqrt(3)
    value = 0.0
    i = 1
    a = 0.0
    error = 1
    while error > tol
        error = 0
        i = i + 1
        h0 = h/i
        b = h0
        a = 0.0
        for _ in 1:i
            t0 = (b-a)/2 * x0 + (b+a)/2
            t1 = (b-a)/2 * x1 + (b+a)/2
            value += (b-a)/2 * (f(t0) + f(t1))
            a += h0
            b += h0
            error = error + abs((b - a)^5 * 16 / (5*24^3) * f4(maximumf(f4,a,h)))
        end
    end
    return i
end

end # module DN2