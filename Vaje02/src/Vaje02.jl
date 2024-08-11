module Vaje02

export Tridiag

import Base: getindex, *, size, \, setindex!

"""
Podatkovna struktura za 3diagonalno matriko.
"""
#Uporabi ] activate Vaje02 in pol using backspace Vaje02

struct Tridiag
  sp # Vektor elementov na spodnji diagonali
  d # Vektor elementov na diagonali
  zg # Vektor elementov na zgornji diagonali
end

"""
Vrni dimenzije tridiagonalne matrike
"""
size(T::Tridiag) = (length(T.d),length(T.d))

"""
  elt = T[i,j]

Vrni element v i-ti vrstici in j-tem stolpcu tridiagonalne matrike `T`
"""

function getindex(T::Tridiag,i,j)
  # TODO: preveri če sta i in j v obsegu veljavnih indeksov
  if i == j-1
    return T.zg[i]
  elseif i == j
    return T.d[i]
  elseif i == j + 1
    return T.sp[j]
  else
    return zero(T.d[1])
  end
end
"""
Nastavi vrednost elementa na mestu i,j
"""

function setindex!(T::Tridiag,x,i,j)
  # TODO: preveri če sta i in j v obsegu veljavnih indeksov
  n,_m = size(T)
  if (i < 1) || (i > n) || (j < 1) || (j > n)
    throw("Indekso so izven obsega matrike")
  end
  if i == j - 1
    T.zg[i] = x
  elseif i == j
    T.d[i] = x
  elseif i == j + 1
    T.sp[j] = x
  else
    throw("Elementa ni mogoče spremeniti")
  end
end

"""
Izračunaj produkt tridiagonalne matrike z vektorjem.
"""

function *(T::Tridiag,x::Vector)
  y = zero(x)
  n = length(T.d)
  y[1] = T[1,1] * x[1] + T[1,2] * x[2]
  for i=2:n-1
    y[i] = T[i,i-1] * x[i-1] + T[i,i] * x[i] + T[i,i+1] * x[i+1]
  end
  y[n] = T[n,n-1] * x[n-1] + T[n,n] * x[n]
  return y
end

#= function \(T::Tridiag,b::Vector)
  d = copy(T.d)
  y = copy(b)
  n = length(d)
  for i=2:n-1
    l = T[i,i-1]/d[i-1]
    d[i] = d[i] - l * d[i-1]
    b[i] = b[i] - l * b[i-1]
  end
end =#

function \(T::Tridiag,b::Vector)
  n,_m = size(T)
  T = Tridiag(T.sp,copy(T.d),T.zg)
  b = copy(b)
  for i=2:n
    l = T[i,i-1] / T[i-1,i-1]
    T[i,i] = T[i,i] - l * T[i-1,i]
    b[i] = b[i] - l * b[i-1]
  end
  # obratno vstavljanje
  x=zero(b)
  x[n] = b[n]/T[n,n]
  for i = (n-1):-1:1
    x[i] = (b[i] - T[i,i+1]*x[i+1]) / T[i,i]
  end
  return x
end


end