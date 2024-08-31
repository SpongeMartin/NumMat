using Main.DN2
#' # Gauss-Legendrove kvadrature
#' Martin Starič
#'
#' Gauss-Legendre kvadratura pravila so namenjena numeričnem integriranju in želijo eksaktno izračunati polinome stopnje $$\leq 2n-1$$.
#' Pri Gauss-Legendre kvadraturah reda 2, želimo aproksimirati integral $$\int_{-1}^1 f(x) dx = w_1f(x_1) + w_2f(x_2)$$, kjer sta $$w_1$$ in $$w_2$$ uteži, $$x_1$$ in $$x_2$$ pa vozla.
#' ## Izpeljava Gauss-Legendre kvadrature reda 2
#' Ta integral izpeljemo na sledeč način:
#' Da določimo uteži $$w_1$$ in $$w_2$$ se osredotočimo na integral polinomov stopnje 1,2,3 in 4.
#'
#' $$\int_{-1}^1 1 dx = 2 = w_1 * 1 + w_2 * 1 = w_1 + w_2$$
#' $$\int_{-1}^1 x dx = 0 = w_1 * x_1 + w_2 * x_2$$
#' $$\int_{-1}^1 x^2 dx = \frac{2}{3} = w_1 * x_1^2 + w_2 * x_2^2$$
#' $$\int_{-1}^1 x^3 dx = 0 = w_1 * x_1^3 + w_2 * x_2^3$$
#'
#' Sedaj če drugo enačbo pomnožimo z $$x_1^2$$ dobimo 
#'
#' $$w_1 * x_1^3 = -w_2x_2x_1^2$$
#'
#' In levi del vstavimo v četrto enačbo dobimo 
#'
#' $$-w_2x_2x_1^2 + w_2 * x_2^3 = 0$$
#'
#' Preoblikujemo jo v $$w_2x_2(x_2^2-x_1^2) = 0$$ in preučimo možnosti za 0.
#' $$w_2 = 0$$ ne velja, ker če si ogledamo drugo in četrto enačbo dobimo protislovje 
#'
#' $$\frac{2}{3} = w_1 * x_1^2 ; w_1 \neq 0, x_1 \neq 0$$ in $$0 = w_1 * x_1^3, w_1 = 0 || x_1 = 0$$.
#'
#' Podobno velja če vzamemo $$x_2 = 0$$
#' Tako nam ostane le še zadnja možnost $$x_2^2 - x_1^2 = 0$$ katero preoblikujemo v $$x_2^2 = x_1^2$$ in iz tretje enačbe dobimo 
#'
#' $$\frac{2}{3} = (w_1 * w_2) * x_1^2$$
#' Iz prve enačbe vemo 
#'
#'$$w_1 * w_2 = 2 \implies x_1 = \pm \frac{1}{\sqrt{3}}$$
#' iz  $$x_2^2 - x_1^2 = 0$$ pa vemo da morata biti $$x_2$$ in $$x_1$$ nasprotno predznačena zato $$x_1 = -\frac{1}{\sqrt{3}}$$ in $$x_2 = \frac{1}{\sqrt{3}}$$.
#' Iz druge enačbe sledi 
#'
#'$$w_1*(-\frac{1}{\sqrt{3}}) + w_2*(\frac{1}{\sqrt{3}}) = 0$$
#' Zato $$w_1 == w_2$$ in iz prve enačbe velja da $$w_1 + w_2 = 2$$ torej sta $$w_1 = w_2 = 1$$.
#' Tako smo izpeljali kvadraturno formulo 
#'
#'$$\int_{-1}^1 f(x) dx = f(-\frac{1}{\sqrt{3}}) + f(\frac{1}{\sqrt{3}})$$
#' Če želimo izpeljati integral z drugimi mejami denimo $$\int_0^h f(x) dx$$ uporabimo linearno preslikavo
#' $$L_2(x) = Ax + B$$, kjer $$L(-1) = 0 = A * (-1) + B$$ in $$L(-1) = h = A * 1 + B$$ iz tega dobimo, da je $$B = \frac{h}{2}$$ in če B vstavimo 
#' v drugo enačbo dobimo $$A = -\frac{h}{2}$$
#' Tedaj rezultat vstavimo 
#'
#' $$\int_0^h f(x) dx = \int_{-1}^1 f(-\frac{h}{2}t + \frac{h}{2}) * - \frac{h}{2}dt$$
#'
#' kjer je
#'
#' $$f(-\frac{h}{2}t + \frac{h}{2}) * - \frac{h}{2} = F(t)$$ 
#'
#' in dobimo aproksimacijo
#'
#' $$F(-\frac{1}{\sqrt{3}})+ F(\frac{1}{\sqrt{3}}) = f(-\frac{h}{2}(-\frac{1}{\sqrt{3}}) + \frac{h}{2}) * - \frac{h}{2} + f(-\frac{h}{2}(\frac{1}{\sqrt{3}}) + \frac{h}{2}) * - \frac{h}{2}$$.
#'
#' Sedaj pa izpeljimo še napako 
#'
#' $$R_f = \frac{(b-a)^{2n+1} * (n!)^4}{(2n+1)[(2n)!]^3} * f^{(2n)}(\epsilon) = $$
#' $$R_f = \frac{(b-a)^5 * 16}{5*24^3} * f^{4}(\epsilon)$$
#'
#' Če želimo uporabiti sestavljeno pravilo, potem moramo integral $$(0,h)$$ preprosto le razdeliti na več delov denimo $$h/n$$ delov, in izračunati integrale z podano kvadraturno formulo in te rezultate sešteti.

#' ## Uporaba izpeljanega pravila v Julia

#' Uporaba pravila za računanje integrala $$sin(x)/x$$ na intervalu $$(0,5)$$
f(x) = sin(x)/x
rezultat = GaussLegendre2(f,5.0,100)

#' Sedaj ročno preverimo napako, četrti odvod funkcije $$f''''(x) = \frac{(x^4 - 12x^2 + 24) * sin(x) + (4x^3 - 24x) * cos(x)}{x^5}$$
f4(x) = ((x^4 - 12x^2 + 24) * sin(x) + (4x^3 - 24x) * cos(x))/ x^5
#' Tedaj izračunajmo koliko približno potrebujemo korakov, da bo integral izračunan na 10 decimalk natančno s pomočjo ocene za napako
tol = 1e-10
GaussLegendre2error(f,f4,5.0,tol)

#' Zgornji rezultat je enak 195, toda temu ni res tako. Oglejmo si rezultat pri točno katerem koraku je natančnost pravila na 10 decimalk.
using QuadGK
result,_ = quadgk(f, 0.0, 5.0)

stevilo_korakov = 2
rezultat = GaussLegendre2(f,5.0,1)
while abs(result - rezultat) - tol > 0
    rezultat = GaussLegendre2(f,5.0,stevilo_korakov)
    stevilo_korakov = stevilo_korakov + 1
end
println(stevilo_korakov)
#' Izkaže se, da je 123 intervalov dovolj.
