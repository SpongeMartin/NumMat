#' # QR razcep zgornje hessenbergove matrike
#' Martin Starič
#'
#' QR razcep zgornje Hessenbergove matrike `H` je najbolj primeren s pomočjo Givensovih rotacij, saj zanj potrebujemo najmanjše število operacij. Givensove rotacije temeljijo na podlagi
#' rotacije vektorja $$v \in R^2$$ s pomočjo rotacijske matrike $$R = G(i,j,\theta)$$. Ko rotacijsko matriko `R_i` pomnožimo z $$H_i-1$$ in dobimo $$H_i$$ (kjer i označuje v kateri iteraciji smo), bo drugi element, ki smo ga izbrali za $$v$$
#' uničen ($$= 0$$). Ta postopek iteriramo dokler vseh poddiagonalnih elementov ne uničimo. QR razcep pa je nakoncu takšen, da je $$H_n = R$$ in zmnožek $$R^T_1 * R^T_2 * ... * R^T_n = Q$$

using Main.DN1
#' Sprva izberimo zgornje Hessenbergovo matriko 
H = DN1.Hessenberg([2 7 2;
                3 2 7;
                0 1 4])

#' Nato izračunamo njen QR razcep s pomočjo Givensovih rotacij, saj je zaradi Hessebergove oblike uporaba Givensovih rotacij najbolj primerna in učinkovita metoda.
Q,R = DN1.qr(H)

#' ## QR iteracija
#' QR iteracija je postopek, kjer matriko A spreminjamo v zgornje trikotno tako, da izračunamo $$A = QR$$ in jo zamenjamo z $$A = RQ$$, to ponavljamo več iteracij in nazadnje preberemo
#' diagonalne elemente matrike `A` - ki so lastne vrednosti matrike `A`. Sledi demonstracija te metode
lastnevrednosti,lastnivektorji = DN1.eigen(H,1000)
#' Da rezultat primerjamo, si pomagamo s knjižnico LinearAlgebra, ki ima funkciji za izračun qr razcepa in lastnih vrednosti
using LinearAlgebra
mt =  [2 7 2;
        3 2 7;
        0 1 4]
QR = LinearAlgebra.qr(mt)
podatki = LinearAlgebra.eigen(mt)
podatki.values

#' Da preverimo čas potreben za izvedbo implementiranih funkcij si pomagamo z dekoratorjem @timed, oglejmo si časovno potratnost metode qr
dn_QR_Result = (@timed DN1.qr(H))
LA_QR_Result = (@timed LinearAlgebra.qr(mt))
dn_QR_Result.time - LA_QR_Result.time
#' Opazimo, da je naša funkcija malenkost hitrejša
#' Sedaj preverimo še rezultate za računanje lastnih vrednosti
dn_eigen_Result = (@timed DN1.eigen(H,10))
LA_eigen_Result = (@timed LinearAlgebra.eigen(mt))
dn_eigen_Result.time - LA_eigen_Result.time

dn_eigen_Result2 = (@timed DN1.eigen(H,100))
dn_eigen_Result2.time - LA_eigen_Result.time
#' Izkaže se, da je vse odvisno od željene natančnosti, v kolikor želimo višjo natančnost metode, je naša implementirana funkcija z qr iteracijo počasnejša