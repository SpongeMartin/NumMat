#' # Aproksimacija z linearnim modelom
#' Podatke ${(x_1,y_1),(x_2,y_2),...,(x_n,y_n)}$ želimo opisati s funkcijo oblike $y = f(x,p)$
#' Delali bomo s podatki CO2 na Mauna Loa
#' Potrebno je prebrati podatke iz podatkovne baze in jih pripraviti za delo.
using FTPClient
using Statistics
url = "ftp://aftp.cmdl.noaa.gov/products/trends/co2/co2_weekly_mlo.txt"
io = download(url)
data = readlines(io)

using Plots
filter!(l->l[1]!='#', data)
data = strip.(data)
data = [split(line, r"\s+") for line in data]
data = [[parse(Float64, x) for x in line] for line in data]
filter!(l->l[5]>0, data)
data_4 = [l[4] for l in data]
t0 = mean(data_4)
t1 = [l[4] for l in data]
t = [l[4] - t0 for l in data] 
co2 = [l[5] for l in data]

#' Graf 1
scatter(t1, co2, title="Atmosferski CO2 na Mauna Loa",
        xlabel="leto", ylabel="parts per milion (ppm)", label="co2",
        markersize=1)

#' Glede na podatke izračunamo letne spremembe CO2
using LinearAlgebra
A = hcat(ones(size(t)), t, t.^2, cos.(2pi*t), sin.(2pi*t))
p = A\co2 #co2 = y
norm(A*p-co2)

#' Graf 2
plot([A[:,1]/norm(A[:,1]), A[:,2]/norm(A[:,2]), A[:,3]/norm(A[:,3])], ylims=[0,0.025], title="Normirani stolpci matrike A")

#' Graf 3
plot(t, p[2].+2*p[3]*(t), title="Letne spremembe CO2", xaxis=[1980,1990,2000,2010,2020])

model(t) =  p[1]+ p[2]*(t-t0) + p[3]*(t-t0)^2 + p[4]*cos(2*pi*t) + p[5]*sin(2*pi*t)
model.([2020, 2030, 2050])

#' Graf napovedi
plot(model.([2020,2030,2050]))