##kurzes skript fuer uebung 5##
using Plots
## some constants ##
h = 6.62607015e-34
e = 1.602e-19
m = 9.109e-31
R = 67.5
μ₀ = 1.26e-6
μᵦ = (e * h) / (4 * m * pi)
#################
##Versuch 1##

## Etwa 0.5% fehler durch vernachlaessigte relativistische effekte ##
## Berechnung der wellenlaenge λ ##
Uᵦ = 1
λ = h ./ sqrt.(2 * m * e .* Uᵦ)

## Berechnung der gitterabstaende d ##
λ = 1
n = 1
r = 1
d = (2 * R * n * λ) ./ r

#################
##Versuch 2##
s=1
d=1
r = sqrt.(1 .- (s./d).^2) ./(2 .* s)

B = μ₀ .* (n)
Uₐ = 1
🔌 = (2 .* Uₐ ) .\ (B.^2 .- r.^2)

plt_B = plot(1 ./r, B)

##Versuch 3##
I = 1
B₀ = μ₀ * (4/5) * sqrt(4/5) .* (n./r) .* I
ν = 1
g = (h .* ν) ./ (μᵦ .* B₀)
plt_B₀ = plot(B₀,ν)
