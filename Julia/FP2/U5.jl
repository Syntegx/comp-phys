##kurzes skript fuer uebung 5##
using Plots, CSV, DataFrames, CurveFit, Measurements
## some constants ##
h = 6.6256e-34
e = 1.602e-19
m = 9.109e-31
R = 67.5 * 1e-3
μ₀ = 1.26e-6
μᵦ = (e * h) / (4 * m * pi)
#################
##Versuch 1##
## Etwa 0.5% fehler durch vernachlaessigte relativistische effekte ##
## Berechnung der wellenlaenge λ ##
function U5_1()
df = CSV.read("C:\\Git\\Julia\\comp-phys\\Julia\\Data\\U5_1.txt", DataFrame, header=1, delim=',')
Uᵦ = measurement.(df[:,3],df[:,3].*0.01) .* 1e3
r = measurement.(df[:,1:2],df[:,1:2].*0.06) .* 0.5e-2
λ = h ./ sqrt.(2 * m * e .* Uᵦ)
## Berechnung der gitterabstaende d ##

linfit = linear_fit(λ,r[:,1])
linfit2 = linear_fit(λ,r[:,2])

d1 = 1/(linfit[2] / (2 * R))
d2 = 1/(linfit2[2] / (4 * R))

λ.*=1e12
r.*=1e2
plt = scatter(λ,r[:,1],
label="Beugungsordnung: 1",
ylabel="r / cm",
xlabel="λ / pm",
title = "Radius in Abhängigkeit von der Wellenlänge")
scatter!(λ,r[:,2],
label="Beugungsordnung: 2")


linfit = linear_fit(λ,r[:,1])
linfit2 = linear_fit(λ,r[:,2])


plot!(Measurements.value.(λ),Measurements.value.(linfit[1].+linfit[2] .* λ),
label="Linearer Fit Beugungsordnung: 1")
plot!(Measurements.value.(λ),Measurements.value.(linfit2[1].+linfit2[2] .* λ),
label="Linearer Fit Beugungsordnung: 2")
display(plt)
savefig(plt,"Elektronenbeugung.png")
println("Gitterabstand d1 beträgt: ", d1,"  ","Gitterabstand d2 beträgt: ", d2)

end
U5_1()

#################
##Versuch 2##
function U5_2()
df = CSV.read("C:\\Git\\Julia\\comp-phys\\Julia\\Data\\U5_2.txt", DataFrame, header=1, delim=',')
d = 135 *1e-3
s = measurement.(df[:,1],2) *1e-3
r = d^2 .* sqrt.(1 .- (s./d).^2) ./(2 .* s)

n = 320
I = measurement.(df[:,2],df[:,2].*0.01 .+ 0.003)

H = 33.8*1e2 .*I
B = μ₀ .* (H)
Uₐ = measurement.(df[:,3], df[:,3].*0.005) .*1e3

linfit = linear_fit(B[1:4],1 ./ r[1:4])
linfit2 = linear_fit(B[5:end],1 ./r[5:end])

d2=-linfit2[2]^2 .* 2 * 2500
d1=-linfit[2]^2 .* 2 *4500
-e/m

B .*= 1e3


plt = scatter(B[1:4], 1 ./r[1:4],
label="Beschleunigungsspannung U₁",
ylabel="r⁻¹ / cm⁻¹",
xlabel="B / mT",
title = "Ablenkung in Abhängigkeit vom Magnetfeld")
scatter!(B[5:end], 1 ./r[5:end],
label="Beschleunigungsspannung U₂")

linfit = linear_fit(B[1:4],1 ./ r[1:4])
linfit2 = linear_fit(B[5:end],1 ./r[5:end])


plot!(Measurements.value.(B),Measurements.value.(linfit[1].+linfit[2] .* B),
label="Linearer Fit für U₁")
plot!(Measurements.value.(B),Measurements.value.(linfit2[1].+linfit2[2] .* B),
label="Linearer Fit für U₂")

display(plt)

savefig(plt,"Elektronenablenkung.png")
println("Messergebnisse e/m : ", d1," und ",d2)

end
 U5_2()

##Versuch 3##
function U5_3()
df = CSV.read("C:\\Git\\Julia\\comp-phys\\Julia\\Data\\U5_3.csv", DataFrame, header=1, delim=',')
I = measurement.(df[:,2],((df[:,2].*0.06 .+ 0.003)))
ν = measurement.((df[:,1]),(df[:,1].* 0.01)).* 1e6 

n = 320
r = 68 * 1e-3

B₀ = μ₀ * (4/5) * sqrt(4/5) .* (n./r) .* I

fit = linear_fit(B₀.*1e3, ν*1e-6)
plt = scatter(B₀.*1e3, ν*1e-6,
label="Daten",
ylabel="ν / MHz",
xlabel="B₀ / mT",
title = "Resonanzfrequenz in Abhängigkeit vom Magnetfeld")
plot!(Measurements.value.(B₀.*1e3),Measurements.value.(fit[1].+fit[2].*B₀*1e3),
label="Linearer Fit")


fit = linear_fit(B₀, ν)
g = fit[2] * (h/μᵦ)

display(plt)

savefig(plt,"ESR.png")
println("g = ", g)
end
U5_3()