using CSV, StatsBase, Distributions, CairoMakie, Statistics, DelimitedFiles, Optim, LinearAlgebra
df = readdlm("mcmc.dat");
x = df[1, :];
y = df[2, :];

#Plot
f = Figure()
ax = Axis(f[1, 1]; xlabel=("Steuergröße x"), ylabel=("Messgröße y"))
#Plotten der Daten
errorbars!(ax, x, y, 0.75, color=:red)
scatter!(ax, x, y, markersize=5, color=:black, label="Data")
#Define function phi, minimizing phi
phi(var) = (mean(y .^ 2) + var[1] .^ 2 + var[2] .^ 2 * mean(x .^ 2) - 2var[1] * mean(y) - 2var[2] * mean(x .* y) + 2var[1] * var[2] * mean(x));
minimizingVars = Optim.minimizer(optimize(phi, [0.0, 0.0]));
#Plot linfit
lines!(ax, x, minimizingVars[1] .+ minimizingVars[2] * x, color= :yellow, linestyle=:dashdot, label="ML");

#Markov Chain
a = [0.0, 0.0]';            #... Startwerte
da = 0.34;           #... maximale Schrittweite
N_mess = 10000;          #... Messungen
N_skip = 20;            #... ausgelassene Werte zwischen den Messungen
N_loop = N_mess * N_skip; #... Schritte insgesamt

N = length(x);
alpha = N / (2 * 0.75^2);

phi_0 = phi(a);
ico_acc = 0;
list = zeros(N_mess, 2); 
list_corr = zeros(N_loop, 2); #... korrelierte Messdaten

for i in 1:N_loop
    a_t = a .+ (rand(1, 2) .- 0.5) .* da
    phi_t = phi(a_t)

    if log.(rand(1))[1] < alpha .* (phi_0 - phi_t)
        global ico_acc = ico_acc + 1
        global a = a_t
        global phi_0 = phi_t
    end
    list_corr[i, :] = a
    list[Int16(ceil(i / N_skip)), :] = a
end


avg_a = mean(list[Int16(ceil(0.1 * N_mess)):20:end, :], dims=1);
C = cov(list[Int16(ceil(0.1 * N_mess)):20:end, :]);
std = sqrt.(diag(C));

lines!(ax, x, avg_a[1] .+ avg_a[2] * x, color= (:blue, 0.5), label = "MCMC");
leg = Legend(f[1,2],ax)
f

