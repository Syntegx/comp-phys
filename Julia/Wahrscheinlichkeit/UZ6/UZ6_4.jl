using Random, Distributions, CairoMakie, StatsBase, Statistics, Weave
#Intialisieren der Variablen und insbesonderen der Speichertypen

N = Int64.([10 100 1000]);
M = Int64.([100 1000 10000]);

# Init der Figures f und g
f = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98),
    resolution=(1000, 700));
g = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98),
    resolution=(1000, 700));
#Init Laufvariable für Indizierung der Subplots
global i = k = 1;
#Random walk - gleiche Wahrsch. um in beide Richtungen zu gehen
#Iterationen über M
for M_i in M
    #Iterationen über N
    for N_i in N
        #Condition fuer subplotting Variablen
        if i == 4
            global i = 1
            global k += 1
        end
        #Definition der Subplotfenster und Label
        ax = Axis(f[i, k], title="N = $N_i", ylabel="M = $M_i")
        ax2 = Axis(g[i, k], title="N = $N_i", ylabel="M = $M_i")
        #N x M faches Sampling aus Pop -1 +1
        global x = sum(sample([-1, 1], (N_i, M_i)), dims=1)
        #Fitting der Normal Distribution
        norm = Normal(mean(x), std(x))
        #Dictionary der Unique Werte von X und deren Vielfachheit
        uniquex = keys(sort(countmap(x)))
        #Normiertes Histrogram vom Sampling x
        hist!(ax, x[:]; bins=round(Int64, length(unique(x))), normalization=:pdf)
        #Fit der Normal Distribution - etwas off
        lines!(ax, LinRange(findmin(x)[1], findmax(x)[1], 100), norm; color="red")
        #Empirical cumulative distribution function wird gefittet
        local c = ecdf(x[:])
        #Plotten der ECDF c
        barplot!(ax2, LinRange(findmin(x)[1], findmax(x)[1], length(x[:])), c)
        #Plotten der CDF von der Normal distribution
        lines!(ax2, LinRange(findmin(x)[1], findmax(x)[1], 100), cdf(norm, LinRange(findmin(x)[1], findmax(x)[1], 100)); color="red")

        global i += 1
    end
end
f
g

