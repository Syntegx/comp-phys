using HypothesisTests, Distributions, CairoMakie, StatsBase
global alpha = Int8.([1, 2, 3]);
global N = Int16.([5, 10, 100, 1000]);
global q = Float16.([0.3, 0.7]);
global j = 1;
global k = 1;
global rgb = ["red", "blue", "green"];
f = Figure()
g = Figure()
for N_i in N

    ax = Axis(f[j, 1], title="Probability at N = $N_i", ylabel = "P(Hα|𝑛,𝑁,𝐵)", xlabel = "𝑛" )
    axg = Axis(g[j, 1], limits = (0,nothing,0,20), title="Odds at N = $N_i",ylabel = "Odds", xlabel = "𝑛")
    MyBinomial1 = Binomial(N_i, q[1])
    MyBinomial2 = Binomial(N_i, q[2])


    MyPdf1 = pdf.(MyBinomial1, 1:1:N_i) ./ 3
    MyPdf2 = pdf.(MyBinomial2, 1:1:N_i) ./ 3
    MyPdf3 = (1 / (1 + N_i)) * ones(N_i) / 3

    lines!(ax, 1:1:N_i, MyPdf1; color=rgb[1], label= "Urne α=1")
    lines!(ax, 1:1:N_i, MyPdf2; color=rgb[2], label= "Urne α=2")
    lines!(ax, 1:1:N_i, MyPdf3; color=rgb[3], label= "Urne α=3")

    lines!(axg, 1:1:N_i, MyPdf1 ./ MyPdf2; color=rgb[1], label= "Urne α=1")
    lines!(axg, 1:1:N_i,  MyPdf1 ./ MyPdf3; color=rgb[2], label= "Urne α=2")
    lines!(axg, 1:1:N_i, MyPdf2 ./ MyPdf3; color=rgb[3], label= "Urne α=3")
    

    Legend(f,ax, "Urnen", orientation = :horizontal ,valign=true)
    
    global j += 1
end

g
f