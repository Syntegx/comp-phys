using DataFrames, CSV, StatsBase, Distributions, CairoMakie, Statistics
df=CSV.read("data.txt",DataFrame; header = false)
df=df[:,1];
x = (LinRange(findmin(df)[1]-1,findmax(df)[1]+1,100))[:]
mu=mean(df)
sig=std(df)

MyNormalPdf = pdf.(Normal(mu,sig),x)
MyCauchyPdf = pdf.(Cauchy(mu,sig),x)
f=Figure()
ax = Axis(f[1,1], limits=(findmin(df)[1]-1,findmax(df)[1]+1 ,nothing,nothing));
hist!(ax,df,x;color="lightblue", normalization = :pdf);
lines!(ax,x,MyNormalPdf;color="red", linestyle= :dot);
lines!(ax,x,MyCauchyPdf;color="green");
h = fit(Histogram,df);

sum(MyNormalPdf);
sum(MyCauchyPdf);
BFNorm = prod(MyNormalPdf);
BFCauchy = prod(MyCauchyPdf);
Odds = BFNorm/BFCauchy
f