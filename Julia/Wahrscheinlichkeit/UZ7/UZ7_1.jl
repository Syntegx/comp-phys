using HypothesisTests, CSV, DataFrames, Distributions
df=CSV.read("Lotto.csv",DataFrame; header = false)
df=df[:,2];
alpha = 0.01;

prob = 1/length(df)
x0n = prob * sum(df)


chi=sum((df[:].-x0n) .^ 2 ./ x0n)
x0 = quantile(Chisq(length(df)-1),1-alpha)