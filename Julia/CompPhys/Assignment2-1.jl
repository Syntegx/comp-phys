using CSV, DataFrames, FFTW, BenchmarkTools
#FOURIER TRANSFORM#

df = CSV.read("E://Julia//comp-phys-1//Julia//Data//single_tone.txt", delim=',', decimal='.', DataFrame, header=false)[1:1000, 2:3];

 function FT(df)
n = length(df[:, 1])
Y = zeros(ComplexF64, n)
M = im * 2 * pi .* collect(0:n-1) ./ n
Y .= sum.(df[:, 1] .* exp.(M * j) for j in 0:n-1)
T = fft(df[:,1])
Bool = all(isapprox.(abs.(Y),abs.(T)))
return Y, Bool
end

function Task_b()
    
end