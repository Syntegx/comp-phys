using Pkg
# Pkg.add CSV
# Pkg.add DataFrames
using CSV, DataFrames, LinearAlgebra, FFTW, TimerOutputs, Plots, WAV, Peaks
Data = CSV.read("C:\\Git\\Julia\\comp-phys\\Julia\\Data\\single_tone.txt", DataFrame, header=false, ignorerepeated=true, delim=' ')
## DISCRETE FOURIER TRANSFORM

function DFT!(Data)
    n = length(Data)
    k = 0:n-1
    j = zeros(ComplexF64, n, n)
    j[:, :] .= k'
    factor = im * 2 * pi * k ./ n
    j .= exp.(factor .* j)
    j[:, 1] .= LinearAlgebra.BLAS.gemv('C', j, Complex.(Data))
end

function checkFT(Data)
    isapprox(DFT!(Data), fft(Data))
end
# checkFT(Data)
##BENCHMARK COMPUTATION TIME FOR EVERY M in 1:100:10000 - takes approx 2 minutes, alloctes approx 51gb##
function timebench(Data)
    m = 1:100:10000
    b = []
    for i in 1:length(m)
        push!(b, (@timed DFT!(Data[1:m[i]])).time)
    end
    b
    plt = plot(m, b,
        title="Computation time as function of segment length",
        xlabel="Segment length m",
        ylabel="Computation time t in s")
    display(plt)
end

##PLOT COMPUTATION TIME OVER M##


##POWER SPECTRAL DENSITY##
function PSD(Data)
    n = length(Data)
    temp = fft(Data)
    temp .= (abs.(temp) .^ 2) ./ n^2
end

function PlotPSD(Data)
    TF = Real.(PSD(Data))
    idx = findmax(TF)[2] * 44100 / length(Data)
    nu = (0:length(Data)-1) * 44100 / length(Data)
    plt = plot(nu, Real.(PSD(Data)),
        xlim=[1, 44100 / 2],
        yscale=:log10,
        xscale=:log10,
        title="Frequency Analysis",
        xlabel="Frequency ν / Hz",
        ylabel="Power Spectral Density",
        label="PSD")
    vline!([idx, idx],
        label="Fundamental Frequency")
    plt2 = plot(nu, Real.(PSD(Data)),
        xlim=[1, 44100 / 2],
        yscale=:log10,
        xscale=:log10,
        title="Frequency Analysis",
        xlabel="Frequency ν / Hz",
        ylabel="Power Spectral Density",
        label="PSD")
    display(plt)

    println("The fundamental frequency is ", idx)
end

PlotPSD(Data[:, 1])

 function InverseFFTUtilized(Data)
    n = length(Data)
    temp = real.(PSD(Data))
    freqz, vals = findmaxima(temp, 100)
    freqz = freqz[1:6]
    m = length(freqz)
    temp_2 = zeros(n)
    temp = fft(Data)
    t = 0:1/44100:(n-1)/44100
    factor = t * freqz' .* 2 * pi *44100 / n

    for i in 1:m
        temp_2 .+= real.(temp[freqz[i]] / n) * cos.(factor[:, i]) + imag.(temp[freqz[i]] / n) * sin.(factor[:, i])
    end
 
    temp_2
 end

audio_1 = InverseFFTUtilized(Data[:, 1])
audio_2 = InverseFFTUtilized(Data[:, 2])

wavwrite([audio_1 audio_2], "test.wav", Fs=44100)
