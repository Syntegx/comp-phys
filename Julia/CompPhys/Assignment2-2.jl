### need these packages in current environment ###
using LinearAlgebra, CSV, Plots, DataFrames, Measurements
### reading Data ###
### input own path for data ### 
Data = CSV.read("C:\\Git\\Julia\\comp-phys\\Julia\\Data\\time_evolution.txt", DataFrame, header=false, delim=',')
## Task (b) - '\' operator from LinearAlgebra.jl uses LU Decomposition, Forward-\Backwardsubstitution    
Datax = Data[:, 1]
Datay = Data[:, 2]
### solver###
function Solve(A, b)
    y = A \ b
end
### Linear fitting routine using polynomial basis functions ### 
function LinearParameterFit(Datax, Datay, m=4)
    weight = 1 ./ Datay

    n = length(Datax)

    phi = 0:m-1

    temp_dat = Datax .^ phi'

    beta = sum(temp_dat .* Datay .* weight, dims=1)[:]

    a = ones(m, m) .* phi
    a .+= phi'
    temp = sum(weight .* Datax .^ a[:]', dims=1)

    temp = reshape(temp, (m, m))

    fit_params = Solve(temp, beta)
    fit = sum((Datax .^ phi') .* fit_params', dims=2)
end

### Wrapper function calling LinearParameterFit function ###
### Plotting fit up to input term m ###
function Task_d(Datax, Datay, m)
    plt = scatter(Datax, Datay,
        markersize=1,
        label="Data",
        title="Linear Fit using polynomial basis-functions")
    for i in 2:m
        fit = LinearParameterFit(Datax, Datay, i)
        plt = plot!(Datax, fit,
            label="m = $i")
    end
    display(plt)
end

Task_d(Datax, Datay, 7)


### Fitting routine specifically implemented for underlying problem of exercise 2###
function NonLinearParameterFit(Datax, Datay, m=3)
    k = (1:m-1)
    n = length(Datax)
    weight = 1 ./ Datay
    phi = [ones(n) 1 .- exp.(-Datax ./ k')]

    beta = sum(Datay .* weight .* phi, dims=1)'

    A = zeros(m, m)

    for i in 1:m
        for j in 1:m
            A[i, j] = sum(weight .* phi[:, i] .* phi[:, j])
        end
    end

    fit_parameter = Solve(A, beta)
    fit = sum(phi .* fit_parameter', dims=2)
    σ = sqrt.(inv(A)[diagind(A)])
    c = fit_parameter[1] - 1
    alpha = sum(-(1 / c) .* fit_parameter[2:m] ./ k)
    fit2 = 1 .+ c .* exp.(-alpha .* Datax)
    fit, fit2, alpha, c, alpha, σ, fit_parameter
end
### Wrapper function used for plotting nonlinear fit of Data using NonLinearParameterFit function ###
### plotting fits up to input term m and calculating standard deviation through error propagation ###
function Task_eh(Datax, Datay, m)
    plt = scatter(Datax, Datay,
        markersize=1,
        label="Data",
        title="Nonlinear Fit using e-basis functions")
    plt2 = scatter(Datax, Datay,
        markersize=1,
        label="Data",
        title="Nonlinear Fit using e-function and parameters")
    for i in 1:m
        fit = NonLinearParameterFit(Datax, 1 ./ Datay, i)
        plot!(plt, Datax, 1 ./ fit[1],
            label=("m = $i"))


        plot!(plt2, Datax, 1 ./ fit[2],
            label=("m = $i"))
    end
    fit = NonLinearParameterFit(Datax, 1 ./ Datay, m)
    sigma = measurement.(fit[7][:], fit[6][:])

    c = sigma[1] - 1
    α = sum(-(1 / c) .* sigma[2:m] ./ (1:m-1))

    plt3 = scatter(Datax, Datay,
        markersize=1,
        label="Data")
    plot!(plt3, Datax, 1 ./ fit[2],
        title="α = $α, c=$c",
        label="(1 + c exp(-α⋅t))⁻¹")

    display(plt)
    display(plt2)
    display(plt3)
end

Task_eh(Datax, Datay, 6)