using LinearAlgebra, CSV, Plots, DataFrames

Data = CSV.read("C:\\Git\\Julia\\comp-phys\\Julia\\Data\\time_evolution.txt", DataFrame, header=false, delim=',')
## Task (b) - '\' operator from LinearAlgebra.jl uses LU Decomposition, Forward-\Backwardsubstitution    
Datax = Data[:, 1]
Datay = Data[:, 2]
function Solve(A, b)
    y = A \ b
end

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
function Task_d(Datax, Datay, m)
    plt = plot(Datax, Datay)
    for i in 2:m
        fit = LinearParameterFit(Datax, Datay, i)
        plt = plot!(Datax, fit)
    end
    display(plt)
end

# Task_d(Datax, Datay, 10)

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
    c = fit_parameter[1] - 1
    alpha = sum(-(1 / c) .* fit_parameter[2:m] ./ k)
    fit2 = 1 .+ c .* exp.(-alpha .* Datax)
    fit, fit2, alpha, c
end


fit =  NonLinearParameterFit(Datax, 1 ./Datay, 5)

 plot(Datax,Datay)
 plot!(Datax,1 ./ fit[1])
 plot!(Datax,1 ./ fit[2])

 function Task_e(Datax, Datay, m)
    plt = plot(Datax, Datay)
    for i in 1:m
        fit = NonLinearParameterFit(Datax, Datay, i)
        plt = plot!(Datax, fit)
    end
    display(plt)
end

# Task_e(Datax,Datay,1)
