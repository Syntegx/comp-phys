using LinearAlgebra, IterativeSolvers, Plots
### DEFINING FUNCTION RING USED TO INITIALIZE USED MATRICES ###
# N specifies the used dimension 
# t specifies the probability for an electron to move from state b to state a
# E is the energy of the system
# delta is the value of the energy broadening of entering/leaving the system
# alpha/beta are the sites of entering/leaving the system
function ring(N, t, E, delta, alpha, beta)
    Matrix = complex(zeros(N, N))
    Matrix .= SymTridiagonal(fill(-Float64(E), N), fill(t, N))
    Matrix[CartesianIndex.([1, N], [N, 1])] .= t
    Matrix[CartesianIndex.([alpha, beta], [alpha, beta])] .+= delta * 1im
    Matrix
end

# for the Gauss Seidel Method beeing able to converge, the matrix A has to be diagonally dominant
# this function checks for diagonal dominance
function CheckDiagDominance(A::Matrix{Any})
    n = size(A, 1)
    x = zeros(n)
    for i in 1:n
        x[i] = sum(abs.(A[i, :])) .- (2 * abs.(A[i, i]))
    end
    # if all(x .< 0)
    #     println("Matrix is diagonal dominant")
    # end
    x
end

# this function alters the matrix A and vector b by adding a term k to the diagonal elements 
# this makes sure the Gauss Seidel method converges
function MakeDiagDominant(A::Matrix{Any}, b::Vector{Any})
    bool = true
    iter = 0
    while bool
        x = CheckDiagDominance(A)
        if iter > 100000
            bool = false
        elseif any(x .> 0)
            k = maximum(x)
            A[diagind(A)] .+= k
            b .+= k
        else
            bool = false
        end
        iter += 1
    end
    A
end

# defining a function with two methods, one for solving a complex system
function MyGS(A::Matrix{ComplexF64}, b::Vector{ComplexF64}, maxitr=10)
    n = length(b)
    x = zeros(Complex{Float64}, n)
    k = 0
    while k < maxitr
        for col in 1:n
            @simd for row = 1:col-1
                @inbounds x[row] -= A[row, col] * x[col]
            end
            x[col] = b[col]
        end

        for col = 1:n
            @inbounds x[col] /= A[col, col]
            @simd for row = col+1:n
                @inbounds x[row] -= A[row, col] * x[col]
            end
        end
        k += 1
    end
    return x
end

# one for solving a real system
function MyGS(A::Matrix{Float64}, b::Vector{Float64}, maxitr=10)
    n = length(b)
    x = zeros(Float64, n)
    k = 0
    while k < maxitr
        for col in 1:n
            @simd for row = 1:col-1
                @inbounds x[row] -= A[row, col] * x[col]
            end
            x[col] = b[col]
        end

        for col = 1:n
            @inbounds x[col] /= A[col, col]
            @simd for row = col+1:n
                @inbounds x[row] -= A[row, col] * x[col]
            end
        end
        k += 1
    end
    return x
end


# a function for plotting all steps between iterations of Gauss Seidel method
function TestMyGS(A, b, maxitr=10)
    MakeDiagDominant(A, b)
    n = length(b)
    temp = [deepcopy(b), zeros(n)]
    k = 0
    plt = plot(1,
        legend=false,
        xlabel=("Iteration p"),
        ylabel=("Absolute Difference"),
        title="Step Difference using MyGaussSeidel")
    while k < maxitr
        temp[2] .= MyGS(A, temp[1], 1)
        plot!(1:n, abs.(temp[2] .- temp[1]))
        temp[1] .= temp[2]
        k += 1
    end
    display(plt)
end
# a wrapper function to call TestMyGS function for different random matrices
function Task_c(first_n::Int64, second_n::Int64, maxitr::Int64)
    TestMyGS(rand(first_n, first_n), rand(first_n), maxitr)
    TestMyGS(rand(second_n, second_n), rand(second_n), maxitr)
end

function Task_d(N, t, E, delta, alpha, beta)
    A = ring(N, t, E, delta, alpha, beta)
    G_R = zeros(Complex{Float64}, N, N)
    B = deepcopy(G_R)
    B[diagind(B)] .= 1
    for i in 1:N
        G_R[:, i] .= IterativeGreenFunction(A, B[:, i])
    end
    G_R
end

function IterativeGreenFunction(A, b, maxitr=10)
    
    k = 0
    while k < maxitr
        G_p = MyGS(A, b)
        b = G_p
        k += 1
    end
    b
end

A = ring(6, -2.6, -6, 0.5, 1, 3);
b = zeros(Complex{Float64}, 6);
b[1] = 1
