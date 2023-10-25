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
function CheckDiagDominance(A)
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
function MakeDiagDominant!(A)
    bool = true
    k = 0
        x = CheckDiagDominance(A)
        if any(x .> 0)
            k += maximum(x)
            if typeof(A[1]) == Complex{Float64}
                A[diagind(A)] .+= (k * im)
            else
                A[diagind(A)] .+= (k)
            end
        end
    return k * im
end

# defining a function with two methods, one for solving a complex system
function MyGS(A::Matrix{ComplexF64}, b::Vector{ComplexF64}, maxitr=30)
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
function MyGS(A::Matrix{Float64}, b::Vector{Float64}, maxitr=30)
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
    MakeDiagDominant!(A, b)
    n = length(b)
    temp = [deepcopy(b), zeros(n)]
    k = 0
    storage = []
    plt = plot(1,
        legend=false,
        xlabel=("Iteration p"),
        ylabel=("Absolute difference of sum"),
        title="Step Difference using MyGaussSeidel")
    while k < maxitr
        temp[2] .= MyGS(A, temp[1], 1)
        # plot!(1:n, abs.(temp[2] .- temp[1]))
        push!(storage, sum(abs.(temp[2] .- temp[1])))
        temp[1] .= temp[2]
        k += 1
    end
    plot!(1:maxitr, storage)
    display(plt)

end
# a wrapper function to call TestMyGS function for different random matrices
function Task_c(first_n::Int64, second_n::Int64, maxitr::Int64)
    TestMyGS(rand(first_n, first_n), rand(first_n), maxitr)
    TestMyGS(rand(second_n, second_n), rand(second_n), maxitr)
end


function IterativeGreenFunction(A, b, maxitr=30)
    scale_factor = 5*im
    A[diagind(A)].+=scale_factor
    k = 0
    temp = deepcopy(b)
    G_p = zeros(ComplexF64,length(b))
    while k < maxitr
        #   G_p .= A\temp
        G_p .= MyGS(A, temp)
        temp .= b .+ (scale_factor .* G_p)
        k += 1
    end
    G_p
end


function Task_d_part1(N, t, E, delta, alpha, beta)
    n = length(E)
    storage = []
    for j in 1:n
        A = ring(N, t, E[j], delta, alpha, beta)
        G_R = zeros(Complex{Float64}, N, N)
        B = deepcopy(G_R)
        B[diagind(B)] .= 1
        # if any(CheckDiagDominance(A) .> 0)
        for i in 1:N
            # G_R[:, i] .= A\B[:,i]  
            # G_R[:, i] .= MyGS(A,B[:,i])
            G_R[:, i] .= IterativeGreenFunction(A, B[:, i])
        end
        # else
        #     for i in 1:N
        #         G_R[:, i] .= MyGS(A, B[:, i])
        #     end
        # end
        push!(storage, G_R)
    end
    storage
end

function Task_d_part2(x, alpha, beta)
    storage = []
    for i in 1:length(x)
        push!(storage, [x[i][alpha, alpha], x[i][beta, beta], x[i][alpha, beta], x[i][beta, alpha]])
    end
    storage
end

function Task_d(N, t, E, delta, alpha, beta)

    storage = Task_d_part1(N, t, E, delta, alpha, beta)
    storage = Task_d_part2(storage, alpha, beta)
    storage
end

function Task_e(N, t, E, delta, alpha, beta)
Hamiltonian = complex(zeros(N, N))
Hamiltonian .= SymTridiagonal(fill(0.0, N), fill(t, N))
Hamiltonian[CartesianIndex.([1, N], [N, 1])] .= t
eigenvalues =eigen(Hamiltonian).values

storage = []
plt = plot(
        legend = false,
        title = "Transmission probability",
        xlabel = "Energy",
        ylabel = "Probability"
            )
for i in 1:length(E)
   push!(storage, Task_d(N, t, E, delta, alpha[1], beta[1])[i][3])
end
plot!(E, abs2.(storage))
vline!([eigenvalues])

storage = []
for i in 1:length(E) 
    push!(storage, Task_d(N, t, E, delta, alpha[2], beta[2])[i][3])
end
plot!(E, abs2.(storage))
display(plt)
eigenvalues
end



Task_e(6, -2.6, LinRange(-6,6,1000), 0.5, [1,1], [3,4])
# A =ring(6, -2.6,0, 0.5, 1, 3);
# b=zeros(ComplexF64,6);
# b[5]=1;
# N=6;
# t=-2.6;
# E=LinRange(-6,6,10)
# delta=0.5
# alpha = 1
# beta = 3