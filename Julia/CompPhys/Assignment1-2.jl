using LinearAlgebra, Plots
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
# outputs a vector x containing the difference of diagonal elements to the sum of the row elements
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

# this function alters the matrix A by adding a term k to the diagonal elements (mutating - alters the input matrix A)
# uses the CheckDiagDominance function to determine the greatest difference of diag elements to sum of elements
# this makes sure the Gauss Seidel method converges
function MakeDiagDominant!(A)
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

# defining a function utilizing gauss seidel algorithm with two methods, one for solving a complex system
# julias complex.jl package manages all complex number operations
function GaussSeidel(A::Matrix{ComplexF64}, b::Vector{ComplexF64}, maxitr=30)
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
function GaussSeidel(A::Matrix{Float64}, b::Vector{Float64}, maxitr=30)
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


# a wrapper function for plotting all steps between iterations of Gauss Seidel method
# tests the defined function GaussSeidel with input matrix A and vector b, maximum of iterations is optional input
function TestGaussSeidel(A, b, maxitr=10)
    MakeDiagDominant!(A)

    n = length(b)
    temp = [deepcopy(b), zeros(n)]
    k = 0
    storage = []
    while k < maxitr
        temp[2] .= GaussSeidel(A, temp[1], 1)
        # plot!(1:n, abs.(temp[2] .- temp[1]))
        push!(storage, sum(abs.(temp[2] .- temp[1])))
        temp[1] .= temp[2]
        k += 1
    end
    storage
end
# a wrapper function to call TestGaussSeidel function for different random matrices of dimension n
# the inputs first_n and second_n dictate the used dimensions
function Task_c(first_n=100::Int64, second_n=1000::Int64, maxitr=10::Int64)

    storage = TestGaussSeidel(rand(first_n, first_n), rand(first_n), maxitr)
    plt = plot(1:maxitr, storage,
        xlabel=("Iteration p"),
        ylabel=("Absolute difference of sum"),
        title="Step Difference using GaussSeidel")
    png(plt, "Matrix100")
    display(plt)
    
    storage = TestGaussSeidel(rand(second_n, second_n), rand(second_n), maxitr)
    plt = plot(1:maxitr, storage,
        xlabel=("Iteration p"),
        ylabel=("Absolute difference of sum"),
        title="Step Difference using GaussSeidel")
    display(plt)
    png(plt, "Matrix1000")
end

# a function implementing the iterative solution of greens functions using hint 3 from the assignment
function IterativeGreenFunction(A, b, maxitr=30)
<<<<<<< HEAD
    scale_factor = 0
    A[diagind(A)].+=scale_factor
=======
    #scale_factor is constant imaginary as suggested by tutors
    scale_factor = 5 * im
    A[diagind(A)] .+= scale_factor
>>>>>>> 1c2cfdf977ec921c5509003bb0738ab39438ee2e
    k = 0
    temp = deepcopy(b)
    G_p = zeros(ComplexF64, length(b))
    while k < maxitr
<<<<<<< HEAD
        G_p .= A\temp
        #  G_p .= MyGS(A, temp)
=======
        G_p .= GaussSeidel(A, temp)
>>>>>>> 1c2cfdf977ec921c5509003bb0738ab39438ee2e
        temp .= b .+ (scale_factor .* G_p)
        k += 1
    end
    G_p
end

# a function computing a part of the problem in Task_d
# calculates the greens function solving the linear system of different matrices with different energies
# for now it uses the inbuild solver of julia since my implementation of the iterative solution
# for greens function wont work as intended

function Task_d_part1(N, t, E, delta, alpha, beta)
    n = length(E)
    storage = []
    for j in 1:n
        A = ring(N, t, E[j], delta, alpha, beta)
        G_R = zeros(Complex{Float64}, N, N)
        B = deepcopy(G_R)
        B[diagind(B)] .= 1
        for i in 1:N
            G_R[:, i] .= A \ B[:, i]
            # G_R[:, i] .= GaussSeidel(A,B[:,i])
            # G_R[:, i] .= IterativeGreenFunction(A, B[:, i])
        end
        push!(storage, G_R)
    end
    storage
end

# a function pushing all GR_aa, GR_bb, GR_ab and GR_ba into a storage array
function Task_d_part2(x, alpha, beta)
    storage = []
    for i in 1:length(x)
        push!(storage, [x[i][alpha, alpha], x[i][beta, beta], x[i][alpha, beta], x[i][beta, alpha]])
    end
    storage
end

# a wrapper function utilizing part 1 and part 2 and part 3of the Task_d functions to complete the whole Task_d
# stores all calculated GR_aa, GR_bb, GR_ab and GR_ba in a storage array
function Task_d(N=6, t=-2.6, E=LinRange(-6, 6, 100), delta=0.5, alpha=[1, 1], beta=[3, 4])
    storage = []
    for i in 1:length(alpha)
        storage = Task_d_part1(N, t, E, delta, alpha[i], beta[i])
        storage = Task_d_part2(storage, alpha[i], beta[i])
    end
    storage
end

# a function thats calculating all eigenvalues of the Hamiltonian H_R using an inbuild eigen solver 
# plotting all GR_ab of both systems used in (c) as probability of transmission for different energy levels Energy
function Task_e(N=6, t=-2.6, E=LinRange(-6, 6, 100), delta=0.5, alpha=[1, 1], beta=[3, 4])
    Hamiltonian = complex(zeros(N, N))
    Hamiltonian .= SymTridiagonal(fill(0.0, N), fill(t, N))
    Hamiltonian[CartesianIndex.([1, N], [N, 1])] .= t
    eigenvalues = eigen(Hamiltonian).values

    storage = []

    plt = plot(
        legend=:top,
        title="Transmission probability",
        xlabel="Energy",
        ylabel="Probability"
    )


    for i in 1:length(E)
        push!(storage, Task_d_part3(N, t, E, delta, alpha[1], beta[1])[i][3])
    end
    plot!(E, abs2.(storage),
        label="System 1"
    )

    storage = []
    for i in 1:length(E)
        push!(storage, Task_d_part3(N, t, E, delta, alpha[2], beta[2])[i][3])
    end
    plot!(E, abs2.(storage),
        label="System 2"
    )

    vline!([eigenvalues],
        label="Eigenvalues")
    png(plt, "Transmission")
    display(plt)
end

### RUNNING FUNCTIONS Task_c AND Task_d AND Task_e ###

<<<<<<< HEAD

Task_e(6, -2.6, LinRange(-6,6,100), 0.5, [1,1], [3,4])
# A =ring(6, -2.6,0, 0.5, 1, 3);
# b=zeros(ComplexF64,6);
# b[5]=1;
# N=6;
# t=-2.6;
# E=LinRange(-6,6,10)
# delta=0.5
# alpha = 1
# beta = 3
=======
Task_c()
Task_d(6, -2.6, LinRange(-6, 6, 100), 0.5, 1, 3);
Task_d(6, -2.6, LinRange(-6, 6, 100), 0.5, 1, 4);
Task_e()
>>>>>>> 1c2cfdf977ec921c5509003bb0738ab39438ee2e
