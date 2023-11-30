using CSV, DataFrames, LinearAlgebra, Distances, Plots
function read_data(Path="C:\\Git\\Julia\\comp-phys\\Julia\\Data\\xyzm_dna.txt")
    R_0 = Matrix(CSV.read("C:\\Git\\Julia\\comp-phys\\Julia\\Data\\xyzm_dna.txt", DataFrame, header=false, delim=','))
    n = length(R_0[:, 1])
    R_0, n
end

R_0 = read_data()
n = R_0[2]
R_0 = R_0[1]
function GaussianHesse(R_0, k=1)
    n = length(R_0[:, 1])
    R = zeros(n, n)
    for i in 1:n
        for j in 1:n
            R[i, j] = ((R_0[i, 1] - R_0[j, 1])^2 + (R_0[i, 2] - R_0[j, 2])^2 + (R_0[i, 3] - R_0[j, 3])^2) < 25
        end
    end
    R .*= (-1)
    R[diagind(R)] .= 0
    R[diagind(R)] .-= sum(R, dims=2).-1
    # R .*= k
    R
end
####
# sum(R_0,dims=2)
# R2 = GaussianHesse2(R_0)

R = GaussianHesse(R_0)


# (R[diagind(R)].==1.0)

function StiffnessMatrix(R_0, R)
    n = length(R_0[:, 1])
    M = zeros(n, n)
    M[diagind(M)] .= 1 ./ sqrt.(R_0[1:n, 4])
    M .= M * R * M
end

K = StiffnessMatrix(R_0, R)

function PowerMethod(K)
    n = length(K[:, 1])
    x = normalize(rand(n))
    x = rand(n)
    Δλ = 1
    λ = 0
    while Δλ > 1e-16
        λ = (x' * K * x) / (x' * x)
        x .= normalize(K * x)
        Δλ = abs(λ - (x' * K * x) / (x' * x))
    end
    λ, x
end


function MatrixDeflation!(K, x, λ)
    K .-= (λ * x * x')
    K
end

function PMMD(K, limit=length(K[:, 1]))
    C = deepcopy(K)
    k = 0
    λ = []
    x = []
    while k < limit
        temp = PowerMethod(C)
        push!(λ, temp[1])
        push!(x, temp[2])
        MatrixDeflation!(C, temp[2], temp[1])
        k += 1
    end
    λ, x
end

function Gershgorin(K)
    C = deepcopy(K)
    temp = maximum(sum(abs.(K), dims=2) + abs.(K[diagind(K)]))
    C ./= temp
    C
end

function Plot_DNA(R_0)
    n = length(R_0[:, 1])
    plt = plot(
        R_0[Int.(1:n/2), 1],
        R_0[Int.(1:n/2), 2],
        R_0[Int.(1:n/2), 3],
        linecolor=:blue
    )

    scatter!(plt,
        R_0[Int.(1:n/2), 1],
        R_0[Int.(1:n/2), 2],
        R_0[Int.(1:n/2), 3],
        markersize=0.5,
        markercolor=:blue)

    plot!(plt,
        R_0[Int.((n/2+1):end), 1],
        R_0[Int.((n/2+1):end), 2],
        R_0[Int.((n/2+1):end), 3],
        linecolor=:red)

    scatter!(plt,
        R_0[Int.((n/2+1):end), 1],
        R_0[Int.((n/2+1):end), 2],
        R_0[Int.((n/2+1):end), 3],
        markersize=0.5,
        markercolor=:red)
end

# Plot_DNA([R_0[:,3] R_0[:,1] R_0[:,2]])

# function InverseMatrix(K)
#     n = length(K[:, 1])
#     ΔK = 1
#     p = 0
#     C = zeros(n, n)
#     while ΔK > 1e-8 && p < 100
#         C .+= (-K)^p
#         ΔK = maximum(abs.(K^p))
#         p += 1
#     end
#     C, p, ΔK
# end


function InverseMatrix(K)
    n = length(K[:, 1])
    ΔK = 1
    p = 0
    C = zeros(n, n)
    k = maximum(abs.(K))
    while ΔK > 1e-8
        p += 1
        ΔK = k^p
    end
    C = sum((-K)^i for i in 0:p)
    C, p, ΔK
end

temp = Gershgorin(K)

m = InverseMatrix(temp)[1]

function WrapperFunction(K, limit)
    tempK = Gershgorin(K)
    tempK .= InverseMatrix(tempK)[1]
    λ = []
    x = []
    k = 0
    while k < limit
        temp = PowerMethod(tempK)
        push!(λ, temp[1])
        push!(x, temp[2])
        MatrixDeflation!(tempK, temp[2], temp[1])
        k += 1
        println("Calculating Eigenmode: $k ")
    end
    # λ = 1 .- λ
    λ = (1 ./ λ) .- 1
    λ, x
end
macro Name(arg)
   x = string(arg)
    quote
    $x
    end
end

tempλ = WrapperFunction(K, 10)[2]
tempλ2 = PMMD(Gershgorin(K), 10)[2]

half = Int.(1:(n/2))
half2 = Int.((n/2)+1:n)
function PlotVector(tempλ, name)
    p = []
    for i in eachindex(tempλ)
        plt = scatter(R_0[half, 3], tempλ[i][half],
            ticks=:none,
            legend=true,
            title="Eigenmode $i",
            marksize=1)
        scatter!(plt, R_0[half2, 3], tempλ[i][half2])
        push!(p, plt)
    end
    plt = p
    plt = plot(plt[1], plt[2], plt[3], plt[4], plt[5], plt[6], plt[7], plt[8], plt[9], plt[10],
        markersize=8,
        legend=false,
        layout=(5, 2),
        size=(2000, 2000))
    display(plt)

     png(name)
end
PlotVector(tempλ,"Niedrig")
