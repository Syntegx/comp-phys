using LinearAlgebra, Plots, BenchmarkTools

e_0 = 16e-20 #elementary charge in coloumb
epsilon_0 = 885e-23 #permittivity of vacuum in farad per nanometer
k_B = 138e-25 #Boltzmann constant in joule per kelvin
c_0 = 6e-3 # concentration of ions per square nm at zero potential
epsilon_r = 80 # the relative dielectric constant of the simulated medium
T = 300 # the temperature of the simulated medium in kelvin
d = 100 # the width of the electrolyte solution in nm


# The Debye-Hueckel equation has the following form
# d^2y(x)/d^2x = k^2 y(x)
# the inverse Debye length k^2 is a constant defined by following term
k = sqrt((2 * c_0 * e_0^2) / (epsilon_0 * epsilon_r * k_B * T)) * 1
# this length is to be divided into a grid of equidistant intervalls N
N = 200;
# using the distance d between electrodes and the number of intervalls N we can calculate the length of a single intervall h
h = d / N;
# x is a vector containing all equidistant steps x_i
x = collect(1:1:N) .* h;
# since its a linear differential equation it can be written in the form
# A y(x) = b where A describes a matrix A, y(x) is a vector containing the elements (y(x_1),...,y(x_N))
# Matrix a has form a_ii = h^2 k^2 + 2, a_i+1 j = -1 = a_i j+1

zero_N = zeros(N);
a_ii = k^2 * h^2 + 2;
A = zeros(N, N) + SymTridiagonal(zero_N .+ a_ii, zero_N .- 1);
# using the boundary conditions y_0 and y_d we can solve this system of linear equations using \
y_0 = -0.5;
y_d = 0.5;
b = vcat(y_0, zeros(N - 2), y_d);
# this results in the numerical solution of the Debye-Hueckel equation(DHE)
y_num = A \ b;
# for comparision we calculate the analytical solution to the DHE
function DebyeHueckelFunction(x, y_0=-0.5, y_d=0.5, d=100)
    y_ana = (y_d .* sinh.(k .* x) .+ y_0 .* sinh.(k .* (d .- x))) ./ sinh(k * d)
    return y_ana
end

plt = plot(x, y_num,
    titlefontsize = 12,
    title="Debye-Hueckel Equation",
    xlabel="Position x/nm",
    ylabel="Electrostatic potential y",
    label="Numerical Solution"
)
plot!(x, DebyeHueckelFunction(x),
    label="Analytical Solution")

# Using LU Decomposition the upper and lower triangular matrices L and U can be found and used to solve the one dimensional Debye-Hueckel equation
function LUdecomp(A)
    # first we determine the size nxn of the matrix A
    n = size(A, 1)
    U = copy(A)
    # here we make sure our matrix contains at least zeros as every element, otherwise the allocation of new values wont work
    for j in 1:n # we loop over every column 1 to n
        # since the first row of our linear system contains only terms of trivial nature we can skip the first row
        @simd for i in j+1:n
            @inbounds U[i, j] = U[i, j] / U[j, j] # since all of the values u_1j are known, theres only a single unknown variables l_n1 in the first column
            @simd for k in j+1:n # after calculating the unknown values l_n1 we can calculate the next unknown variable u_2m since all others are given 
                @inbounds U[i, k] = U[i, k] - U[i, j] * U[j, k]
            end
        end
    end
    return U
end

# forward substitution 
function ForwardSubstitution(L, b)
    n = size(b, 1)
    ForwardVector = deepcopy(b)
    for i in 1:n # loop every row
        @simd for k in 1:i-1 # loop every column except the last one using only indices lower than i
            @inbounds ForwardVector[i] = ForwardVector[i] - L[i, k] * ForwardVector[k]
        end
    end
    return ForwardVector
end

# backward substitution
function BackwardSubstitution(U, b)
    n = size(b, 1)
    BackwardVector = deepcopy(b)
    @simd for i in n:-1:1 # loop every row
        @simd for k in i+1:n # using indices including i and higher than i
            @inbounds BackwardVector[i] = (BackwardVector[i] - U[i, k] * BackwardVector[k])
        end
        @inbounds BackwardVector[i] = BackwardVector[i] / U[i, i]
    end
    return BackwardVector
end

# wrapper function using all other functions created for solving lin equation systems
function SolveLinSys(A, b)
    U = LUdecomp(A)
    ForwardVector = ForwardSubstitution(U, b)
    BackwardVector = BackwardSubstitution(U, ForwardVector)
    return BackwardVector
end
plot!(x, SolveLinSys(A, b),
    label="LU Solution")
png(plt, "Debye-Hueckel Equation")
# function that takes the first numerical derivate using the finite difference method
function FirstDerivative(y, x, h)
    Derivative = (y[x] - y[x-1]) / h
    return Derivative
end

# function using the given equation to calculate the total negative charge of the capacitor
function NegativeTotalCharge(y, k, h)
    q = ((2 * c_0) / k^2) .* (FirstDerivative(y, length(y), h) - FirstDerivative(y, 2, h))
    return q
end

### Defining a function to calculate the capacitance per unit area of the supercapacitor###
# inputs are variable properties of the capacitator such as
# temperature,
# epsilon of electrolyte,
# device width,
# and a boolean which is used to switch the function between a mode which variates device width d,
# and a mode which variated temperature T and prints the chosen plot of capacitance as function of d/t

function Capacitance(
    bool=true,
    c_0=6e-3,
    epsilon_r=80,
    T=300,
    d=100,
    Total_Its=100)

    e_0 = 16e-20 #elementary charge in coloumb
    epsilon_0 = 885e-23 #permittivity of vacuum in farad per nanometer
    k_B = 138e-25 #Boltzmann constant in joule per kelvin
    # initializing storage arrays
    Charge_Vec = zeros(Total_Its)
    Cap_Vec = zeros(Total_Its)
    y_0 = collect(LinRange(-5, -0.5, Total_Its))
    if bool
        d = collect(LinRange(1, d, Total_Its))
    else
        T = LinRange(T - 300, T + 300, Total_Its)
    end


    for j in 1:Total_Its
        if bool
            local h = d[j] / N
            local k = sqrt((2 * c_0 * e_0^2) / (epsilon_0 * epsilon_r * k_B * T))
        else
            local h = d / N
            local k = sqrt.((2 * c_0 * e_0^2) ./ (epsilon_0 * epsilon_r * k_B .* T[j]))
        end

        for i in 1:Total_Its
            local y_d = 0
            if bool
                local a_ii = k .^ 2 * h^2 + 2
            else
                local a_ii = k .^ 2 * h^2 + 2
            end
            local A = zeros(N, N) + SymTridiagonal(zero_N .+ a_ii, zero_N .- 1)
            local b = vcat(y_0[i], zeros(N - 2), y_d)
            y = SolveLinSys(A, b)
            Charge_Vec[i] = NegativeTotalCharge(y, k, h)
        end
        Cap_Vec[j] = sum(Charge_Vec ./ y_0) / Total_Its
    end

    plt1 = plot(y_0, Charge_Vec,
        titlefontsize = 12,
        title="Capacitance per unit area",
        ylabel="Charge density dq",
        xlabel="Change of voltage y(0)",
        label = "dq/dy(0)"
    )
    png(plt1,"Capacitance")
    if bool
        plt2 = plot(d, Cap_Vec,
            titlefontsize = 12,
            title="Capacitance per unit area as function of device width d",
            ylabel="C(d) / F",
            xlabel="d / nm",
            label = "C(d)"
            )
            png(plt2,"Capacitance per device width")
    else
        plt2 = plot(T, Cap_Vec,
        titlefontsize = 12,
            title="Capacitance per unit area as function of temperature",
            ylabel="C(T) / F",
            xlabel="T / K",
            label = "C(T)"
            )
            png(plt2,"Capacitance per temperature")
    end
    display(plt1)
    display(plt2)
end
Capacitance(true)
Capacitance(false)
