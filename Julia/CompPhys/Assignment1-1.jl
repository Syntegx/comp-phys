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
N = 1000;
# using the distance d between electrodes and the number of intervalls N we can calculate the length of a single intervall h
h = d / N;
# x is a vector containing all equidistant steps x_i
x = collect(1:1:N) .* h;
# since its a linear differential equation it can be written in the form
# A y(x) = b where A describes a matrix A, y(x) is a vector containing the elements (y(x_1),...,y(x_N))
# Matrix a has form a_ii = h^2 k^2 + 2, a_i+1 j = -1 = a_i j+1
zero_N = zeros(N);
a_ii = k^2 * h^2 + 2;
A = zeros(N, N) + SymTridiagonal(zero_N .- a_ii, zero_N .+ 1);
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
plot(x, y_num);
plot!(x, DebyeHueckelFunction(x));

# Using LU Decomposition the upper and lower triangular matrices L and U can be found and used to solve the one dimensional Debye-Hueckel equation
function LUdecomp!(A)
    # first we determine the size nxn of the matrix A
    n = size(A, 1)
    A = zeros(n, n) + A
    # here we make sure our matrix contains at least zeros as every element, otherwise the allocation of new values wont work
    for j in 1:n # we loop over every column 1 to n
        # since the first row of our linear system contains only terms of trivial nature we can skip the first row
        for i in j+1:n
            A[i, j] = A[i, j] / A[j, j] # since all of the values u_1j are known, theres only a single unknown variables l_n1 in the first column
            for k in j+1:n # after calculating the unknown values l_n1 we can calculate the next unknown variable u_2m since all others are given 
                A[i, k] = A[i, k] - A[i, j] * A[j, k]
            end
        end
    end
    return A
end

# forward substitution 
function ForwardSubstitution!(L, b)
    n = size(b, 1)
    for i in 1:n # loop every row
        for k in 1:i-1 # loop every column except the last one using only indices lower than i
            b[i] = b[i] - L[i, k] * b[k]
        end
    end
    return b
end

# backwardsubstitution
function BackwardSubstitution!(U, b)
    n = size(b, 1)
    for i in n:-1:1 # loop every row
        for k in i+1:n # using indices including i and higher than i
            b[i] = (b[i] - U[i, k] * b[k])
        end
        b[i] = b[i] / U[i, i]
    end
    return b
end

# function using all other functions created for solving lin equation systems
function SolveLinSys(A, b)
    A = LUdecomp!(A)
    ForwardSubstitution!(A, b)
    BackwardSubstitution!(A, b)
    return b
end

function NegativeTotalCharge_q(y, k, h)
    q = ((2 * c_0) / k^2) .* (((y[N] - y[N-1]) ./ h) - ((y[2] - y[1]) ./ h))
    return q
end

Charge_Vec = zeros(N);
Cap_Vec = zeros(N);
d = collect(LinRange(1, 100, N))
y_0 = collect(LinRange(-5, -0.5, N))
    
for j in 1:N
    h = d[j] / N
    k = sqrt((2 * c_0 * e_0^2) / (epsilon_0 * epsilon_r * k_B * T)) * 1
    for i in 1:N
        local y_d = 0;
        local a_ii = k^2 * h^2 + 2
        local A = zeros(N, N) + SymTridiagonal(zero_N .+ a_ii, zero_N .- 1)
        local b = vcat(y_0[i], zeros(N - 2), y_d)
        SolveLinSys(A, b)
        Charge_Vec[i] = NegativeTotalCharge_q(b, k, h)
    end

    Cap_Vec[j] = (Charge_Vec[N] - Charge_Vec[1]) ./ (h*N)
end

plot(d, Charge_Vec)