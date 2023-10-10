using LinearAlgebra, Plots

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
A = SymTridiagonal(zero_N .+ a_ii, zero_N .- 1);
# using the boundary conditions y_0 and y_d we can solve this system of linear equations using \
y_0 = -0.5;
y_d = 0.5;
b = vcat(y_0,zeros(N-2),y_d);
# this results in the numerical solution of the Debye-Hueckel equation(DHE)
y_num = A\b;
# for comparision we calculate the analytical solution to the DHE
y_ana = (y_d.* sinh.(k.*x) .+ y_0.* sinh.(k.*(d.-x))) ./ sinh(k*d) ;
plot(x,y_num)
plot!(x,y_ana)
# Using LU Decomposition the matrix A can found and used to solve the one dimensional Debye-Hueckel equation
function LUdecomp(A)
    n = size(A,1)
    for i in 1:n #columns of A
        for j in 1:(i-1)
            A[[i,j]] = A[[i,j]] - 
        end
    end
    return LU
end