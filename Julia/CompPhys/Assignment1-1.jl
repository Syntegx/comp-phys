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

k_sqr = (2 * c_0 * e_0^2) / (epsilon_0 * epsilon_r * k_B * T)

# since its a linear differential equation it can be written in the form
# A y(x) = b where A describes a matrix A, y(x) is a vector containing the elements (y(x_1),...,y(x_N))



# Using LU Decomposition the matrix A can found and used to solve the one dimensional Debye-Hueckel equation
# before implementing LU Decomposition by myself I will use to LinearAlgebra package to have a solution I can compare to

# following boundary conditions are given

y_0 = -.5;
y_d = .5;
