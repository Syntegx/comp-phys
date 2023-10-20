using LinearAlgebra
### DEFINING FUNCTION RING USED TO INITIALIZE USED MATRICES ###
# N specifies the used dimension 
# t specifies the probability for an electron to move from state b to state a
# E is the energy of the system
# delta is the value of the energy broadening of entering/leaving the system
# alpha/beta are the sites of entering/leaving the system
function ring(N, t, E, delta, alpha, beta)
    Matrix = complex(zeros(N, N));
    Matrix .= SymTridiagonal(fill(-Float64(E), N), fill(t, N));
    Matrix[CartesianIndex.([1, N], [N, 1])] .= t;
    Matrix[CartesianIndex.([alpha, beta], [alpha, beta])] .+= delta * 1im;
    Matrix
end

function GaussSeidel(A, x, b)
    delta_x .= x
    for  in
        x[:] .-= delta_x[:] 
        delta_x[:] .+= sum.()
    
        x[:] .-= delta_x[:]
    end
end