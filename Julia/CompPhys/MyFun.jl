module MyFun
using LinearAlgebra

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
    L = zeros(n,n) + LowerTriangular(L)
    L[diagind(L)].=1

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
    U = zeros(n,n) + UpperTriangular(U)
   
   
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
    LUdecomp!(A)
    ForwardSubstitution!(A, b)
    BackwardSubstitution!(A, b)
    return b
end

end