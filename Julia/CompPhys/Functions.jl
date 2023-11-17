module Functions
function LUdecomp(A)
    # first we determine the size nxn of the matrix A
    n = size(A, 1)
    U = copy(A)
    # here we make sure our matrix contains at least zeros as every element, otherwise the allocation of new values wont work
    for j in 1:n # we loop over every column 1 to n
        # since the first row of our linear system contains only terms of trivial nature we can skip the first row
        for i in j+1:n
            U[i, j] = U[i, j] / U[j, j] # since all of the values u_1j are known, theres only a single unknown variables l_n1 in the first column
            for k in j+1:n # after calculating the unknown values l_n1 we can calculate the next unknown variable u_2m since all others are given 
                U[i, k] = U[i, k] - U[i, j] * U[j, k]
            end
        end
    end
    return U
end

# forward substitution 
function ForwardSubstitution(L, b)
    n = size(b, 1)
    ForwardVector = deepcopy(b);
    for i in 1:n # loop every row
        for k in 1:i-1 # loop every column except the last one using only indices lower than i
            ForwardVector[i] = ForwardVector[i] - L[i, k] * ForwardVector[k]
        end
    end
    return ForwardVector
end

# backwardsubstitution
function BackwardSubstitution(U, b)
    n = size(b, 1)
    BackwardVector = deepcopy(b);
    for i in n:-1:1 # loop every row
        for k in i+1:n # using indices including i and higher than i
            BackwardVector[i] = (BackwardVector[i] - U[i, k] * BackwardVector[k])
        end
        BackwardVector[i] = BackwardVector[i] / U[i, i]
    end
    return BackwardVector
end

# function using all other functions created for solving lin equation systems
function SolveLinSys(A, b)
    U = LUdecomp(A)
   ForwardVector = ForwardSubstitution(U, b)
   BackwardVector = BackwardSubstitution(U, ForwardVector)
    return BackwardVector
end
end