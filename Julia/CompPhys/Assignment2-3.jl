using CSV, DataFrames, LinearAlgebra, MyFunctions

R_0 = CSV.read("C:\\Git\\Julia\\comp-phys\\Julia\\Data\\xyzm_dna.txt", DataFrame, header=false, delim=',')
M = R_0[:,4]
R_0 = R_0[:,1:3]
