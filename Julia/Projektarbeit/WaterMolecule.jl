using Unitful
#Defining Struct storing some physical properties of Atoms
mutable struct Atom
Dichte::Float64
Radius::Float64
Masse::Float64
end

#Defining variable Wasserstoff containing density, mass, radius of the Atom
Wasserstoff = Atom(0.0709,25,1.008);


function ()
    
end