module CompChem

include("BasisModule.jl")
include("GeometryModule.jl")
include("IntegralModule.jl")
include("SCFModule.jl")

# Write your package code here.

function RHF_full(num_iterations, tolerance, basis_name, basis_file, structure_file)

# Set up basis and geometry
basis_type, basis = BasisModule.getBasisDict(basis_name, basis_file)
total_electrons, total_atoms, sys_name, atoms = GeometryModule.getAtomsFromFile(basis_type, structure_file)

# Set up all the orbitals
total_orbitals, cumulative_orbitals, orbitals = BasisModule.generateOrbitalsFromAtoms(atoms, basis, total_atoms, basis_type)
#print(orbitals)

# Do all the integrals
overlap_matrix, kinetic_matrix, nuclear_matrix, eri_matrix = IntegralModule.calIntegrals(total_orbitals, orbitals, atoms)

# Solve the Roothaan-Hall/Pople-Nesbit equation(s)
vectors, levels, total_energy = SCFModule.doSCF(num_iterations, overlap_matrix, kinetic_matrix, nuclear_matrix, 
eri_matrix, total_orbitals, total_electrons, atoms, total_atoms, tolerance, 1)

return vectors, levels, total_energy 

end


function RHF_energy(num_iterations, tolerance, basis_name, basis_file, structure_file)

# Set up basis and geometry
basis_type, basis = BasisModule.getBasisDict(basis_name, basis_file)
total_electrons, total_atoms, sys_name, atoms = GeometryModule.getAtomsFromFile(basis_type, structure_file)
    
# Set up all the orbitals
total_orbitals, cumulative_orbitals, orbitals = BasisModule.generateOrbitalsFromAtoms(atoms, basis, total_atoms, basis_type)
#print(orbitals)
    
# Do all the integrals
overlap_matrix, kinetic_matrix, nuclear_matrix, eri_matrix = IntegralModule.calIntegrals(total_orbitals, orbitals, atoms)
    
# Solve the Roothaan-Hall/Pople-Nesbit equation(s)
vectors, levels, total_energy = SCFModule.doSCF(num_iterations, overlap_matrix, kinetic_matrix, nuclear_matrix, 
eri_matrix, total_orbitals, total_electrons, atoms, total_atoms, tolerance, 0)
    
return total_energy 
    
end

end