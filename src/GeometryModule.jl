using LinearAlgebra, SpecialFunctions

module GeometryModule

mutable struct Atom
    name::String
    atomic_number::Int64
    num_electrons::Int64
    pos::Array{Float64}
    num_orbitals::Int64
end

function getNumOrbitals(basis_type, num_electrons)
    num_orbitals = 0
    if (basis_type == "STO-NG")
      if (num_electrons <= 2)
         num_orbitals = 1
      elseif (num_electrons <= 10)
         num_orbitals = 5
      elseif (num_electrons <= 12)
         num_orbitals = 6
      elseif (num_orbitals <= 18)
         num_orbitals = 9
      end
    elseif (basis_type == "X-YZG")
        if (num_electrons <= 2)
            num_orbitals = 2
        elseif (num_electrons <= 10)
            num_orbitals = 9
        elseif (num_electrons <= 18)
            num_orbitals = 13
        end
    elseif (basis_type == "X-YZG**" || basis_type == "cc-pVDZ")
        if (num_electrons <= 2)
            num_orbitals = 5
        elseif (num_electrons <= 10)
            num_orbitals = 15
        elseif (num_electrons <= 18)
            num_orbitals = 18
        end
    end
    return num_orbitals
end

function getAtomicNumber(name)
    if (name == "H")
        return 1, 1
    elseif (name == "He")
        return 2, 2
    elseif (name == "C")
        return 6, 6
    elseif (name == "N")
        return 7, 7
    elseif (name == "O")
        return 8, 8
    elseif (name == "F")
        return 9, 9
    end
end

function getAtomsFromFile(basis_type, filename)
    f = open(filename, "r")
    posdata = readlines(f)
    close(f)
    num_atoms = parse(Int64, posdata[1])
    sys_name = posdata[2]
    Atoms = Array{Atom}(undef, num_atoms)
    for i in 1:num_atoms
        temp = split(posdata[2+i])
        atom_name = temp[1]
        atomic_number, num_electrons = getAtomicNumber(atom_name)
        pos = [parse(Float64, temp[2])/0.5291, parse(Float64, temp[3])/0.5291, parse(Float64, temp[4])/0.5291]
        num_orbitals = getNumOrbitals(basis_type, num_electrons)
        Atoms[i] = Atom(atom_name, atomic_number, num_electrons, pos, num_orbitals)
    end
    
    # Get the number of electrons (this has to be even, otherwise no point to doing RHF) 
    total_electrons = 0.0
    for i in 1:num_atoms
        total_electrons = total_electrons + Atoms[i].num_electrons
    end
    return total_electrons, num_atoms, sys_name, Atoms
end

function returnPositionVector(total_atoms, atoms)
    pos_vector = zeros(Float64, 3*total_atoms)
    for i in 1:total_atoms
        for j in 1:3
            pos_vector[3*(i-1)+j] = atoms[i].pos[j]
        end
    end
    return pos_vector
end

function modifyPosition(total_atoms, atoms, shift_vector)
    new_atoms = deepcopy(atoms)
    for i in 1:total_atoms
        for j in 1:3
            new_atoms[i].pos[j] = atoms[i].pos[j] + shift_vector[3*(i-1)+j]
        end
    end
    return new_atoms
end

function computeSVector(total_atoms, atoms, new_atoms)
    s_k = zeros(Float64, 3*total_atoms)
    for i in 1:total_atoms
        for j in 1:3
            s_k[3*(i-1)+j] = new_atoms[i].pos[j] - atoms[i].pos[j]
        end
    end
    return s_k
end

function computeYVector(total_atoms, old_gradient, new_gradient)
    y_k = zeros(Float64, 3*total_atoms)
    for i in 1:3*total_atoms
        y_k[i] = new_gradient[i] - old_gradient[i]
    end
    return y_k
end

end