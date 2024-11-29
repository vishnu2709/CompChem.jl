using LinearAlgebra, SpecialFunctions, Plots

module BasisModule

mutable struct orbitalCG
    overall_norm::Float64
    pos::Array{Float64}
    o::Array{Int64}
    alphas::Array{Float64}
    coeffs::Array{Float64}
    norms::Array{Float64}
    num_gaussians::Int64
end


function getOrbitalTypeFromNumber(basis_type, orbital_number, atomic_number)
    if (basis_type == "STO-NG")
        if (orbital_number == 1)
            return "1s", [0, 0, 0]
        elseif (orbital_number == 2)
            return "2s", [0, 0, 0]
        elseif (orbital_number == 3)
            return "2p", [1, 0, 0]
        elseif (orbital_number == 4)
            return "2p", [0, 1, 0]
        elseif (orbital_number == 5)
            return "2p", [0, 0, 1]
        end
    elseif (basis_type == "X-YZG")
        if (atomic_number <= 2)
            if (orbital_number == 1)
                return "1s_1", [0,0,0]
            elseif (orbital_number == 2)
                return "1s_2", [0,0,0]
            end
        elseif (atomic_number <= 10)
            if (orbital_number == 1)
                return "1s", [0, 0, 0]
            elseif (orbital_number == 2)
                return "2s_1", [0, 0, 0]
            elseif (orbital_number == 3)
                return "2s_2", [0, 0, 0]
            elseif (orbital_number == 4)
                return "2p_1", [1, 0, 0]
            elseif (orbital_number == 5)
                return "2p_2", [1, 0, 0]
            elseif (orbital_number == 6)
                return "2p_1", [0, 1, 0]
            elseif (orbital_number == 7)
                return "2p_2", [0, 1, 0]
            elseif (orbital_number == 8)
                return "2p_1", [0, 0, 1]
            elseif (orbital_number == 9)
                return "2p_2", [0, 0, 1]
            end
        end
    elseif (basis_type == "X-YZG**" || basis_type == "cc-pVDZ")
        if (atomic_number <= 2)
            if (orbital_number == 1)
                return "1s_1", [0,0,0]
            elseif (orbital_number == 2)
                return "1s_2", [0,0,0]
            elseif (orbital_number == 3)
                return "1p", [1,0,0]
            elseif (orbital_number == 4)
                return "1p", [0, 1, 0]
            elseif (orbital_number == 5)
                return "1p", [0, 0, 1]
            end
        elseif (atomic_number <= 10)
            if (orbital_number == 1)
                return "1s", [0, 0, 0]
            elseif (orbital_number == 2)
                return "2s_1", [0, 0, 0]
            elseif (orbital_number == 3)
                return "2s_2", [0, 0, 0]
            elseif (orbital_number == 4)
                return "2p_1", [1, 0, 0]
            elseif (orbital_number == 5)
                return "2p_2", [1, 0, 0]
            elseif (orbital_number == 6)
                return "2p_1", [0, 1, 0]
            elseif (orbital_number == 7)
                return "2p_2", [0, 1, 0]
            elseif (orbital_number == 8)
                return "2p_1", [0, 0, 1]
            elseif (orbital_number == 9)
                return "2p_2", [0, 0, 1]
            elseif (orbital_number == 10)
                return "2d", [2, 0, 0]
            elseif (orbital_number == 11)
                return "2d", [0, 2, 0]
            elseif (orbital_number == 12)
                return "2d", [0, 0, 2]
            elseif (orbital_number == 13)
                return "2d", [1, 1, 0]
            elseif (orbital_number == 14)
                return "2d", [1, 0, 1]
            elseif (orbital_number == 15)
                return "2d", [0, 1, 1]
            end
        end
    end
end

function calNormGaussian(o, alpha)
    factorial_term = 1.0
    for i in 1:3
        for m in 0:(o[i]-1)
            factorial_term =  factorial_term*(-0.5 - m)
        end
    end
    factorial_term = factorial_term * pi^1.5 * (-1.0)^(o[1] + o[2] + o[3])
    norm = factorial_term * (2*alpha)^(-1.5 - o[1] - o[2] - o[3])
    return 1.0/sqrt(norm)
end

# Contracted Gaussian Norm
function calNormCG(orbital::orbitalCG)
    factorial_term = 1.0
    for i in 1:3
        for m in 0:(orbital.o[i]-1)
            factorial_term =  factorial_term*(-0.5 - m)
        end
    end
    factorial_term = factorial_term * pi^1.5 * (-1.0)^(orbital.o[1] + orbital.o[2] + orbital.o[3])
    contraction_sum = 0.0
    for i in 1:orbital.num_gaussians
        for j in 1:orbital.num_gaussians
            contraction_sum = contraction_sum + orbital.coeffs[i]*orbital.coeffs[j]*orbital.norms[i]*orbital.norms[j]*(orbital.alphas[i] 
            + orbital.alphas[j])^(-1.5 - orbital.o[1] - orbital.o[2] - orbital.o[3])
        end
    end
    norm = factorial_term * contraction_sum
    return 1.0/sqrt(norm)
end

function determineAtomofOrbital(orbital::orbitalCG, atoms, num_atoms)
    atom_index = 1
    for i in 1:num_atoms
        diff = norm(atoms[i].pos - orbital.pos)
        if (diff < 1e-3)
            atom_index = i 
            break
        end
    end
    return atom_index
end

function getBasisDict(basis_name, basis_file)
    if (contains(basis_name, "STO"))
        basis_type = "STO-NG"
        ng = parse(Int64, basis_name[5])
        basis = getBasisDict_STONG(basis_file, ng)
    elseif (basis_name == "6-31G")
        basis_type = "X-YZG"
        ncore = 6
        nvalence_1 = 3
        nvalence_2 = 1
        basis = getBasisDict_XYZG(basis_file, ncore, nvalence_1, nvalence_2)
    elseif (basis_name == "3-21G")
        basis_type = "X-YZG"
        ncore = 3
        nvalence_1 = 2
        nvalence_2 = 1
        basis = getBasisDict_XYZG(basis_file, ncore, nvalence_1, nvalence_2)
    elseif (basis_name == "6-31G(d,p)")
        basis_type = "X-YZG**"
        ncore = 6
        nvalence_1 = 3
        nvalence_2 = 1
        basis = getBasisDict_XYZGdp(basis_file, ncore, nvalence_1, nvalence_2)
    elseif (basis_name == "cc-pVDZ")
        basis_type = "cc-pVDZ"
        basis = getBasisDict_ccpvdz(basis_file)
    else
        basis = Dict()
        basis_type = nothing
    end
    return basis_type, basis
end

function ccpvdz_helper_function_light(basisdata, num_basis_lines, name)
    nvalence_1 = 4
    nvalence_2 = 1
    for line in 15:num_basis_lines
        if (contains(basisdata[line], name))
            data1s = basisdata[line+1:line+nvalence_1]
            exponent1s = zeros(Float64, nvalence_1)
            coeff1s_1 = zeros(Float64, nvalence_1)
            coeff1s_2 = zeros(Float64, nvalence_1)
            for i in 1:nvalence_1
                temp = split(data1s[i])
                exponent1s[i] = parse(Float64, temp[1])
                coeff1s_1[i] = parse(Float64, temp[2])
                coeff1s_2[i] = parse(Float64, temp[3])
            end
            data1p = basisdata[line+nvalence_1+2:line+nvalence_1+nvalence_2+1]
            exponent1p = zeros(Float64, nvalence_2)
            coeff1p = zeros(Float64, nvalence_2)
            temp = split(data1p[1])
            exponent1p[1] = parse(Float64, temp[1])
            coeff1p[1] = parse(Float64, temp[2])
            return [[exponent1s, coeff1s_1], [exponent1s, coeff1s_2], [exponent1p, coeff1p]]
        end
    end         
end

function ccpvdz_helper_function_heavy(basisdata, num_basis_lines, name)
    ncore = 9
    nvalence_1 = 4
    nvalence_2 = 1
    for line in 15:num_basis_lines
        if (contains(basisdata[line], name))
            data_s = basisdata[line+1:line+ncore]
            exponent1s = zeros(Float64, ncore)
            coeff1s = zeros(Float64, ncore)
            coeff2s_1 = zeros(Float64, ncore)
            coeff2s_2 = ones(Float64, nvalence_2)
            for i in 1:ncore
                temp = split(data_s[i])
                exponent1s[i] = parse(Float64, temp[1])
                coeff1s[i] = parse(Float64, temp[2])
                coeff2s_1[i] = parse(Float64, temp[3])
                #if (i == ncore)
                #    coeff2s_2[nvalence_2] = parse(Float64, temp[4])
                #end
            end
            data_p = basisdata[line+ncore+2:line+ncore+nvalence_1+1]
            exponent2p = zeros(Float64, nvalence_1)
            coeff2p_1 = zeros(Float64, nvalence_1)
            coeff2p_2 = ones(Float64, nvalence_2)
            for i in 1:nvalence_1
                temp = split(data_p[i])
                exponent2p[i] = parse(Float64, temp[1])
                coeff2p_1[i] = parse(Float64, temp[2])
                #if (i == nvalence_1)
                #    coeff2p_2[nvalence_2] = parse(Float64, temp[3])
                #end
            end
            data_d = basisdata[line+ncore+nvalence_1+3:line+ncore+nvalence_1+3]
            temp = split(data_d[1])
            exponent2d = [parse(Float64, temp[1])]
            coeff2d = [parse(Float64, temp[2])]
            return [[exponent1s, coeff1s], [exponent1s, coeff2s_1], [[exponent1s[ncore]], coeff2s_2], [exponent2p, coeff2p_1], [[exponent2p[nvalence_1]], coeff2p_2], [exponent2d, coeff2d]]
        end
    end
end

function STONG_helper_function_light(basisdata, num_basis_lines, name, ng)
    for line in 15:num_basis_lines
        if (contains(basisdata[line], name))
            dataH1s = basisdata[line+1:line+ng]
            coeffs = zeros(Float64, ng)
            exponents = zeros(Float64, ng)
            for i in 1:ng
                temp = split(dataH1s[i])
                exponents[i] = parse(Float64, temp[1])
                coeffs[i] = parse(Float64, temp[2])
            end
            return [[exponents, coeffs]]
        end
    end
end

function STONG_helper_function_heavy(basisdata, num_basis_lines, name, ng)
    for line in 15:num_basis_lines
        if (contains(basisdata[line], name))
            dataO1s = basisdata[line+1:line+ng]
            coeff1s = zeros(Float64, ng)
            coeff2s = zeros(Float64, ng)
            coeff2p = zeros(Float64, ng)
            exponent1s = zeros(Float64, ng)
            exponent2s = zeros(Float64, ng)
            for i in 1:ng
                temp = split(dataO1s[i])
                exponent1s[i] = parse(Float64, temp[1])
                coeff1s[i] = parse(Float64, temp[2])
            end
            #basis["O1s"] = [exponent1s, coeff1s]
            dataO2sp = basisdata[line+ng+2:line+ng+ng+1]
            for i in 1:ng
                temp = split(dataO2sp[i])
                exponent2s[i] = parse(Float64, temp[1])
                coeff2s[i] = parse(Float64, temp[2])
                coeff2p[i] = parse(Float64, temp[3])
            end
            #basis["O2s"] = [exponent2s, coeff2s]
            #basis["O2p"] = [exponent2s, coeff2p]
            return [[exponent1s, coeff1s], [exponent2s, coeff2s], [exponent2s, coeff2p]]
        end
    end
end

function XYZG_helper_function_light(basisdata, num_basis_lines, name, nvalence_1, nvalence_2, is_polar)
    for line in 15:num_basis_lines
    if (contains(basisdata[line], name))
        dataH1s_1 = basisdata[line+1:line+nvalence_1]
        coeff1s_1 = zeros(Float64, nvalence_1)
        exponent1s_1 = zeros(Float64, nvalence_1)
        coeff1s_2 = zeros(Float64, nvalence_2)
        exponent1s_2 = zeros(Float64, nvalence_2)
        for i in 1:nvalence_1
            temp = split(dataH1s_1[i])
            exponent1s_1[i] = parse(Float64, temp[1])
            coeff1s_1[i] = parse(Float64, temp[2])
        end
        #basis["H1s_1"] = [exponent1s_1, coeff1s_1]
        dataH1s_2 = basisdata[line+nvalence_1+2:line+nvalence_1+nvalence_2+1]
        for i in 1:nvalence_2
            temp = split(dataH1s_2[i])
            exponent1s_2[i] = parse(Float64, temp[1])
            coeff1s_2[i] = parse(Float64, temp[2])
        end
        #basis["H1s_2"] = [exponent1s_2, coeff1s_2]
        if (is_polar == 1)
            temp = split(basisdata[line+nvalence_1+nvalence_2+3:line+nvalence_1+nvalence_2+3][1])
            exponent1p = [parse(Float64, temp[1])]
            coeff1p = [parse(Float64, temp[2])]
            #basis["H1p"] = [exponent1p, coeff1p]
            return [[exponent1s_1, coeff1s_1], [exponent1s_2, coeff1s_2], [exponent1p, coeff1p]]
        else
            return [[exponent1s_1, coeff1s_1], [exponent1s_2, coeff1s_2]]
        end
    end
    end
end

function XYZG_helper_function_heavy(basisdata, num_basis_lines, name, ncore, nvalence_1, nvalence_2, is_polar)
    for line in 15:num_basis_lines
    if (contains(basisdata[line], name))
        data1s = basisdata[line+1:line+ncore]
        #println(data1s)
        coeff1s = zeros(Float64, ncore)
        exponent1s = zeros(Float64, ncore)
        for i in 1:ncore
            temp = split(data1s[i])
            exponent1s[i] = parse(Float64, temp[1])
            coeff1s[i] = parse(Float64, temp[2])
        end
        #basis["N1s"] = [exponent1s, coeff1s]
        dataN2sp_1 = basisdata[line+ncore+2:line+ncore+nvalence_1+1]
        coeff2s_1 = zeros(Float64, nvalence_1)
        coeff2p_1 = zeros(Float64, nvalence_1)
        exponent2sp_1 = zeros(Float64, nvalence_1)
        for i in 1:nvalence_1
            temp = split(dataN2sp_1[i])
            exponent2sp_1[i] = parse(Float64, temp[1])
            coeff2s_1[i] = parse(Float64, temp[2])
            coeff2p_1[i] = parse(Float64, temp[3])
        end
        #basis["N2s_1"] = [exponent2sp_1, coeff2s_1]
        #basis["N2p_1"] = [exponent2sp_1, coeff2p_1]
        dataN2sp_2 = basisdata[line+ncore+nvalence_1+3:line+ncore+nvalence_1+nvalence_2+2]
        coeff2s_2 = zeros(Float64, nvalence_2)
        coeff2p_2 = zeros(Float64, nvalence_2)
        exponent2sp_2 = zeros(Float64, nvalence_2)
        for i in 1:nvalence_2
            temp = split(dataN2sp_2[i])
            exponent2sp_2[i] = parse(Float64, temp[1])
            coeff2s_2[i] = parse(Float64, temp[2])
            coeff2p_2[i] = parse(Float64, temp[3])
        end
        if (is_polar == 1)
            temp = split(basisdata[line+ncore+nvalence_1+nvalence_2+4:line+ncore+nvalence_1+nvalence_2+4][1])
            exponent2d = [parse(Float64, temp[1])]
            coeff2d = [parse(Float64, temp[2])]
            return [[exponent1s, coeff1s], [exponent2sp_1, coeff2s_1], [exponent2sp_1, coeff2p_1], [exponent2sp_2, coeff2s_2], [exponent2sp_2, coeff2p_2], [exponent2d, coeff2d]]
        else
            return [[exponent1s, coeff1s], [exponent2sp_1, coeff2s_1], [exponent2sp_1, coeff2p_1], [exponent2sp_2, coeff2s_2], [exponent2sp_2, coeff2p_2]]
        end
    end
    end
end

# Function that reads 6-31G type basis files
function getBasisDict_XYZG(basis_file, ncore, nvalence_1, nvalence_2)
    f = open(basis_file, "r")
    basisdata = readlines(f)
    num_basis_lines = length(basisdata)
    basis = Dict()
    close(f)
    light_names = ["H", "He"]
    for name in light_names
        basis_tuples = XYZG_helper_function_light(basisdata, num_basis_lines, name, nvalence_1, nvalence_2, 0)
        basis[name*"1s_1"] = basis_tuples[1]
        basis[name*"1s_2"] = basis_tuples[2]
    end

    heavy_names = ["Li", "Be", "B", "C", "N", "O", "F", "Ne"]
    for name in heavy_names
        basis_tuples = XYZG_helper_function_heavy(basisdata, num_basis_lines, name*"    S", ncore, nvalence_1, nvalence_2, 0)
        basis[name*"1s"] = basis_tuples[1]
        basis[name*"2s_1"] = basis_tuples[2]
        basis[name*"2p_1"] = basis_tuples[3]
        basis[name*"2s_2"] = basis_tuples[4]
        basis[name*"2p_2"] = basis_tuples[5]
    end
    return basis
end

# Function that reads 6-31G type basis files
function getBasisDict_ccpvdz(basis_file)
    f = open(basis_file, "r")
    basisdata = readlines(f)
    num_basis_lines = length(basisdata)
    basis = Dict()
    close(f)
    light_names = ["H", "He"]
    for name in light_names
        basis_tuples = ccpvdz_helper_function_light(basisdata, num_basis_lines, name)
        basis[name*"1s_1"] = basis_tuples[1]
        basis[name*"1s_2"] = basis_tuples[2]
        basis[name*"1p"] = basis_tuples[3]
    end

    heavy_names = ["Li", "Be", "B", "C", "N", "O", "F", "Ne"]
    for name in heavy_names
        basis_tuples = ccpvdz_helper_function_heavy(basisdata, num_basis_lines, name*"    S")
        basis[name*"1s"] = basis_tuples[1]
        basis[name*"2s_1"] = basis_tuples[2]
        basis[name*"2s_2"] = basis_tuples[3]
        basis[name*"2p_1"] = basis_tuples[4]
        basis[name*"2p_2"] = basis_tuples[5]
        basis[name*"2d"] = basis_tuples[6]
    end
    return basis
end

# Function that reads 6-31G(d,p) type basis files
function getBasisDict_XYZGdp(basis_file, ncore, nvalence_1, nvalence_2)
    f = open(basis_file, "r")
    basisdata = readlines(f)
    num_basis_lines = length(basisdata)
    basis = Dict()
    close(f)
    light_names = ["H", "He"]
    for name in light_names
        basis_tuples = XYZG_helper_function_light(basisdata, num_basis_lines, name, nvalence_1, nvalence_2, 1)
        basis[name*"1s_1"] = basis_tuples[1]
        basis[name*"1s_2"] = basis_tuples[2]
        basis[name*"1p"] = basis_tuples[3]
    end

    heavy_names = ["Li", "Be", "B", "C", "N", "O", "F", "Ne"]
    for name in heavy_names
        basis_tuples = XYZG_helper_function_heavy(basisdata, num_basis_lines, name*"    S", ncore, nvalence_1, nvalence_2, 1)
        basis[name*"1s"] = basis_tuples[1]
        basis[name*"2s_1"] = basis_tuples[2]
        basis[name*"2p_1"] = basis_tuples[3]
        basis[name*"2s_2"] = basis_tuples[4]
        basis[name*"2p_2"] = basis_tuples[5]
        basis[name*"2d"] = basis_tuples[6]
    end
    return basis
end

# Function that reads STO-NG basis files
function getBasisDict_STONG(basis_file, ng)
    f = open(basis_file, "r")
    basisdata = readlines(f)
    num_basis_lines = length(basisdata)
    basis = Dict()
    close(f)
    light_names = ["H", "He"]
    for name in light_names
        basis_tuples = STONG_helper_function_light(basisdata, num_basis_lines, name, ng)
        basis[name*"1s"] = basis_tuples[1]
    end

    heavy_names = ["Li", "Be", "B", "C", "N", "O", "F", "Ne"]
    for name in heavy_names
        basis_tuples = STONG_helper_function_heavy(basisdata, num_basis_lines, name*"    S", ng)
        basis[name*"1s"] = basis_tuples[1]
        basis[name*"2s"] = basis_tuples[2]
        basis[name*"2p"] = basis_tuples[3]
    end
    return basis
end

function generateOrbitalsFromAtoms(Atoms, basis, num_atoms, basis_type)
    # In a minimal basis, there is one Contracted Gaussian per orbital. First count the total number of orbitals

    total_orbitals = Atoms[1].num_orbitals
    cumulative_orbitals = zeros(Int64, num_atoms)
    for i in 2:num_atoms
        total_orbitals = total_orbitals + Atoms[i].num_orbitals
        cumulative_orbitals[i] = cumulative_orbitals[i-1] + Atoms[i-1].num_orbitals
    end
    
    # Create an array of orbital types
    orbitals = Array{orbitalCG}(undef, total_orbitals)
    for i in 1:num_atoms
        atom_name = Atoms[i].name
        atomic_number = Atoms[i].atomic_number
        position = Atoms[i].pos
        for j in 1:Atoms[i].num_orbitals
            otype, o_r = getOrbitalTypeFromNumber(basis_type, j, atomic_number)
            atom_orbital_name = atom_name*otype
            basis_info = basis[atom_orbital_name]
            ng = length(basis_info[1])
            norms = ones(Float64, ng)
            overall_norm = 1 #temporary
            orbitals[cumulative_orbitals[i] + j] = orbitalCG(overall_norm, position, o_r, basis_info[1], basis_info[2], norms, ng)
            for k in 1:ng
                orbitals[cumulative_orbitals[i] + j].norms[k] = calNormGaussian(orbitals[cumulative_orbitals[i] + j].o, 
                orbitals[cumulative_orbitals[i]+j].alphas[k])
            end
            orbitals[cumulative_orbitals[i] + j].overall_norm = calNormCG(orbitals[cumulative_orbitals[i] + j])
        end
    end
    #for i in 1:total_orbitals
    #    println(orbitals[i])
    #end
    return total_orbitals, cumulative_orbitals, orbitals
end

# Getting real space functions for the molecular orbitals
function contracted_gaussian(r_vec, pos, o, alphas, coefficients, norms, ng)
    cg = 0.0
    diff_vec = r_vec - pos
    diff = norm(r_vec - pos)
    for i in 1:ng
        cg = cg + norms[i]*coefficients[i]*exp(-alphas[i]*diff^2)
    end
    return cg
end

function molecular_orbital(r_vec, converged_coefficients, orbitals, num_atomic_orbitals)
    total = 0.0
    for i in 1:num_atomic_orbitals
        total = total + converged_coefficients[i]*contracted_gaussian(r_vec, orbitals[i].pos, orbitals[i].pos, orbitals[i].alphas, 
        orbitals[i].coeffs, orbitals[i].norms, orbitals[i].num_gaussians)
    end
    return total
end

end