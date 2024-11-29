include("IntegralModule.jl")
using .IntegralModule

module LocalizationModule
# Edmiston-Ruedenberg localization


# Calculate <psi_a psi_b | g | psi_c psi_d>
# given that we already have the integrated elements in terms of basis functions
function calculateInteractionMatrixElement(vector_a, vector_b, vector_c, vector_d, eri_matrix, total_orbitals)
    value = 0.0
    for i in 1:total_orbitals
        for j in 1:total_orbitals
            for k in 1:total_orbitals
                for l in 1:total_orbitals
                    value = value + vector_a[i]*vector_b[j]*vector_c[k]*vector_d[l]*eri_matrix[i,j,k,l]
                end
            end
        end
    end
    return value
end

# calculate Aij as highlighted in paper
function calculateAij(vector_i, vector_j, eri_matrix, total_orbitals)
    # <ii | g | jj >
    term_1 = calculateInteractionMatrixElement(vector_i, vector_i, vector_j ,vector_j, eri_matrix, total_orbitals)
    
    # <ii | g | ii >
    term_2 = calculateInteractionMatrixElement(vector_i, vector_i, vector_i, vector_i, eri_matrix, total_orbitals)
    
    # <ij | g | ij >
    term_3 = calculateInteractionMatrixElement(vector_i, vector_j, vector_i, vector_j, eri_matrix, total_orbitals)

    # <ji | g | ji >
    term_4 = calculateInteractionMatrixElement(vector_j, vector_i, vector_j, vector_i, eri_matrix, total_orbitals)

    # <jj | g || jj >
    term_5 = calculateInteractionMatrixElement(vector_j, vector_j, vector_j, vector_j, eri_matrix, total_orbitals)

    return term_1 - 0.25* (term_2 - term_3 - term_4 + term_5)
end

# calculate Bij as highlighted in paper
function calculateBij(vector_i, vector_j, eri_matrix, total_orbitals)
    # < ii | g | ij>
    term_1 = calculateInteractionMatrixElement(vector_i, vector_i, vector_i, vector_j, eri_matrix, total_orbitals)

    # < ji | g | jj >
    term_2 = calculateInteractionMatrixElement(vector_j, vector_i, vector_j, vector_j, eri_matrix, total_orbitals)

    return term_1 - term_2
end

# Calculate difference between Dmax and D for a given pair
function calculateDiff(vector_i, vector_j, eri_matrix, total_orbitals)
    Aij = calculateAij(vector_i, vector_j, eri_matrix, total_orbitals)
    Bij = calculateBij(vector_i, vector_j, eri_matrix, total_orbitals)
    return Aij + sqrt(Aij^2 + Bij^2)
end

# Calculate new pair from old pair
function calNewPair(vector_i, vector_j, eri_matrix, total_orbitals)
    Aij = calculateAij(vector_i, vector_j, eri_matrix, total_orbitals)
    Bij = calculateBij(vector_i, vector_j, eri_matrix, total_orbitals)
    if (abs(atan(-Bij/Aij)) <= 1e-5)
        alpha = 0.25*pi
    else
        alpha = 0.25*atan(-Bij/Aij)
        if (alpha < 0)
            alpha = alpha + pi/2
        end
    end
    println("Alpha in degrees: ", (alpha/pi) * 180)
    #println(Aij,Bij,alpha)
    new_vector_i = zeros(Float64, total_orbitals)
    new_vector_j = zeros(Float64, total_orbitals)
    for i in 1:total_orbitals
        new_vector_i[i] = cos(alpha)*vector_i[i] + sin(alpha)*vector_j[i]
        new_vector_j[i] = - sin(alpha)*vector_i[i] + cos(alpha)*vector_j[i] 
    end
    return new_vector_i, new_vector_j
end

function normalizationTest(vectors, overlap_matrix, total_orbitals, num_occupied_orbitals)
    normalization_test = zeros(Float64, num_occupied_orbitals)
    for i in 1:num_occupied_orbitals
        norm = 0
        for j in 1:total_orbitals
            for k in 1:total_orbitals
                norm = norm + vectors[i,j]*overlap_matrix[j,k]*vectors[i,k]
            end
        end
        normalization_test[i] = norm
    end
    return normalization_test
end

# Calculate D
function calculateD(vectors, eri_matrix, num_occupied_orbitals, total_orbitals)
    diag = 0.0
    for i in 1:num_occupied_orbitals
        diag = diag + calculateInteractionMatrixElement(vectors[i,:], vectors[i,:], vectors[i,:], vectors[i,:], eri_matrix, total_orbitals)
    end
    return diag
end

# Calculate for all pairs
function performLocalization(vectors, eri_matrix, num_occupied_orbitals, total_orbitals)
    
    max_val = 1000
    while (max_val > 1e-3)
    max_tuple = [0,0]
    max_val = 0
    for i in 1:num_occupied_orbitals
        for j in i+1:num_occupied_orbitals
            diff_ij = calculateDiff(vectors[i,:], vectors[j,:], eri_matrix, total_orbitals)
            if (diff_ij > max_val)
                max_val = diff_ij
                max_tuple = [i,j]
            end
            #println(diff_ij)
        end
    end
    vectors[max_tuple[1],:], vectors[max_tuple[2],:] = calNewPair(vectors[max_tuple[1],:], vectors[max_tuple[2],:], eri_matrix, total_orbitals)
    #println("New diff:",calculateDiff(vectors[max_tuple[1],:], vectors[max_tuple[2],:], eri_matrix, total_orbitals))
    println("Max val:",max_val)
    end
    finalD = calculateD(vectors, eri_matrix, num_occupied_orbitals, total_orbitals)
    println("finalD:",finalD)
    return vectors
end

end