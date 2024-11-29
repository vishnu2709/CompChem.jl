
module SCFModule
using LinearAlgebra, SpecialFunctions
include("GeometryModule.jl")
include("BasisModule.jl")
include("IntegralModule.jl")
import .GeometryModule
import .BasisModule
import .IntegralModule

function constructDensities_UHF(C_up, C_down, total_orbitals, num_occupied_up_orbitals, num_occupied_down_orbitals)
    P_up = zeros(Float64, (total_orbitals, total_orbitals))
    P_down = zeros(Float64, (total_orbitals, total_orbitals))
    for i in 1:total_orbitals
        for j in 1:total_orbitals
            for k in 1:num_occupied_up_orbitals
                P_up[i,j] = P_up[i,j] + C_up[k,i]*C_up[k,j]
            end
            for k in 1:num_occupied_down_orbitals
                P_down[i,j] = P_down[i,j] + C_down[k,i]*C_down[k,j]
            end
        end
    end
    return P_up, P_down
end

function buildFockMatrices_UHF(density_up, density_down, kinetic_matrix, nuclear_matrix, eri_matrix, total_orbitals)
    fock_up = zeros(Float64, (total_orbitals, total_orbitals))
    fock_down = zeros(Float64, (total_orbitals, total_orbitals))
    for i in 1:total_orbitals
        for j in 1:total_orbitals
            fock_up[i,j] = kinetic_matrix[i,j] + nuclear_matrix[i,j]
            fock_down[i,j] = kinetic_matrix[i,j] + nuclear_matrix[i,j]
            for k in 1:total_orbitals
                for l in 1:total_orbitals
                    fock_up[i,j] = fock_up[i,j] + density_up[k,l]*(eri_matrix[i,k,j,l] - eri_matrix[i,k,l,j]) + density_down[k,l]*eri_matrix[i,k,j,l]
                    fock_down[i,j] = fock_down[i,j] + density_down[k,l]*(eri_matrix[i,k,j,l] - eri_matrix[i,k,l,j]) + density_up[k,l]*eri_matrix[i,k,j,l]
                end
            end
        end
    end
    return fock_up, fock_down
end

function constructDensity(C, total_orbitals, num_occupied_orbitals)
    P = zeros(Float64, (total_orbitals, total_orbitals))
    for i in 1:total_orbitals
        for j in 1:total_orbitals
            for k in 1:num_occupied_orbitals
                P[i,j] = P[i,j] + 2*C[k,i]*C[k,j]
            end
        end
    end
    return P
end

function buildFockMatrix(density, kinetic_matrix, nuclear_matrix, eri_matrix, total_orbitals)
    fock = zeros(Float64, (total_orbitals, total_orbitals))
    for i in 1:total_orbitals
        for j in 1:total_orbitals
            fock[i,j] = kinetic_matrix[i,j] + nuclear_matrix[i,j]
            for k in 1:total_orbitals
                for l in 1:total_orbitals
                    fock[i,j] = fock[i,j] + density[k,l]*(eri_matrix[i,k,j,l] - 0.5*eri_matrix[i,k,l,j])
                end
            end
        end
    end
    return fock
end

function normalizeVector(vector, overlap_matrix, total_orbitals)
    magnitude = 0.0
    for i in 1:total_orbitals
        for j in 1:total_orbitals
            magnitude = magnitude + vector[i]*overlap_matrix[i,j]*vector[j]
        end
    end
    return (1.0/sqrt(magnitude))*vector
end

function calTotalEnergy(density, kinetic_matrix, nuclear_matrix, fock_matrix, total_orbitals)
    total_energy = 0.0
    core_matrix = kinetic_matrix + nuclear_matrix
    for i in 1:total_orbitals
        for j in 1:total_orbitals
            total_energy = total_energy + 0.5*density[i,j]*(core_matrix[i,j] + fock_matrix[i,j])
        end
    end
    return total_energy
end

function calTotalEnergy_UHF(density_up, density_down, kinetic_matrix, nuclear_matrix, fock_matrix_up, fock_matrix_down, total_orbitals)
    total_energy = 0.0
    core_matrix = kinetic_matrix + nuclear_matrix
    for i in 1:total_orbitals
        for j in 1:total_orbitals
            total_energy = total_energy + 0.5*density_up[i,j]*(core_matrix[i,j] + fock_matrix_up[i,j]) + 0.5*density_down[i,j]*(core_matrix[i,j] + fock_matrix_down[i,j])
        end
    end
    return total_energy
end

function calNuclearRepulsion(atoms, total_atoms)
    repulsion_energy = 0.0
    for i in 1:total_atoms
        for j in 1:total_atoms
            if (i!=j)
                rij = norm(atoms[i].pos - atoms[j].pos)
                repulsion_energy = repulsion_energy + 0.5*atoms[i].atomic_number*atoms[j].atomic_number/rij
            end
        end
    end
    return repulsion_energy 
end

function calError(vals, old_vals, total_orbitals)
    error = 0.0
    for i in 1:total_orbitals
        error = error + abs(vals[i] - old_vals[i])
    end
    return error
end

function doSCF_UHF(num_iterations, overlap_matrix, kinetic_matrix, nuclear_matrix, eri_matrix, total_orbitals, 
    total_electrons, atoms, total_atoms, tolerance)
    num_occupied_up_orbitals = Int64(ceil(total_electrons/2))
    num_occupied_down_orbitals = Int64(total_electrons - num_occupied_up_orbitals)

    C_up = Matrix{Float64}(I ,total_orbitals, total_orbitals)
    C_down = (1.0/sqrt(total_orbitals))*zeros(Float64, (total_orbitals, total_orbitals))
    error_up = 1.0
    error_down = 1.0
    eigenvalues_up = zeros(Float64, total_orbitals)
    eigenvalues_down = zeros(Float64, total_orbitals)
    println("Starting Run")
    density_up, density_down = constructDensities_UHF(C_up, C_down, total_orbitals, num_occupied_up_orbitals, num_occupied_down_orbitals)
    fock_up, fock_down = buildFockMatrices_UHF(density_up, density_down, kinetic_matrix, nuclear_matrix, eri_matrix, total_orbitals)
    solve_system_up = eigen(fock_up, overlap_matrix)
    solve_system_down = eigen(fock_down, overlap_matrix)
    old_eigenvalues_up = solve_system_up.values
    old_eigenvalues_down = solve_system_down.values
    for k in 1:total_orbitals
        C_up[k,:] = normalizeVector(solve_system_up.vectors[:,k], overlap_matrix, total_orbitals)
        C_down[k,:] = normalizeVector(solve_system_down.vectors[:,k], overlap_matrix, total_orbitals)
    end
    for i in 1:num_iterations
        println("Iteration: ",i)
        density_up, density_down = constructDensities_UHF(C_up, C_down, total_orbitals, num_occupied_up_orbitals, num_occupied_down_orbitals)
        fock_up, fock_down = buildFockMatrices_UHF(density_up, density_down, kinetic_matrix, nuclear_matrix, eri_matrix, total_orbitals)
        solve_system_up = eigen(fock_up, overlap_matrix)
        solve_system_down = eigen(fock_down, overlap_matrix)
        eigenvalues_up = solve_system_up.values
        eigenvalues_down = solve_system_down.values
        for k in 1:total_orbitals
            C_up[k,:] = normalizeVector(solve_system_up.vectors[:,k], overlap_matrix, total_orbitals)
            C_down[k,:] = normalizeVector(solve_system_down.vectors[:,k], overlap_matrix, total_orbitals)
        end
        
        error_up = calError(eigenvalues_up, old_eigenvalues_up, total_orbitals)
        error_down = calError(eigenvalues_down, old_eigenvalues_down, total_orbitals)
        println("Spin Up:", error_up)
        println("Spin Down:", error_down)
        old_eigenvalues_up = eigenvalues_up
        old_eigenvalues_down = eigenvalues_down
        if (abs(error_up) < tolerance && abs(error_down) < tolerance)
            println("Converged")
            break
        end    
    end
    electronic_energy = calTotalEnergy_UHF(density_up, density_down, kinetic_matrix, nuclear_matrix, fock_up, fock_down, total_orbitals)
    nuclear_repulsion = calNuclearRepulsion(atoms, total_atoms)
    total_energy = electronic_energy + nuclear_repulsion
    return C_up, C_down, eigenvalues_up, eigenvalues_down, total_energy
end

function doSCF(num_iterations, overlap_matrix, kinetic_matrix, nuclear_matrix, eri_matrix, total_orbitals, 
    total_electrons, atoms, total_atoms, tolerance, write_status)
    num_occupied_orbitals = Int64(total_electrons/2)
    C = (1.0/sqrt(total_orbitals))*zeros(Float64, (total_orbitals, total_orbitals))
    error = 1.0
    eigenvalues = zeros(Float64, total_orbitals)
    if (write_status == 1)
        println("Starting Run")
    end
    density = constructDensity(C, total_orbitals, num_occupied_orbitals)
    fock = buildFockMatrix(density, kinetic_matrix, nuclear_matrix, eri_matrix, total_orbitals)
    solve_system = eigen(fock, overlap_matrix)
    old_eigenvalues = solve_system.values
    for k in 1:total_orbitals
        C[k,:] = normalizeVector(solve_system.vectors[:,k], overlap_matrix, total_orbitals)
    end

    for i in 1:num_iterations
        if (write_status == 1)
            println("Iteration: ",i)
        end
        density = constructDensity(C, total_orbitals, num_occupied_orbitals)
        fock = buildFockMatrix(density, kinetic_matrix, nuclear_matrix, eri_matrix, total_orbitals)
        solve_system = eigen(fock, overlap_matrix)
        eigenvalues = solve_system.values
        for k in 1:total_orbitals
            C[k,:] = normalizeVector(solve_system.vectors[:,k], overlap_matrix, total_orbitals)
        end
        error = calError(eigenvalues, old_eigenvalues, total_orbitals)
        if (write_status == 1)
            println(error)
        end
        old_eigenvalues = eigenvalues
        if (abs(error) < tolerance)
            if (write_status == 1)
                println("Converged")
            end
            break
        end    
    end
    electronic_energy = calTotalEnergy(density, kinetic_matrix, nuclear_matrix, fock, total_orbitals)
    nuclear_repulsion = calNuclearRepulsion(atoms, total_atoms)
    total_energy = electronic_energy + nuclear_repulsion
    return C, eigenvalues, total_energy
end

end