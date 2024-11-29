using CompChem
using Test

@testset "CompChem.jl" begin
    # Write your tests here.
    @test abs(CompChem.RHF_energy(100, 1e-5, "6-31G", "../basis/6-31G", "../teststruct/water.xyz") + 75.98456188139424) < 1e-3
    @test abs(CompChem.RHF_energy(100, 1e-5, "STO-3G", "../basis/STO-3G", "../teststruct/water.xyz") + 74.96254413336464) < 1e-3
    @test abs(CompChem.RHF_energy(100, 1e-5, "cc-pVDZ", "../basis/cc-pVDZ", "../teststruct/water.xyz") + 76.02617082940125) < 1e-3
end
