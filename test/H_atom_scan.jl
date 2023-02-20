using PWDFT
using Test
using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")

function init_Ham_H_atom_SCAN()

    atoms = Atoms(xyz_string="""
   1

   H  0.0  0.0  0.0
   """, LatVecs=gen_lattice_sc(16.0))

    pspfiles = [joinpath(DIR_PWDFT, "pseudopotentials", "scan_upf", "H.SCAN.UPF")]
    ecutwfc = 20.0

    return Hamiltonian(atoms, pspfiles, ecutwfc, xcfunc="SCAN", meshk=[10, 10, 10])
end

const Etot_H_atom_PBE = -0.4528600294 # ABINIT


@testset "H atom SCF Rhoe mix simple" begin
    Random.seed!(1234)
    Ham = init_Ham_H_atom_SCAN()
    KS_solve_SCF!(Ham, betamix=0.7, verbose=true)
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_H_atom_PBE atol = 1e-5
end

@testset "H atom SCF Rhoe mix linear adaptive" begin
    Random.seed!(1234)
    Ham = init_Ham_H_atom_PBE()
    KS_solve_SCF!(Ham, betamix=0.1, mix_method="linear_adaptive", verbose=true)
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_H_atom_PBE atol = 1e-5
end

@testset "H atom SCF Rhoe mix anderson" begin
    Random.seed!(1234)
    Ham = init_Ham_H_atom_PBE()
    KS_solve_SCF!(Ham, betamix=0.7, mix_method="anderson", verbose=true)
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_H_atom_PBE atol = 1e-5
end

@testset "H atom SCF Rhoe mix pulay" begin
    Random.seed!(1234)
    Ham = init_Ham_H_atom_PBE()
    KS_solve_SCF!(Ham, betamix=0.5, mix_method="pulay", verbose=true)
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_H_atom_PBE atol = 1e-5
end

@testset "H atom SCF Rhoe mix ppulay" begin
    Random.seed!(1234)
    Ham = init_Ham_H_atom_PBE()
    KS_solve_SCF!(Ham, betamix=0.5, mix_method="ppulay", verbose=true)
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_H_atom_PBE atol = 1e-5
end

@testset "H atom SCF Rhoe mix rpulay" begin
    Random.seed!(1234)
    Ham = init_Ham_H_atom_PBE()
    KS_solve_SCF!(Ham, betamix=0.5, mix_method="rpulay", verbose=true)
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_H_atom_PBE atol = 1e-5
end

@testset "H atom SCF Rhoe mix broyden" begin
    Random.seed!(1234)
    Ham = init_Ham_H_atom_PBE()
    KS_solve_SCF!(Ham, betamix=0.5, mix_method="broyden", verbose=true)
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_H_atom_PBE atol = 1e-5
end


@testset "H atom SCF potential mix simple" begin
    Random.seed!(1234)
    Ham = init_Ham_H_atom_PBE()
    KS_solve_SCF_potmix!(Ham, betamix=0.5, mix_method="simple", verbose=true)
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_H_atom_PBE atol = 1e-5
end

@testset "H atom SCF potential mix linear_adaptive" begin
    Random.seed!(1234)
    Ham = init_Ham_H_atom_PBE()
    KS_solve_SCF_potmix!(Ham, betamix=0.1, mix_method="linear_adaptive", verbose=true)
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_H_atom_PBE atol = 1e-5
end

@testset "H atom SCF potential mix broyden" begin
    Random.seed!(1234)
    Ham = init_Ham_H_atom_PBE()
    KS_solve_SCF_potmix!(Ham, betamix=0.5, mix_method="broyden", verbose=true)
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_H_atom_PBE atol = 1e-5
end


@testset "H atom Emin PCG Fletcher-Reeves" begin
    Random.seed!(1234)
    Ham = init_Ham_H_atom_PBE()
    KS_solve_Emin_PCG!(Ham, i_cg_beta=1, verbose=true)
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_H_atom_PBE atol = 1e-5
end

@testset "H atom Emin PCG Polak-Ribiere" begin
    Random.seed!(1234)
    Ham = init_Ham_H_atom_PBE()
    KS_solve_Emin_PCG!(Ham, i_cg_beta=2, verbose=true)
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_H_atom_PBE atol = 1e-5
end

@testset "H atom Emin PCG Hestenes-Stiefeld" begin
    Random.seed!(1234)
    Ham = init_Ham_H_atom_PBE()
    KS_solve_Emin_PCG!(Ham, i_cg_beta=3, verbose=true)
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_H_atom_PBE atol = 1e-5
end

@testset "H atom Emin PCG Dai-Yuan" begin
    Random.seed!(1234)
    Ham = init_Ham_H_atom_PBE()
    KS_solve_Emin_PCG!(Ham, i_cg_beta=4, verbose=true)
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_H_atom_PBE atol = 1e-5
end
