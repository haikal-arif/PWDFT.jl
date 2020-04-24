using LinearAlgebra
using Printf
using PWDFT
using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("smearing.jl")
include("create_Ham.jl")
include("ElecVars.jl")
include("test_ElecVars.jl")
include("calc_energies_grad.jl")
include("emin_smearing.jl")
include("linmin_grad.jl")
include("KS_solve_Emin_PCG_Haux_v1.jl")

function print_vec_mat( v::Vector{Matrix{ComplexF64}} )
    Nkspin = length(v)
    for i in 1:Nkspin
        println("Real part of ikspin = ", i)
        display(real(v[i])); println()
        println("Imag part of ikspin = ", i)
        display(imag(v[i])); println()    
    end
end

function main()

    Random.seed!(1234)

    kT = 0.01
    Ham = create_Ham_atom_Al_smearing()
    #Ham = create_Ham_Al_fcc_smearing()
    #Ham = create_Ham_atom_Pt_smearing()
    #Ham = create_Ham_Pt_fcc_smearing()
    #println(Ham)

    #test_ElecVars(Ham)

    evars = ElecVars(Ham)
    KS_solve_Emin_PCG_Haux_v1!( Ham, evars, NiterMax=30 )

    println("Pass here")
end

main()