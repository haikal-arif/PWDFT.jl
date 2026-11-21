function test_ortho_01(Ham)

    @assert Ham.need_overlap == false

    for isp in 1:Ham.atoms.Nspecies
        println(Ham.pspots[1])
    end

    ik = 1
    ispin = 1

    Ham.ik = ik
    Ham.ispin = ispin

    Nbasis = Ham.pw.gvecw.Ngw[ik]
    Nstates = Ham.electrons.Nstates
    psi = rand(ComplexF64, Nbasis, Nstates)

    ortho_sqrt!(psi)

    O = psi' * psi
    
    println("O should be an identity matrix\n")
    println("Real O")
    display(real(O))

    println("Imag O")
    display(imag(O))

    ortho_check(psi)

    return
end


function test_ortho_02(Ham)

    for isp in 1:Ham.atoms.Nspecies
        println(Ham.pspots[1])
    end

    ik = 1
    ispin = 1

    Ham.ik = ik
    Ham.ispin = ispin

    Nbasis = Ham.pw.gvecw.Ngw[ik]
    Nstates = Ham.electrons.Nstates
    psi = rand(ComplexF64, Nbasis, Nstates)

    ortho_sqrt!(Ham, psi)

    O = psi' * op_S(Ham, psi)
    
    println("O should be an identity matrix\n")
    println("Real O")
    display(real(O))

    println("Imag O")
    display(imag(O))

    ortho_check(Ham, psi)

    return
end

@testset "orthonormalization without op_S" begin
    for f in [create_Ham_Si_fcc_oncv]
        @test test_ortho_01(f()) == nothing
    end
end


@testset "orthonormalization with op_S" begin
    for f in [create_Ham_Si_fcc_oncv,
              create_Ham_Si_fcc_gbrv]
        @test test_ortho_02(f()) == nothing
    end
    # PAW case needs PSLIB_US_PAW_LDA, which is not included in the repo
    paw_file = joinpath(DIR_PWDFT, "pseudopotentials",
        "PSLIB_US_PAW_LDA", "Si.pz-n-kjpaw_psl.1.0.0.UPF")
    if isfile(paw_file)
        @test test_ortho_02(create_Ham_Si_fcc_paw_pslib()) == nothing
    else
        @test_skip test_ortho_02(create_Ham_Si_fcc_paw_pslib()) == nothing
    end
end
