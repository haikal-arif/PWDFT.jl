import Base: show
function show(io::IO, Ham::Hamiltonian; header=true)
    if header
        @printf("\n")
        @printf("                                  -----------\n")
        @printf("                                  Hamiltonian\n")
        @printf("                                  -----------\n")
        @printf("\n")
    end
    @printf(io, "size (MiB) = %18.5f\n", Base.summarysize(Ham) / 1024 / 1024)
    println(io, "")
    println(io, "xcfunc       = ", Ham.xcfunc)
    
    if typeof(Ham.xc_calc) == LibxcXCCalculator
        println(io, "xc_calc.x_id = ", Ham.xc_calc.x_id)
        println(io, "xc_calc.c_id = ", Ham.xc_calc.c_id)
    elseif typeof(Ham.xc_calc) == XCCalculator
        println("Using internal implementation")
    end


    println(io, "")
    show(io, Ham.atoms)
    show(io, Ham.pw)
    show(io, Ham.pw.gvecw.kpoints)
    show(io, Ham.electrons)
    for isp = 1:Ham.atoms.Nspecies
        show(io, Ham.pspots[isp])
    end
end
show(Ham::Hamiltonian; header=true) = show(stdout, Ham, header=header)