# Apply kinetic operator to wave function in reciprocal space

function op_K( Ham::Hamiltonian, psi::Array{ComplexF64,2} )
    #
    ik = Ham.ik

    Nstates = size(psi)[2]

    pw = Ham.pw
    Ngwx = pw.gvecw.Ngwx
    Ngw = pw.gvecw.Ngw
    idx_gw2g = pw.gvecw.idx_gw2g[ik]
    k = pw.gvecw.kpoints.k[:,ik]

    Gw = zeros(3)
    out = Array{ComplexF64}(undef,size(psi))

    for ist = 1:Nstates
        for igk = 1:Ngw[ik]
            ig = idx_gw2g[igk]
            Gw[:] = pw.gvec.G[:,ig] + k[:]
            Gw2 = Gw[1]^2 + Gw[2]^2 + Gw[3]^2
            out[igk,ist] = psi[igk,ist]*Gw2
        end
    end

    return 0.5*out # two minus signs -> positive
end

# This function is used by CheFSI
function op_K( Ham::Hamiltonian, psi::Array{ComplexF64,1} )
    #
    ik = Ham.ik

    pw = Ham.pw
    Ngwx = pw.gvecw.Ngwx
    Ngw = pw.gvecw.Ngw
    idx_gw2g = pw.gvecw.idx_gw2g[ik]
    k = pw.gvecw.kpoints.k[:,ik]

    Gw = zeros(3)
    out = Array{ComplexF64}(undef,size(psi))

    for igk = 1:Ngw[ik]
        ig = idx_gw2g[igk]
        Gw[:] = pw.gvec.G[:,ig] + k[:]
        Gw2 = Gw[1]^2 + Gw[2]^2 + Gw[3]^2
        out[igk] = psi[igk]*Gw2
    end

    return 0.5*out # two minus signs -> positive
end

