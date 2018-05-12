"""
GVectors for density and potentials
"""
struct GVectors
    Ng::Int64
    G::Array{Float64,2}
    G2::Array{Float64,1}
    idx_g2r::Array{Int64,1}
end


"""
GVectors for wave function
"""
struct GVectorsW
    Ngwx::Int64          # maximum(Ngk)
    Ngw::Array{Int64,1}  # no of GvectorsW for each kpoints
    idx_gw2g::Array{Array{Int64,1},1}
    idx_gw2r::Array{Array{Int64,1},1}
    kpoints::KPoints
end

struct PWGrid
    ecutwfc::Float64
    ecutrho::Float64
    Ns::Tuple{Int64,Int64,Int64}
    LatVecs::Array{Float64,2}
    RecVecs::Array{Float64,2}
    Ω::Float64
    r::Array{Float64,2}
    gvec::GVectors
    gvecw::GVectorsW
    planfw::Base.DFT.FFTW.cFFTWPlan{Complex{Float64},-1,false,3}
    planbw::Base.DFT.ScaledPlan{Complex{Float64},Base.DFT.FFTW.cFFTWPlan{Complex{Float64},1,false,3},Float64}
end

"""
LatVecs defines three lattice vectors: v1, v2, and v3 arranged by column.
"""
function PWGrid( ecutwfc::Float64, LatVecs::Array{Float64,2}; kpoints=nothing )

    ecutrho = 4.0*ecutwfc
    #
    RecVecs = 2*pi*inv(LatVecs')
    Ω = det(LatVecs)
    #
    LatVecsLen = Array{Float64}(3)
    LatVecsLen[1] = norm(LatVecs[:,1])
    LatVecsLen[2] = norm(LatVecs[:,2])
    LatVecsLen[3] = norm(LatVecs[:,3])

    Ns1 = 2*round( Int64, sqrt(ecutrho/2)*LatVecsLen[1]/pi ) + 1
    Ns2 = 2*round( Int64, sqrt(ecutrho/2)*LatVecsLen[2]/pi ) + 1
    Ns3 = 2*round( Int64, sqrt(ecutrho/2)*LatVecsLen[3]/pi ) + 1

    Ns1 = good_fft_order(Ns1)
    Ns2 = good_fft_order(Ns2)
    Ns3 = good_fft_order(Ns3)

    Ns = (Ns1,Ns2,Ns3)

    Npoints = prod(Ns)
    r = init_grid_R( Ns, LatVecs )
    
    gvec = init_gvec( Ns, RecVecs, ecutrho )

    if kpoints==nothing
        kpoints = KPoints( 1, zeros(3,1), [1.0], RecVecs )
    end

    gvecw = init_gvecw( ecutwfc, gvec, kpoints )

    planfw = plan_fft( zeros(Ns) )
    planbw = plan_ifft( zeros(Ns) )

    return PWGrid( ecutwfc, ecutrho, Ns, LatVecs, RecVecs, Ω, r, gvec, gvecw,
                   planfw, planbw )
end

function mm_to_nn(mm::Int64,S::Int64)
    if mm > S/2
        return mm - S
    else
        return mm
    end
end


function calc_Ng( Ns, RecVecs, ecutrho )
    ig = 0
    Ng = 0
    #
    G = zeros(Float64,3)
    #
    for k in 0:Ns[3]-1
    for j in 0:Ns[2]-1
    for i in 0:Ns[1]-1
        ig = ig + 1
        gi = mm_to_nn( i, Ns[1] )
        gj = mm_to_nn( j, Ns[2] )
        gk = mm_to_nn( k, Ns[3] )
        G[1] = RecVecs[1,1]*gi + RecVecs[1,2]*gj + RecVecs[1,3]*gk
        G[2] = RecVecs[2,1]*gi + RecVecs[2,2]*gj + RecVecs[2,3]*gk
        G[3] = RecVecs[3,1]*gi + RecVecs[3,2]*gj + RecVecs[3,3]*gk
        G2 = G[1]^2 + G[2]^2 + G[3]^2
        if 0.5*G2 < ecutrho
            Ng = Ng + 1
        end
    end
    end
    end
    return Ng
end

function init_gvec( Ns, RecVecs, ecutrho )

    Ng = calc_Ng( Ns, RecVecs, ecutrho )

    G_temp = zeros(Float64,3)

    G  = Array{Float64}(3,Ng)
    G2 = Array{Float64}(Ng)
    idx_g2r = Array{Int64}(Ng)

    ig = 0
    ip = 0
    for k in 0:Ns[3]-1
    for j in 0:Ns[2]-1
    for i in 0:Ns[1]-1
        ip = ip + 1
        gi = mm_to_nn( i, Ns[1] )
        gj = mm_to_nn( j, Ns[2] )
        gk = mm_to_nn( k, Ns[3] )
        G_temp[1] = RecVecs[1,1]*gi + RecVecs[1,2]*gj + RecVecs[1,3]*gk
        G_temp[2] = RecVecs[2,1]*gi + RecVecs[2,2]*gj + RecVecs[2,3]*gk
        G_temp[3] = RecVecs[3,1]*gi + RecVecs[3,2]*gj + RecVecs[3,3]*gk
        G2_temp = G_temp[1]^2 + G_temp[2]^2 + G_temp[3]^2
        if 0.5*G2_temp < ecutrho
            ig = ig + 1
            @inbounds G[:,ig] = G_temp[:]
            @inbounds G2[ig] = G2_temp
            @inbounds idx_g2r[ig] = ip
        end
    end
    end
    end

    return GVectors( Ng, G, G2, idx_g2r )
end


function init_gvecw( ecutwfc, gvec::GVectors, kpoints::KPoints )
    G = gvec.G
    Ng = gvec.Ng
    idx_g2r = gvec.idx_g2r
    #
    kpts = kpoints.k
    Nkpt = kpoints.Nkpt
    #
    Gk2 = zeros(Float64,Ng)
    Gk = zeros(Float64,3)
    idx_gw2g = Array{Array{Int64,1},1}(Nkpt)
    idx_gw2r = Array{Array{Int64,1},1}(Nkpt)
    Ngw = Array{Int64,1}(Nkpt)
    #
    for ik = 1:Nkpt
        for ig = 1:Ng
            Gk[:] = G[:,ig] .+ kpts[:,ik]
            Gk2[ig] = Gk[1]^2 + Gk[2]^2 + Gk[3]^2
        end
        idx_gw2g[ik] = findn( 0.5*Gk2 .< ecutwfc )
        idx_gw2r[ik] = idx_g2r[idx_gw2g[ik]]
        Ngw[ik] = length(idx_gw2g[ik])
        #@printf("ik = %8d, k = [%10.5f,%10.5f,%10.5f]\n",
        #        ik, kpts[1,ik], kpts[2,ik], kpts[3,ik])
        #@printf("Ngw = %8d\n\n", Ngw[ik])
    end

    # print out memory information
    #Nelems = sum(Ngw)
    #memMiB = 2*Nelems*sizeof(Int64)/1024.0/1024.0
    #@printf("mem = %f MiB\n\n", memMiB)
    
    Ngwx = maximum(Ngw)

    return GVectorsW( Ngwx, Ngw, idx_gw2g, idx_gw2r, kpoints )

end


function init_grid_R( Ns, LatVecs )
    #
    Npoints = prod(Ns)
    #
    R = Array{Float64}(3,Npoints)
    ip = 0
    for k in 0:Ns[3]-1
    for j in 0:Ns[2]-1
    for i in 0:Ns[1]-1
        ip = ip + 1
        @inbounds R[1,ip] = LatVecs[1,1]*i/Ns[1] + LatVecs[1,2]*j/Ns[2] + LatVecs[1,3]*k/Ns[3]
        @inbounds R[2,ip] = LatVecs[2,1]*i/Ns[1] + LatVecs[2,2]*j/Ns[2] + LatVecs[2,3]*k/Ns[3]
        @inbounds R[3,ip] = LatVecs[3,1]*i/Ns[1] + LatVecs[3,2]*j/Ns[2] + LatVecs[3,3]*k/Ns[3]
    end
    end
    end
    #
    return R
end


# Overloaded println

import Base.println

function println( pw::PWGrid )
    @printf("\n")
    @printf("                                     ------\n")
    @printf("                                     PWGrid\n")
    @printf("                                     ------\n")
    LatVecs = pw.LatVecs
    RecVecs = pw.RecVecs
    @printf("\n")
    @printf("Direct lattice vectors:\n")
    @printf("\n")
    for i = 1:3
        @printf("%18.10f %18.10f %18.10f\n", LatVecs[i,1], LatVecs[i,2], LatVecs[i,3])
    end
    @printf("\n")
    @printf("Reciprocal lattice vectors:\n")
    @printf("\n")
    for i = 1:3
        @printf("%18.10f %18.10f %18.10f\n", RecVecs[i,1], RecVecs[i,2], RecVecs[i,3])
    end
    @printf("\n")
    @printf("Direct lattive volume = %18.10f bohr^3\n", pw.Ω )
    @printf("ecutwfc               = %18.10f Ha\n", pw.ecutwfc)
    @printf("ecutrho               = %18.10f Ha\n", pw.ecutrho)    
    @printf("Sampling points       = (%5d,%5d,%5d)\n", pw.Ns[1], pw.Ns[2], pw.Ns[3])
    #
    println( pw.gvec )
    println( pw.gvec, pw.gvecw )
end

function println( gvec::GVectors )
    Ng = gvec.Ng
    G = gvec.G
    G2 = gvec.G2
    
    @printf("\n")
    @printf("                                    --------\n")
    @printf("                                    GVectors\n")
    @printf("                                    --------\n")
    @printf("\n")
    @printf("Ng = %12d\n", Ng)
    @printf("\n")
    for ig = 1:3
        @printf("%8d [%18.10f,%18.10f,%18.10f] : %18.10f\n", ig, G[1,ig], G[2,ig], G[3,ig], G2[ig])        
    end
    @printf(" ....... \n")
    for ig = Ng-3:Ng
        @printf("%8d [%18.10f.%18.10f,%18.10f] : %18.10f\n", ig, G[1,ig], G[2,ig], G[3,ig], G2[ig])
    end
    @printf("\n")
    @printf("Max G2 = %18.10f\n", maximum(G2))
end

function println( gvec::GVectors, gvecw::GVectorsW )
    G = gvec.G
    G2 = gvec.G2

    Ngwx = gvecw.Ngwx
    Ngw = gvecw.Ngw
    k = gvecw.kpoints.k
    Nkpt = gvecw.kpoints.Nkpt

    @printf("\n")
    @printf("                                    ---------\n")
    @printf("                                    GVectorsW\n")
    @printf("                                    ---------\n")
    @printf("\n")
    @printf("Ngwx = %12d\n", Ngwx)
    
    #ik = 1  # sample for first kpoint
    
    for ik = 1:Nkpt

        idx_gw2g = gvecw.idx_gw2g[ik]

        Gw = zeros(3,Ngw[ik])
        Gw2 = zeros(Ngw[ik])
        for igk = 1:Ngw[ik]
            ig = idx_gw2g[igk]
            Gw[:,igk] = G[:,ig] + k[:,ik]
            Gw2[igk] = Gw[1,igk]^2 + Gw[2,igk]^2 + Gw[3,igk]^2
        end

        @printf("\n")
        @printf("Several GvectorsW for k = [%f,%f,%f]\n", k[1,ik], k[2,ik], k[3,ik])
        @printf("\n")
        for igk = 1:3
            @printf("%8d [%18.10f,%18.10f,%18.10f] : %18.10f\n",
                     igk, Gw[1,igk], Gw[2,igk], Gw[3,igk], Gw2[igk])
        end
        @printf(" ....... \n")
        for igk = Ngw[ik]-3:Ngw[ik]
            @printf("%8d [%18.10f.%18.10f,%18.10f] : %18.10f\n",
                    igk, Gw[1,igk], Gw[2,igk], Gw[3,igk], Gw2[igk])
        end
        @printf("\n")
        @printf("Max Gw2[%d] = %18.10f\n", ik, maximum(Gw2))    
    end
end

