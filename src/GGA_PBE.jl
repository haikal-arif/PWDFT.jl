using Libxc

# TODO cumulate spin 1,2 as one function or dispatch by julia's type system

function calc_epsxc_PBE( pw::PWGrid, Rhoe::Array{Float64,1} )

    Npoints = size(Rhoe)[1]
    Nspin = 1

    # calculate gRhoe2
    gRhoe = op_nabla( pw, Rhoe )
    gRhoe2 = zeros( Float64, Npoints )
    for ip = 1:Npoints
        gRhoe2[ip] = dot( gRhoe[:,ip], gRhoe[:,ip] )
    end

    eps_x = zeros(Float64,Npoints)
    eps_c = zeros(Float64,Npoints)

    ptr = Libxc.xc_func_alloc()
    # exchange part
    Libxc.xc_func_init(ptr, Libxc.GGA_X_PBE, Nspin)
    Libxc.xc_gga_exc!(ptr, Npoints, Rhoe, gRhoe2, eps_x)
    Libxc.xc_func_end(ptr)

    #
    # correlation part
    Libxc.xc_func_init(ptr, Libxc.GGA_C_PBE, Nspin)
    Libxc.xc_gga_exc!(ptr, Npoints, Rhoe, gRhoe2, eps_c)
    Libxc.xc_func_end(ptr)

    #
    Libxc.xc_func_free(ptr)

    return eps_x + eps_c

end

function calc_epsxc_PBE( pw::PWGrid, Rhoe::Array{Float64,2} )

    Nspin = size(Rhoe)[2]
    if Nspin == 1
        return calc_epsxc_PBE( pw, Rhoe[:,1] )
    end

    Npoints = size(Rhoe)[1]

    # calculate gRhoe2
    gRhoe_up = op_nabla( pw, Rhoe[:,1] )
    gRhoe_dn = op_nabla( pw, Rhoe[:,2] )
    gRhoe2 = zeros( Float64, 3, Npoints )
    for ip = 1:Npoints
        gRhoe2[1,ip] = dot( gRhoe_up[:,ip], gRhoe_up[:,ip] )
        gRhoe2[2,ip] = dot( gRhoe_up[:,ip], gRhoe_dn[:,ip] )
        gRhoe2[3,ip] = dot( gRhoe_dn[:,ip], gRhoe_dn[:,ip] )
    end

    Rhoe_tmp = zeros(2,Npoints)
    Rhoe_tmp[1,:] = Rhoe[:,1]
    Rhoe_tmp[2,:] = Rhoe[:,2]

    eps_x = zeros(Float64,Npoints)
    eps_c = zeros(Float64,Npoints)

    ptr = Libxc.xc_func_alloc()
    # exchange part
    Libxc.xc_func_init(ptr, Libxc.GGA_X_PBE, Nspin)
    Libxc.xc_gga_exc!(ptr, Npoints, Rhoe, gRhoe2, eps_x)
    Libxc.xc_func_end(ptr)

    #
    # correlation part
    Libxc.xc_func_init(ptr, Libxc.GGA_C_PBE, Nspin)
    Libxc.xc_lda_exc!(ptr, Npoints, Rhoe, gRhoe2, eps_c)
    Libxc.xc_func_end(ptr)

    #
    Libxc.xc_func_free(ptr)

    return eps_x + eps_c

end

function calc_Vxc_PBE( pw::PWGrid, Rhoe::Array{Float64,1} )

    Npoints = size(Rhoe)[1]
    Nspin = 1

    # calculate gRhoe2
    gRhoe = op_nabla( pw, Rhoe )
    gRhoe2 = zeros( Float64, Npoints )
    for ip = 1:Npoints
        gRhoe2[ip] = dot( gRhoe[:,ip], gRhoe[:,ip] )
    end

    V_x = zeros(Float64,Npoints)
    V_c = zeros(Float64,Npoints)
    V_xc = zeros(Float64,Npoints)

    Vg_x = zeros(Float64,Npoints)
    Vg_c = zeros(Float64,Npoints)
    Vg_xc = zeros(Float64,Npoints)

    ptr = Libxc.xc_func_alloc()
    # exchange part
    Libxc.xc_func_init(ptr, Libxc.GGA_X_PBE, Nspin)
    Libxc.xc_gga_vxc!(ptr, Npoints, Rhoe, gRhoe2, V_x, Vg_x)
    Libxc.xc_func_end(ptr)

    #
    # correlation part
    Libxc.xc_func_init(ptr, Libxc.GGA_C_PBE, Nspin)
    Libxc.xc_gga_vxc!(ptr, Npoints, Rhoe, gRhoe2, V_c, Vg_c)
    Libxc.xc_func_end(ptr)

    V_xc = V_x + V_c
    Vg_xc = Vg_x + Vg_c

    # gradient correction
    h = zeros(Float64,3,Npoints)
    divh = zeros(Float64,Npoints)

    for ip = 1:Npoints
        h[1,ip] = Vg_xc[ip] * gRhoe[1,ip]
        h[2,ip] = Vg_xc[ip] * gRhoe[2,ip]
        h[3,ip] = Vg_xc[ip] * gRhoe[3,ip]
    end
    # div ( vgrho * gRhoe )
    divh = op_nabla_dot( pw, h )
    #
    for ip = 1:Npoints
        V_xc[ip] = V_xc[ip] - 2.0*divh[ip]
    end


    return V_xc

end

function calc_Vxc_PBE( pw::PWGrid, Rhoe::Array{Float64,2} )

    Nspin = size(Rhoe)[2]
    if Nspin == 1
        return calc_Vxc_PBE( pw, Rhoe[:,1] )
    end

    Npoints = size(Rhoe)[1]

    # calculate gRhoe2
    gRhoe_up = op_nabla( pw, Rhoe[:,1] ) # gRhoe = ∇⋅Rhoe
    gRhoe_dn = op_nabla( pw, Rhoe[:,2] )

    gRhoe2 = zeros( Float64, 3, Npoints )
    for ip = 1:Npoints
        gRhoe2[1,ip] = dot( gRhoe_up[:,ip], gRhoe_up[:,ip] )
        gRhoe2[2,ip] = dot( gRhoe_up[:,ip], gRhoe_dn[:,ip] )
        gRhoe2[3,ip] = dot( gRhoe_dn[:,ip], gRhoe_dn[:,ip] )
    end

    V_xc = zeros(Float64, Npoints, 2)
    V_x  = zeros(Float64, Npoints*2)
    V_c  = zeros(Float64, Npoints*2)

    Vg_xc = zeros(Float64, 3, Npoints)
    Vg_x  = zeros(Float64, 3, Npoints)
    Vg_c  = zeros(Float64, 3, Npoints)

    Rhoe_tmp = zeros(2,Npoints)
    Rhoe_tmp[1,:] = Rhoe[:,1]
    Rhoe_tmp[2,:] = Rhoe[:,2]

    ptr = Libxc.xc_func_alloc()
    # exchange part
    Libxc.xc_func_init(ptr, Libxc.GGA_X_PBE, Nspin)
    Libxc.xc_gga_vxc!(ptr, Npoints, Rhoe, gRhoe2, V_x, Vg_x)
    Libxc.xc_func_end(ptr)

    #
    # correlation part
    Libxc.xc_func_init(ptr, Libxc.GGA_C_PBE, Nspin)
    Libxc.xc_gga_vxc!(ptr, Npoints, Rhoe, gRhoe2, V_c, Vg_c)
    Libxc.xc_func_end(ptr)

    ipp = 0
    for ip = 1:2:2*Npoints
        ipp = ipp + 1
        V_xc[ipp,1] = V_x[ip] + V_c[ip]
        V_xc[ipp,2] = V_x[ip+1] + V_c[ip+1]
    end

    Vg_xc = Vg_x + Vg_c

    h = zeros(Float64,3,Npoints)
    divh = zeros(Float64,Npoints)

    #
    # spin up
    #
    for ip = 1:Npoints
        h[1,ip] = 2*Vg_xc[1,ip]*gRhoe_up[1,ip] + Vg_xc[2,ip]*gRhoe_dn[1,ip]
        h[2,ip] = 2*Vg_xc[1,ip]*gRhoe_up[2,ip] + Vg_xc[2,ip]*gRhoe_dn[2,ip]
        h[3,ip] = 2*Vg_xc[1,ip]*gRhoe_up[3,ip] + Vg_xc[2,ip]*gRhoe_dn[3,ip]
    end

    divh = op_nabla_dot( pw, h )

    # spin up
    for ip = 1:Npoints
        V_xc[ip,1] = V_xc[ip,1] - divh[ip]
    end

    #
    # Spin down
    #
    for ip = 1:Npoints
        h[1,ip] = 2*Vg_xc[3,ip]*gRhoe_dn[1,ip] + Vg_xc[2,ip]*gRhoe_up[1,ip]
        h[2,ip] = 2*Vg_xc[3,ip]*gRhoe_dn[2,ip] + Vg_xc[2,ip]*gRhoe_up[2,ip]
        h[3,ip] = 2*Vg_xc[3,ip]*gRhoe_dn[3,ip] + Vg_xc[2,ip]*gRhoe_up[3,ip]
    end

    #
    divh = op_nabla_dot( pw, h )
    # spin down
    for ip = 1:Npoints
        V_xc[ip,2] = V_xc[ip,2] - divh[ip]
    end

    return V_xc

end
