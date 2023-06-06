using Test
using PWDFT

include("../src/XC_funcs/XC_x_scan.jl")
include("../src/XC_funcs/XC_c_scan.jl")
include("../src/XC_funcs/XC_c_lsda.jl")
include("../src/Libxc_old.jl")
include("../src/XC_funcs/XC_x_slater.jl")

rho = [0.1, 0.2, 0.3, 0.4, 0.5]
sigma = [0.2, 0.3, 0.4, 0.5, 0.6]
tau = [0.2, 0.3, 0.4, 0.5, 0.6]
lapl = [0.0, 0.0, 0.0, 0.0, 0.0]


@testset "LDA_VWN xc" begin
    @test calc_epsxc_VWN(LibxcXCCalculator(), rho) ≈ [-0.396206, -0.490557, -0.556226, -0.608272, -0.652089] atol = 1e-5
    @test calc_Vxc_VWN(LibxcXCCalculator(), rho) ≈ [-0.51789, -0.64225, -0.728923, -0.797671, -0.85558] atol = 1e-5
end

# @testset "GGA_PBE xc" begin
#     @test calc_epsxc_PBE(rho, sigma) ≈ [-0.459808, -0.5073, -0.562455, -0.611123, -0.653534] atol = 1e-5
# end

@testset "MGGA_SCAN xc" begin
    exc_libxc = zeros(5)
    ex_libxc = zeros(5)
    ec_libxc = zeros(5)
    xc_func_ptr = Libxc_xc_func_alloc()
    Libxc_xc_func_init(xc_func_ptr, 263, 1)
    Libxc_xc_mgga_exc!(xc_func_ptr, 5, rho, sigma, lapl, tau, ex_libxc)
    Libxc_xc_func_init(xc_func_ptr, 267, 1)
    Libxc_xc_mgga_exc!(xc_func_ptr, 5, rho, sigma, lapl, tau, ec_libxc)
    Libxc_xc_func_end(xc_func_ptr)

    exc_libxc = ex_libxc + ec_libxc

    ex, v1, v2, v3 = XC_x_scan(rho, sigma, tau)
    @test ex ≈ ex_libxc atol = 1e-5
    ec, v1, v2, v3 = XC_c_scan(rho, sigma, tau)
    @test ex + ec ≈ exc_libxc atol = 1e-5
end

