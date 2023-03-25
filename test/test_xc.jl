rho = [0.1, 0.2, 0.3, 0.4, 0.5]
sigma = [0.2, 0.3, 0.4, 0.5, 0.6]
tau = [0.2, 0.3, 0.4, 0.5, 0.6]

@testset "LDA_VWN xc" begin
    @test calc_epsxc_VWN(LibxcXCCalculator(), rho) ≈ [-0.396206, -0.490557, -0.556226, -0.608272, -0.652089] atol = 1e-5
    @test calc_Vxc_VWN(LibxcXCCalculator(), rho) ≈ [-0.51789, -0.64225, -0.728923, -0.797671, -0.85558] atol = 1e-5
end

# @testset "GGA_PBE xc" begin
#     @test calc_epsxc_PBE(rho, sigma) ≈ [-0.459808, -0.5073, -0.562455, -0.611123, -0.653534] atol=1e-5
# end
@testset "MGGA_SCAN xc" begin
    # From Libxc
    @test XC_x_scan(rho, sigma, tau) ≈ [-0.408368, -0.471540, -0.533396, -0.591370, -0.642808] atol = 1e-5
end

