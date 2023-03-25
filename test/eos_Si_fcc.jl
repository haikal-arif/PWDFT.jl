using PWDFT
using Test
using Dates
using LsqFit
using Printf

function init_Ham_Si_fcc_SCAN(lattice_param::Float64)::Hamiltonian
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0   0.0   0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(lattice_param))

    pspfiles = [joinpath(@__DIR__, "../pseudopotentials", "scan_upf", "Si.SCAN.UPF2")]
    ecutwfc = 25.0
    Hamiltonian(atoms, pspfiles, ecutwfc, meshk=[10, 10, 10], xcfunc="SCAN")
end

"""
This function will generate points for eos calculation
- system_creator : a function to generate hamiltonian, Its signature must be f(Float64) -> Hamiltonian
- lower_bound : fraction of middle_lattice_param that will be its lower bound
- higher_bound : fraction of middle_lattice_param that will be its upper bound
- num_points : Number of points that will be generated
"""
function generate_data_point(
    system_creator::Function,
    lower_bound::Float64=0.8,
    upper_bound::Float64=1.2,
    num_points::Int64=16)::Array{Float64,2}

    num_points = num_points
    calculation_points = range(lower_bound, upper_bound, num_points)
    data_points = zeros(Float64, num_points, 2)
    x_column = 1
    y_column = 2
    Threads.@threads for (index, value) in collect(enumerate(calculation_points))
        data_points[index, x_column] = value
        system::Hamiltonian = system_creator(value)
        println("Starting calculation for thread $index")
        KS_solve_SCF!(system, verbose=false)
        println("Finished calculation for thread $index")
        energies = sum(system.energies)
        data_points[index, y_column] = energies
    end
    curr_day = today()
    println(data_points)
    open(joinpath(@__DIR__, "eos_data_point_Si_fcc_$curr_day.dat"), "w") do file
        println(file, data_points)
    end

    data_points
end

function do_fit_sjeos(volumes, energies)
    # 
    @. model(x, p) = p[1] + p[2] * x^(-1 / 3) + p[3] * x^(-2 / 3) + p[4] * x^(-1)

    p0 = [0.5, 0.5, 0.5, 0.5]

    fit = curve_fit(model, volumes, energies, p0, show_trace=true)

    println("Residuals:")
    println(fit.resid)

    dump(fit)

    p = fit.param
    Ndata = length(volumes)
    for i in 1:Ndata
        V = volumes[i]
        E = energies[i]
        E_m = model(volumes[i], p)
        diffE = abs(E - E_m)  # this is already calculated in residuals
        @printf("%18.10f %18.10f %18.10f %18.10e\n", V, E, E_m, diffE)
    end
    return p
end

function do_fit_parabolic(volume, energies)
    volume * energies
end

@testset "Equation of State Test" begin
    expected_points = 16
    number_of_col = 2
    Silicon_lattice_param = 10.2481
    data_points = zeros(Float64, expected_points, number_of_col)
    @testset "Data points creation" begin
        datas = generate_data_point(init_Ham_Si_fcc_SCAN, Silicon_lattice_param)
        @test size(data_points) == size(datas)
        data_points = datas
    end
    @testset "Equation of state fitting" begin
        lattice_params = data_points[:, 1]
        energies = data_points[:, 2]
        volumes = 0.25 * (BOHR2ANG * lattice_params) .^ 3
        result = do_fit_sjeos(data_points)
    end
end


generate_data_point()
