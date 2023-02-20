using PWDFT
using Test
using Dates
using Threads

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

function generate_data_point()::Array{Float64, 2}
    primary_guess = 10.2481
    lower_bound = primary_guess * 0.8
    upper_bound = primary_guess * 1.2
    num_points = 15
    calculation_points = range(lower_bound, upper_bound, num_points)
    data_points = zeros(Float64, num_points, 2)
    x_column = 1
    y_column = 2
    Threads.@threads for (index, value) in enumerate(calculation_points)
        data_points[index, x_column] = value
        system::Hamiltonian = init_Ham_Si_fcc_SCAN(value)
        KS_solve_SCF!(system)
        energies = sum(system.energies)
        data_points[index, y_column] = energies
    end
    curr_time = now()
    open(joinpath(@__DIR__, "eos_data_point_Si_fcc_$curr_time.dat"), "w") do file
        println(file, data_points)
    end


end

generate_data_point()