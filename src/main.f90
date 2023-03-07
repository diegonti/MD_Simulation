program main
    use, intrinsic :: iso_fortran_env, only: DP => real64, I64 => int64, i32 => int32, input_unit, output_unit, error_unit
    ! Module definitions
    use            :: initialization, only: getInitialParams, initializePositions, initializeVelocities, initializeGeneral, &
    initBimodal
    use            :: testing
    use            :: readers_m,      only: read_nml
    implicit none

    ! ~ Memory definition ~
    ! Array variables
    double precision, allocatable,dimension(:,:) :: r, v

    ! Scalar variables
    integer(kind=i32) :: status_cli
    integer(kind=i64) :: M, N, n_steps, write_file, write_stats, gdr_num_bins, write_frame, log_unit
    real(kind=dp)     :: init_time, end_time, density, L, a, T, lj_epsilon, lj_sigma, mass, dt, &
    andersen_nu

    ! String variables
    character(len=2048)           :: nml_path, sim_name, log_name
    character(len=:), allocatable :: cell_type
    character(len=:), allocatable :: init_vel


    !!! ~ MAIN PROGRAM ~ !!!

    call cpu_time(init_time)
    write(output_unit, '(A)') '~ Welcome to MDEMI!! ~'
    write(output_unit, '(A)') 'Reading input file...'

    call get_command_argument(1, nml_path, status=status_cli)
    if (status_cli /= 0_i32) then
        write(error_unit, '(A)') "Cannot understand the input nml file. Please, provide-it as the first argument"
        stop 1
    end if

    call read_nml(param_file=nml_path, lj_epsilon=lj_epsilon, lj_sigma=lj_sigma, mass=mass, timestep=dt, &
    density=density, andersen_nu=andersen_nu, n_particles=N, n_steps=n_steps, write_file=write_file, &
    write_stats=write_stats, gdr_num_bins=gdr_num_bins, write_frame=write_frame, sim_name=sim_name, &
    cell_type=cell_type, init_velocities=init_vel, temperature=T)

    ! Trimming character variables in order to avoid blank spaces
    sim_name = trim(sim_name)
    cell_type = trim(cell_type)
    init_vel = trim(init_vel)

    write(output_unit, '(A)') 'Successfully loaded parameter file, starting simulation'

    ! ~ Memmory allocation ~
    allocate(r(3,N))
    allocate(v(3,N))

    ! Opening files
    ! log_unit -> file where the simulation time, energy, instant temperature, etc.. will be placed

    log_name = trim(sim_name) // "_logfile.log"
    open(newunit=log_unit, file=trim(log_name), access='sequential', action='write', &
    status='replace', form='formatted')

    ! ~ Initialization of the system ~
    call getInitialParams(cell_type,N,density,M,L,a)
    call initializePositions(M,a,r,cell_type)
    call initializeVelocities(T,v,init_vel)

    !call testMatrix(r)


    ! ~ Closing files ~
    close(log_unit)
    
    ! ~ Program finalization ~
    call cpu_time(end_time)
    write(output_unit, '(A,F12.8,A)') 'Execution done in: ', end_time - init_time, ' seconds.'

end program main