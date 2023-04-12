program main
    use, intrinsic :: iso_fortran_env, only: DP => real64, I64 => int64, i32 => int32, input_unit, output_unit, error_unit
    ! Module definitions
    use            :: mpi
    use            :: initialization, only: changeIUnits, getInitialParams, initializePositions, initializeVelocities, &
                                            divide_positions
    use            :: testing
    use            :: readers_m,      only: read_nml
    use            :: integrators,    only: mainLoop
    implicit none

    ! ~ Memory definition ~
    ! Array variables
    double precision, allocatable,dimension(:,:) :: r, v
    integer, allocatable, dimension(:)           :: state

    ! Scalar variables
    integer(kind=i32) :: status_cli
    integer(kind=i64) :: M, N, n_steps, write_file, write_stats, gdr_num_bins, n_sweeps, &
                         write_frame, log_unit, traj_unit, rdf_unit, msd_unit
    real(kind=dp)     :: init_time, end_time, density, L, a, T, lj_epsilon, lj_sigma, mass, dt, &
                         andersen_nu, cutoff, vcutoff
    integer, parameter:: seed_number = 165432156
    integer           :: state_size

    ! String variables
    character(len=2048) :: nml_path, sim_name, log_name, traj_name, rdf_name, msd_name
    character(len=2048) :: cell_type, init_vel

    ! MPI memory definition
    integer, parameter                  :: MASTER = 0
    integer                             :: taskid, ierror, numproc
    integer(kind=i64)                   :: imin, imax, local_N
    integer,  allocatable, dimension(:) :: sendcounts, displs


    !!! ~ MAIN PROGRAM ~ !!!

    ! MPI initialization
    call MPI_INIT(ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numproc, ierror)

    ! Setting random seed of RNG
    call random_seed( size=state_size )
    allocate(state(state_size))
    state = seed_number
    call random_seed( put=state )

    if (taskid == MASTER) then

        init_time = MPI_Wtime()

        write(output_unit, '(A)') '~ Welcome to MDEMI!! ~'
        write(output_unit, '(A)') 'Reading input file...'

        call get_command_argument(1, nml_path, status=status_cli)
        if (status_cli /= 0_i32) then
            write(error_unit, '(A)') "Cannot understand the input nml file. Please, provide-it as the first argument"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierror)
        end if

        call read_nml(param_file=nml_path, lj_epsilon=lj_epsilon, lj_sigma=lj_sigma, mass=mass, timestep=dt, &
        density=density, andersen_nu=andersen_nu, n_particles=N, n_steps=n_steps, write_file=write_file, &
        write_stats=write_stats, gdr_num_bins=gdr_num_bins, n_sweeps=n_sweeps, write_frame=write_frame,  &
        sim_name=sim_name, cell_type=cell_type, init_velocities=init_vel, temperature=T, cutoff=cutoff, &
        vlcutoff=vcutoff)
        
        write(output_unit, '(A)') 'Successfully loaded parameter file, starting simulation'

        ! Opening files
        ! log_unit  -> file where the simulation time, energy, instant temperature, etc.. will be placed
        ! traj_name -> trajectory file, where the xyz of each snapshot is placed 
        ! rdf_name -> RDF file, where the RDF will be written 
        ! msd_name -> MSD file, where the MSD will be written 

        log_name = trim(sim_name) // "_logfile.log"
        open(newunit=log_unit, file=trim(log_name), access='sequential', action='write', &
        status='replace', form='formatted')

        traj_name = trim(sim_name) // "_trajectory.xyz"
        open(newunit=traj_unit, file=trim(traj_name), access='sequential', action='write', &
        status='replace', form='formatted')

        rdf_name = trim(sim_name) // "_rdf.log"
        open(newunit=rdf_unit, file=trim(rdf_name), access='sequential', action='write', &
        status='replace', form='formatted')
        
        msd_name = trim(sim_name) // "_msd.log"
        open(newunit=msd_unit, file=trim(msd_name), access='sequential', action='write', &
        status='replace', form='formatted')

    end if

    ! INT64 Broadcasting
    call MPI_Bcast(N,            1, MPI_Int, MASTER, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(n_steps,      1, MPI_Int, MASTER, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(write_file,   1, MPI_Int, MASTER, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(write_stats,  1, MPI_Int, MASTER, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(write_frame,  1, MPI_Int, MASTER, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(n_sweeps,     1, MPI_Int, MASTER, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(gdr_num_bins, 1, MPI_Long_Int, MASTER, MPI_COMM_WORLD, ierror)

    ! REAL64 Broadcasting
    call MPI_Bcast(lj_epsilon,  1, MPI_DOUBLE_PRECISION, MASTER, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(lj_sigma,    1, MPI_DOUBLE_PRECISION, MASTER, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(mass,        1, MPI_DOUBLE_PRECISION, MASTER, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(dt,          1, MPI_DOUBLE_PRECISION, MASTER, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(density,     1, MPI_DOUBLE_PRECISION, MASTER, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(andersen_nu, 1, MPI_DOUBLE_PRECISION, MASTER, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(T,           1, MPI_DOUBLE_PRECISION, MASTER, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(cutoff,      1, MPI_DOUBLE_PRECISION, MASTER, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(vcutoff,     1, MPI_DOUBLE_PRECISION, MASTER, MPI_COMM_WORLD, ierror)
    
    ! CHARACTER broadcasting
    call MPI_Bcast(cell_type, 2048, MPI_CHARACTER, MASTER, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(init_vel, 2048, MPI_CHARACTER, MASTER, MPI_COMM_WORLD, ierror)

    ! ~ Memmory allocation ~
    allocate(r(3,N))
    allocate(v(3,N))
    allocate(sendcounts(numproc))
    allocate(displs(numproc))

    ! ~ Initialization of the system ~
    call changeIUnits(lj_epsilon,lj_sigma,mass,density,dt,T)
    call getInitialParams(trim(cell_type),N,density,M,L,a)

    call divide_positions(taskid,numproc,N, sendcounts,displs,imin,imax,local_N)
    call initializePositions(M,a,r,trim(cell_type),imin,imax,sendcounts,displs)
    call initializeVelocities(T,v,init_vel,imin,imax,sendcounts,displs)


    ! ~ Starting the trajectory of the system ~

    call mainLoop(log_unit, traj_unit, rdf_unit, msd_unit, lj_epsilon, lj_sigma, mass, &
    n_steps, dt, L, T, andersen_nu, cutoff*L, gdr_num_bins, n_sweeps, r, v, write_stats, &
    write_frame, taskid, imin, imax, sendcounts, displs, local_N, vcutoff*L)    


    ! ~ Memmory deallocation ~
    deallocate(r)
    deallocate(v)
    deallocate(sendcounts)
    deallocate(displs)

    if (taskid == MASTER) then

        ! ~ Closing files ~ 
        close(log_unit)
        close(traj_unit)
        close(rdf_unit)
        close(msd_unit)

        ! ~ Program finalization ~
        end_time = MPI_Wtime()
        write(output_unit,'(A)') ''
        write(output_unit, '(A,ES18.8e4,A)') 'Execution done in: ', end_time - init_time, ' seconds.'
    
    end if

    call MPI_FINALIZE(ierror)
    
end program main
