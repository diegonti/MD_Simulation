module initialization
    ! Module with the initialization of the position and velocity matrices.
    ! Contains two general functions that choose the desired initialization method.
    ! 3 types of initial positions are implemented -> SC, BCC, FCC.
    ! 2 types of initial velocities are implemented -> bimodal, zero.

    use, intrinsic :: iso_fortran_env, only: DP => real64, I64 => int64, input_unit, output_unit
    use mpi
    contains

    subroutine getInitialParams(cell,N,density,M,L,a)
        ! Returns the initial cell parameters of the simulation.
        !
        ! Args: 
        !   cell    (CHARACTER) : String with the name of the cell (sc, bcc, fcc).
        !   N       (INT64)     : Number of particles.
        !   density (REAL64)    : Density of the system.
        !
        ! Returns:
        !   M       (INT64)     : Number of cells in one direction.
        !   L       (REAL64)    : Length of the simulation box.
        !   a       (REAL64)    : Lattice parameter of the unit cell (L/M).

        character(*), intent(in) :: cell
        integer(kind=i64), intent(in) :: N
        double precision, intent(in) :: density
        integer(kind=i64), intent(out) :: M
        double precision, intent(out) :: L, a
        
        L = (real(N,kind=dp)/density)**(1.d0/3.d0)
        if (cell=="sc") then; M = nint((real(N,kind=dp))**(1.d0/3.d0),kind=i64)
        else if (cell=="fcc") then; M = nint((real(N,kind=dp)/4.d0)**(1.d0/3.d0), kind=i64)
        else if (cell=="bcc") then; M = nint((real(N,kind=dp)/2.d0)**(1.d0/3.d0), kind=i64)
        else; print*, "Select a valid initial position: sc, bcc, fcc."
        end if
    
        a = L/real(M, kind=dp)            ! Lattice parameter
        
    end subroutine getInitialParams


    subroutine changeIUnits(lj_epsilon,lj_sigma,mass,density,dt,T)
        ! Changes user Input units to reduced units.
        !
        ! Args:
        !   lj_epsilon  (REAL64) : Lennard Jones epsilon parameter for the gas (in kJ/mol).
        !   lj_sigma    (REAL64) : Lennard Jones sigma parameter for the gas (in Ang).
        !   mass        (REAL64) : Molar Mass of the gas (g/mol).
        !
        ! Inout:
        !   density     (REAL64) : Density of the system (in g/mL).
        !   dt          (REAL64) :  Time-step of the simulation (in ps).
        !   T           (REAL64) : Temperature (in K).

        double precision, intent(in) :: lj_epsilon,lj_sigma,mass
        double precision, intent(inout) :: density,dt,T
        double precision :: ru_time,ru_dens,ru_temp
        double precision, parameter :: Na = 6.0221408d23
        double precision, parameter :: kb = 1.380649d-23


        ! Conversion factors between reduced and real units
        ru_time = sqrt(mass*(lj_sigma*1d-10)**2_i64 / (lj_epsilon*1d6))*1d12    ! time in picoseconds
        ru_dens = 1d24 * mass / (Na*lj_sigma**3_i64)                            ! density in g/mL
        ru_temp = 1d3*lj_epsilon/(kb*Na)                                        ! temperature in K


        T = T / ru_temp
        density = density / ru_dens
        dt = dt / ru_time

    end subroutine changeIUnits

    
    subroutine divide_positions(taskid,numproc,N, sendcounts,displs, imin,imax,local_N)
        ! Divides the whole system by particles and their indices, into chunks of N/numproc atoms. 
        ! It creates also the sendcounts and displs vectors (needed for different MPI subroutines),
        ! which contain the local_N of each process (number of atoms per rank), and the displacements
        ! for each process (the imin for each rank), respectively, from which the imin and imax can be derived. 
        ! Each process will look at the imin_th to the imax_th particle of the full matrix.
        ! 
        ! Args:
        !   taskid    (INT)   : Index of the process (given by MPI_RANK())
        !   numproc   (INT)   : Number of cores (processors) used.
        !   N         (INT64) : Number of particles in the system.
        !
        ! Inout:
        !   sendcounts    (INT[numproc]) : Array containing the number of particles (local_N) of each process (index).
        !   displs        (INT[numproc]) : Array containing the index of the first particle of each process (imin).
        !
        ! Returns:
        !   imin    (INT64) : Index of the first particle of each process.
        !   imax    (INT64) : Index of the last particle of each process.
        !   local_N (INT64) : Number of particles that each process has. 

        implicit none
        integer(kind=i64), intent(in) :: N
        integer, intent(in) :: taskid,numproc
        integer, dimension(:), allocatable,intent(inout):: sendcounts,displs
        integer(kind=i64), intent(out) :: imin,imax,local_N
        integer :: chunk,chunk_extra, k,i


        chunk = int(N)/numproc                  ! Number of particles per task (N/numproc) (int)
        chunk_extra = mod(int(N),numproc)       ! Extra particles if division is not exact 

        ! Fills the sendcounts and displs vectors (needed for different MPI subroutines)
        k = 0
        do I = 1, numproc
            ! The extra particles (if any) will go to the first counts
            if (I > (numproc - chunk_extra)) then
                sendcounts(I) = chunk + 1
            else
                sendcounts(I) = chunk
            end if
            displs(I) = k
            k = k + sendcounts(I)
        end do

        local_N = int(sendcounts(taskid+1),kind=i64)
        imin = int(displs(taskid+1),kind=i64) + 1_i64
        imax = imin + local_N - 1_i64

        ! print*,imin,imax
        ! if (taskid==0) then
        !     print*,displs
        !     print*,sendcounts
        ! end if
        
    end subroutine divide_positions


    subroutine initializePositions(M,a,r,cell,imin,imax,sendcounts,displs)
        ! Chooses and runs the specified initial position.
        ! For the given cell, selects the unit cell matrix and initializes the positions.
        !
        ! Args:
        !   M           (INT64)        : Number of cells in one direction.
        !   a           (REAL64)       : Lattice parameter of the unit cell (L/M)
        !   cell        (CHARACTER)    : String with the name of the cell (sc, bcc, fcc).
        !   imin        (INT64)        : Index of the first particle of each process.
        !   imax        (INT64)        : Index of the last particle of each process.
        !   sendcounts  (INT[numproc]) : Array containing the number of particles (local_N) of each process (index).
        !   displs      (INT[numproc]) : Array containing the index of the first particle of each process (imin).
        !
        ! Inout:
        !   r       (REAL64[3,N]) : 3xN Positions matrix.  

        implicit none
        integer(kind=i64), intent(in) :: M, imin,imax
        double precision, intent(in) :: a
        character(*), intent(in) :: cell
        integer, dimension(:), intent(in):: sendcounts,displs
        double precision, intent(inout), dimension(:,:) :: r
        double precision, allocatable, dimension(:,:) :: ucell


        if ((index(cell,"fcc")==1) .OR. (index(cell,"FCC")==1))  then
            allocate(ucell(3,4))
            ucell = reshape( (/0.d0,0.d0,0.d0, 0.d0,0.5d0,0.5d0,  0.5d0,0.d0,0.5d0,  0.5d0,0.5d0,0.d0/),shape(ucell), order=(/1,2/))
        else if((index(cell,"sc")==1) .OR. (index(cell,"SC")==1)) then
            allocate(ucell(3,1))
            ucell = reshape( (/0.d0,0.d0,0.d0/),shape(ucell), order=(/1,2/))
        else if ((index(cell,"bcc")==1) .OR. (index(cell,"BCC")==1))  then
            allocate(ucell(3,2))
            ucell = reshape( (/0.d0,0.d0,0.d0, 0.5d0,0.5d0,0.5d0/),shape(ucell), order=(/1,2/))
        else
            print*, "Select a valid initial position: sc, bcc, fcc."
        end if

        call initializeGeneral(M,a,ucell,r, imin,imax,sendcounts,displs)     ! Creaties position matrix with the selected unit cell

        deallocate(ucell)

    end subroutine initializePositions


    subroutine initializeVelocities(T,v,init_vel,imin,imax,sendcounts,displs)
        ! Chooses and runs the specified initial velocities.
        !
        ! Args:
        !   T           (REAL64)       : Temperature.
        !   init_vel    (CHARACTER)    : String with the initialization method (bimodal, zero).
        !   imin        (INT64)        : Index of the first particle of each process.
        !   imax        (INT64)        : Index of the last particle of each process.
        !   sendcounts  (INT[numproc]) : Array containing the number of particles (local_N) of each process (index).
        !   displs      (INT[numproc]) : Array containing the index of the first particle of each process (imin).
        !
        ! Inout:
        !   v           (REAL64[3,N]) : 3xN Velocity matrix.

        double precision,intent(in) :: T
        double precision,intent(inout),dimension(:,:) :: v
        character(*), intent(in) :: init_vel
        integer(kind=i64), intent(in) :: imin, imax
        integer, dimension(:), intent(in):: sendcounts,displs 
        integer(kind=i64) :: local_N,d

        local_N = imax - imin + 1_i64

        if ((index(init_vel,"bim")==1) .OR. (index(init_vel,"BIM")==1))  then
            call initBimodal(T,v,imin,imax)
        else if((index(init_vel,"zero")==1) .OR. (index(init_vel,"ZERO")==1)) then
            v(:,imin:imax) = 0d0
        else
            print*, "Select a valid initial velocity: bimodal, zero."
        end if

        do d=1,3
            call MPI_ALLGATHERV(v(d,imin:imax), int(local_N), MPI_DOUBLE_PRECISION, v(d,:), sendcounts, displs, &
                    MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
        end do

    end subroutine initializeVelocities


    subroutine initializeGeneral(M,a,ucell,r, imin,imax,sendcounts,displs)
        ! Initializes the position matrix for a given unit cell.
        ! It does this by replicating the unit cell in all 3 dimensions.
        !
        ! Args:
        !   M           (INT64)        : Number of cells in one direction.
        !   a           (REAL64)       : Lattice parameter of the unit cell (L/M).
        !   ucell       (REAL[:,:])    : Matrix with the unit cell atom's positions.
        !   imin        (INT64)        : Index of the first particle of each process.
        !   imax        (INT64)        : Index of the last particle of each process.
        !   sendcounts  (INT[numproc]) : Array containing the number of particles (local_N) of each process (index).
        !   displs      (INT[numproc]) : Array containing the index of the first particle of each process (imin).
        !
        ! Inout:
        !   r       (REAL64[3,N]) : 3xN Positions matrix.   
            
        implicit none
        integer(kind=i64), intent(in) :: M, imin, imax
        double precision, intent(in) :: a
        double precision, dimension(:,:), intent(in) :: ucell
        integer, dimension(:), intent(in):: sendcounts,displs
        double precision, intent(inout), dimension(:,:) :: r
        integer(kind=i64) :: i,j,k,d, p,at,f, local_N
        integer :: ierror

        r = 0
        local_N = imax - imin + 1_i64

        f = size(ucell, dim=2,kind=i64)     ! Number of whole atoms in a unit cell

        do p=imin,imax
            at = mod(p-1,f) + 1                         ! Current unit cell atom index
            k = mod((p-at)/f,M)                         ! Current Z displacement (in unit cells)
            j = mod(((p-at)/f-k)/M,M)                   ! Current Y displacement (in unit cells)
            i = int((((p-at)/f-k)/M - j)/M,kind=i64)    ! Current X displacement (in unit cells)
            ! The p atom position will be the unit cell one displacez by i,j,k units
            r(:, p) = ucell(:,at) + (/real(i,kind=dp),real(j,kind=dp),real(k,kind=dp)/)
        end do
        
        r(:,imin:imax) = r(:,imin:imax)*a           ! Rescaling the cells with the unit cell latice lenght
        
        ! r(:,imin:imax) = r(:,imin:imax) - L/2     ! Moving particles to have a (0,0,0) centered box

        ! Combine and update all positions matrix for each worker
        do d=1,3
            call MPI_ALLGATHERV(r(d,imin:imax), int(local_N), MPI_DOUBLE_PRECISION, r(d,:), sendcounts, displs, &
                    MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
        end do

    end subroutine initializeGeneral


    subroutine initBimodal(T,v,imin,imax)
        ! Initial velocities as Bimodal distribution.
        !
        ! Args:
        !   T           (REAL64)      : Temperature.
        !   imin        (INT64)       : Index of the first particle of each process.
        !   imax        (INT64)       : Index of the last particle of each process.
        !
        ! Inout:
        !   v           (REAL64[3,N]) : 3xN Velocity matrix.

        implicit none
        double precision,intent(in) :: T
        double precision,intent(inout),dimension(:,:) :: v
        integer(kind=i64), intent(in) :: imin, imax
        double precision :: vi
        double precision :: sign

        if (mod(imin,2_i64)==0) then
            sign = 1.d0
        else
            sign = -1.d0
        end if

        vi = sqrt(T)
        v(:,imin:imax:2) = sign*vi
        v(:,imin+1_i64:imax:2) = -sign*vi

    end subroutine initBimodal


end module initialization