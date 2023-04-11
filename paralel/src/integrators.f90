module integrators
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64, output_unit, i32 => INT32
    use :: mpi
    use :: periodic_bc, only: PBC
    use :: potential_m, only: calc_KE, calc_pressure, calc_vdw_force, calc_vdw_pbc, calc_Tinst,&
           compute_com_momenta, compute_vlist, update_vlist
    use :: simulation,  only: MSD, g_r
    use :: writers_m,   only: writePositions, writeRdf, writeMSD, writeSystem
    use :: testing

    implicit none
    public :: verlet_step, vv_integrator,vel_andersen
    private :: boxmuller

contains

    subroutine mainLoop(log_unit,traj_unit,rdf_unit,msd_unit,lj_epsilon,lj_sigma,mass,N_steps,dt,L,T,nu,cutoff,&
        gdr_num_bins,n_sweeps,r,v,write_log, write_pos, irank, imin, imax, sendcounts, displs, local_N, vcutoff)
        ! Main Simulation Loop
        !
        ! Args:
        !   *_unit      (INT64)     : Integer unit of the lof, trajectory and rdf files. 
        !   lj_*        (REAL64)    : Lennard Jones parameters, epsilon (kJ/mol) and sigma (Ang).
        !   mass        (REAL64)    : Molar mass (g/mol) for the selected gas.
        !   N_steps     (INT64)     : Number of simulation steps
        !   dt          (REAL64)    : Timestep (reduced units).
        !   L           (REAL64)    : Length of the simulation box. (N/dens)^(1/3).
        !   T           (REAL64)    : Selected temperature (reduced units).
        !   nu          (REAL64)    : Probability for the Andersen Thermostat.
        !   cutoff      (REAL64)    : Cutoff distance for the interaction (<0.5L).
        !   gdr_num_bins(INT64)     : Number of bins for calculating RDF.
        !   n_sweeps    (INT64)     : Number of sweeps to compute for the MSD calculation.
        !   vcutoff     (REAL64)    : Cuttof for the verlet list buffer zone.
        ! 
        ! Inout:
        !    r      (REAL64[3,N]) : positions of all N particles, in reduced units.
        !    v      (REAL64[3,N]) : velocities of all N particles, in reduced units.

        implicit none
        integer(kind=i64), intent(in)                :: log_unit, traj_unit, rdf_unit, msd_unit, N_steps, gdr_num_bins, &
                                                        n_sweeps, write_log, write_pos, imin, imax, local_N
        integer(kind=i32), intent(in)                :: irank
        real(kind=dp), intent(in)                    :: lj_epsilon,lj_sigma,mass, L,cutoff,T,nu,dt,vcutoff
        real(kind=dp), dimension(:,:), intent(inout) :: r,v
        integer(kind=i32), dimension(:), intent(in)  :: sendcounts, displs

        ! local variables
        real(kind=dp), dimension(:,:), allocatable   :: r0, rold, F, displacement!, rnew
        real(kind=dp)                                :: time, Etot,Epot,Ekin,Tinst,press,p_com_t,dr
        real(kind=dp), dimension(3)                  :: p_com, local_p_com
        real(kind=dp), dimension(:,:), allocatable   :: gdr, local_gdr
        real(kind=dp), dimension(:), allocatable     :: v_MSD, local_MSD
        real(kind=dp), dimension(:,:,:), allocatable :: time_r
        integer(kind=i64), dimension(:), allocatable :: vlist
        integer(kind=i64)                            :: i, N, d
        real(kind=dp)                                :: local_Ekin, local_Epot, local_virial, virial
        integer                                      :: ierror

        N = size(r,dim=2,kind=i64)
        dr = cutoff /dble(gdr_num_bins)

        allocate(r0(3,N))
        allocate(rold(3,N))
        allocate(F(3,N))
        allocate(gdr(2,gdr_num_bins))
        allocate(local_gdr(2, gdr_num_bins))
        allocate(vlist(N * local_N))

        allocate(displacement(3, local_N))
        allocate(time_r(n_sweeps,3,N))
        allocate(v_MSD(N_steps))
        allocate(local_MSD(N_steps))

        F = 0.0_dp
        rold = r
        !r(:, imin:imax) = r(:, imin:imax) - (L / 2.0_dp)

        local_MSD = 0.0d0
        displacement = 0.0_DP

        call compute_vlist(L, r, vcutoff, imin, imax, vlist)

        !call MPI_Barrier(MPI_COMM_WORLD, ierror)
        call g_r(local_gdr, r, 1_i64, gdr_num_bins, L, cutoff, N_steps, imin, imax, vlist)

        if (irank == 0) then
            write(log_unit, '(A)') "time  Etot  Epot  Ekin  Tinst  Pinst  Pt"
            write(output_unit,"(a)",advance='no') "Completed (%): "
        end if


        do i = 1, N_steps
            
            
            !choose integrator depending on user?
            ! call verlet_step(rnew, r, rold, v, F, dt, L, cutoff)
            call vv_integrator(r, v, cutoff, L, dt, imin, imax, vlist, displacement)
            ! call euler()
            call vel_Andersen(v,nu,T, imin, imax)

            
            if (mod(i, write_log) == 0) then
                
                ! Computation of local variables
                local_Ekin = calc_KE(v, imin, imax)
                local_Epot = calc_vdw_pbc(r, cutoff, L, imin, imax, vlist)
                local_virial = calc_pressure(L, r, cutoff, imin, imax, vlist)
                call compute_com_momenta(v, local_p_com, imin, imax)
                !local_p_com_t = dsqrt(dot_product(p_com,p_com))

                call MPI_Barrier(MPI_COMM_WORLD, ierror)
                call MPI_Reduce(local_Ekin, Ekin,     1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
                call MPI_Reduce(local_Epot, Epot,     1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
                call MPI_Reduce(local_virial, virial, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
                call MPI_Reduce(local_p_com, p_com,   3, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)

                if (irank == 0) then
                    
                    time = real(i, kind=dp)*dt
                    Ekin = Ekin * 0.5_DP
                    Epot = Epot * 0.5_DP  ! To account for the double countiung because of verlet lists
                    p_com_t = sqrt(dot_product(p_com, p_com))

                    Etot = Epot + Ekin
                    Tinst = calc_Tinst(Ekin, N)
                    press = (real(N, kind=dp)*Tinst / (L**3)) + ((1.0_dp / (3.0_dp * (L**3))) * virial)
                    
                    call writeSystem(log_unit,lj_epsilon,lj_sigma,mass, time,Etot,Epot,Ekin,Tinst,&
                    press,p_com_t)

                end if
            end if

            
            if (mod(i, write_pos) == 0) then 
                if (irank == 0) call writePositions(r, traj_unit)
            end if

            if (mod(i, N_steps/10) == 0) then
                if (irank == 0) write(output_unit, '(1x,i0)', advance='no') (100*i)/N_steps
            end if


            if (update_vlist(displacement, vcutoff)) then
                call compute_vlist(L, r, vcutoff, imin, imax, vlist)
                !write(output_unit, '(A,I2,A,I9,A,F12.8)') 'Worker ', irank, ' updating verlet list at step', i, &
                !' Mean number of neighbours per atom: ', &
                !real(count(vlist > 0_I64, kind=i64) - local_N, kind=dp) / real(local_N, kind=dp)
                displacement = 0.0_DP
            end if

            call MPI_Barrier(MPI_COMM_WORLD, ierror)
            
            do d = 1, 3
                call MPI_ALLGATHERV(r(d,imin:imax), int(local_N), MPI_DOUBLE_PRECISION, r(d,:), sendcounts, displs, &
                MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
                
                call MPI_ALLGATHERV(v(d,imin:imax), int(local_N), MPI_DOUBLE_PRECISION, v(d,:), sendcounts, displs, &
                MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
            end do

            if (i <= n_sweeps) then
                time_r(i,:,imin:imax) = r(:,imin:imax)
            end if

            call MSD(r, time_r, L, i, n_sweeps, local_MSD, imin, imax)
            call g_r(local_gdr, r, 2_i64, gdr_num_bins, L, cutoff, N_steps, imin, imax, vlist)
            ! call MPI_Barrier(MPI_COMM_WORLD, ierror)

        end do

        call MPI_Barrier(MPI_COMM_WORLD, ierror)
        call MPI_REDUCE(local_MSD, v_MSD, int(N_steps, kind=i32), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
        call MPI_REDUCE(local_gdr(2,:), gdr(2,:), int(gdr_num_bins, kind=i32), &
                        MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)

        if (irank == 0) then
            gdr(1,:) = local_gdr(1,:)
            call g_r(gdr, r, 3_i64, gdr_num_bins, L, cutoff, N_steps, imin, imax, vlist)
            call writeRDF(gdr, rdf_unit, lj_sigma)
            call writeMSD(v_MSD, msd_unit, n_sweeps, write_log, lj_sigma, lj_epsilon, dt, mass)
        end if

        deallocate(r0)
        deallocate(rold)
        deallocate(F)
        deallocate(vlist)
        deallocate(v_MSD)
        deallocate(local_MSD)
        deallocate(time_r)

    end subroutine mainLoop


    subroutine verlet_step(r_new, r, r_old, v, F, dt, L, cutoff, imin, imax, vlist)
        ! Performs one step of the integration using the Verlet Algorithm.
        !
        ! Args:
        !
        !   dt          (REAL64)    : Timestep (reduced units).
        !   L           (REAL64)    : Length of the simulation box. (N/dens)^(1/3).
        !   cutoff      (REAL64)    : Cutoff distance for the interaction (<0.5L).
        !
        ! Inout:
        !    r      (REAL64[3,N]) : positions of all N particles of the current iteration, in reduced units.
        !    rold   (REAL64[3,N]) : positions of all N particles of the last iteration, in reduced units.
        !    v      (REAL64[3,N]) : velocities of all N partciles, in reduced units.
        !    F      (REAL64[3,N]) : forces of all N partciles, in reduced units.
        ! Returns:
        !    rnew   (REAL64[3,N]) : positions of all N particles of the new iteration, in reduced units.

        ! In/Out variables
        real(kind=dp), dimension(:,:), intent(inout) :: r_old, v, F
        real(kind=dp), dimension(:,:), intent(inout) :: r
        real(kind=dp), dimension(:,:), intent(out)   :: r_new
        real(kind=dp), intent(in)                    :: dt, L, cutoff
        integer(kind=i64), intent(in)                :: imin, imax
        integer(kind=i64), dimension(:), intent(in)  :: vlist
        
        ! Computing F(t)
        call calc_vdw_force(r, cutoff, L, F, imin, imax, vlist)
        
        ! Setting r(t+dt)
        r_new = (2.0_dp * r) - r_old + (F * dt * dt)

        ! Appplying periodic boundary conditions
        call PBC(r_new, L)

        ! Computing velocities
        v = (r_new - r_old) / (2.0_dp * dt)

        ! Setting r(t-dt) -> r(t)
        r_old(:, :) = r(:, :)

    end subroutine verlet_step


    subroutine vv_integrator(positions, velocities, cutoff, L, dt, imin, imax, vlist, dsp)
        !
        !  Subroutine to update the positions of all particles using the Velocity Verlet
        ! algorithm. 

        ! Args:
        !    positions  (REAL64[3,N]) : positions of all N particles, in reduced units.
        !    velocities (REAL64[3,N]) : velocities of all N partciles, in reduced units.
        !    cutoff          (REAL64) : cutoff value of the interaction.
        !    L               (REAL64) : length of the sides of the box.
        !    dt              (REAL64) : value of the integration timestep.
        !    dsp  (REAL64[3,local_N]) : Accumulated displacement for each particle in charge for each rank
        
        ! Returns:
        !    positions  (REAL64[3,N]) : positions of all N particles, in reduced units.
        !    velocities (REAL64[3,N]) : velocities of all N partciles, in reduced units.

        implicit none
        double precision, dimension(:,:), intent(inout)      :: positions, velocities, dsp
        double precision, intent(in)                         :: cutoff, L, dt
        integer(kind=i64), intent(in)                        :: imin, imax
        integer(kind=i64),dimension(:), intent(in)           :: vlist
        ! local variables
        double precision, dimension(3, size(positions(1,:))) :: forces

        call calc_vdw_force(positions, cutoff, L, forces, imin, imax, vlist)
        
        dsp(:, :) = dsp(:, :) - positions(:, imin:imax)
        positions(:, imin:imax) = positions(:, imin:imax) + (dt*velocities(:, imin:imax)) + (0.5d0*dt*dt*forces(:, imin:imax))
        dsp(:, :) = dsp(:, :) + positions(:, imin:imax)

        call PBC(positions, L)

        velocities(:, imin:imax) = velocities(:, imin:imax) + (0.5d0*dt*forces(:, imin:imax))

        call calc_vdw_force(positions, cutoff, L, forces, imin, imax, vlist)
        velocities(:, imin:imax) = velocities(:, imin:imax) + (0.5d0*dt*forces(:, imin:imax))

    end subroutine vv_integrator



    subroutine BoxMuller(sgm, x1, x2, xout)
        ! Generates a random number that follows a normal distribution
        ! Args:
        !   sgm     (REAL64) : Normal distirbution sigma.
        !   x1, x2  (REAL64) : Random numbers.
        ! Returns:
        ! xout      (REAL64) : Random number from ND.

        implicit none
            real(kind=dp), intent(in) :: sgm, x1, x2
            real(kind=dp), intent(out) :: xout
            real(kind=dp), parameter :: pi = 4d0*datan(1d0)

            xout=sgm*dsqrt(-2d0*(dlog(1d0-x1)))*dcos(2d0*pi*x2)

    end subroutine BoxMuller


    subroutine vel_Andersen(vel,nu,temp, imin, imax)
        ! Andersen thermostat, changes the velocities in a system with 
        ! a certain probability that depends on the temperature
        !
        ! Args:
        !   nu      (REAL64)      : Probability for the Andersen Thermostat.
        !   temp    (REAL64)      : Temperature.
        ! Inout:
        !   vel     (REAL64[3,N]) : velocities of all N partciles, in reduced units.

        implicit none
        real(kind=dp), intent(in)                    :: nu,temp
        real(kind=dp), dimension(:,:), intent(inout) :: vel
        real(kind=dp)                                :: sig, nurand, x1, x2
        integer(kind=i64)                            :: i,k,N
        integer(kind=i64), intent(in)                :: imin, imax
        

        N = size(vel,dim=2,kind=i64)
        sig = dsqrt(temp) !temperature t is a parameter defined in parameters.f90
        
        do i=imin, imax
        ! a random number is generated for every particle,
        ! only if this number < nu, the particle's velocity is changed

            call random_number(nurand)
            if (nurand < nu) then
                do k=1,3
                    call random_number(x1)
                    call random_number(x2)

                    call boxmuller(sig, x1, x2, vel(k,i))
                enddo
            endif
        enddo
    end subroutine vel_Andersen

end module integrators

