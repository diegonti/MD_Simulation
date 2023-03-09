module integrators
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use periodic_bc, only: PBC
    use potential_m, only: calc_KE, calc_pressure, calc_vdw_force, calc_vdw_pbc, calc_Tinst, compute_com_momenta
    use simulation,  only: MSD, g_r
    use writers_m,   only: writePositions, writeRdf, writeSystem
    use testing

    implicit none
    public :: verlet_step, vv_integrator,vel_andersen
    private :: boxmuller

contains

    subroutine mainLoop(log_unit,traj_unit,rdf_unit,lj_epsilon,lj_sigma,mass,N_steps,dt,L,T,nu,cutoff,gdr_num_bins,r,v)
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
        ! 
        ! Inout:
        !    r      (REAL64[3,N]) : positions of all N particles, in reduced units.
        !    v      (REAL64[3,N]) : velocities of all N partciles, in reduced units.

        implicit none
        integer(kind=i64), intent(in) :: log_unit,traj_unit,rdf_unit, N_steps, gdr_num_bins
        real(kind=dp), intent(in) :: lj_epsilon,lj_sigma,mass, L,cutoff,T,nu,dt
        real(kind=dp), dimension(:,:), intent(inout) :: r,v

        real(kind=dp), dimension(:,:), allocatable :: r0, rold, rnew, F
        real(kind=dp) :: time, Etot,Epot,Ekin,Tinst,press,rMSD,p_com_t, dr
        real(kind=dp), dimension(3) :: p_com
        real(kind=dp), dimension(:,:), allocatable :: gdr
        integer(kind=i64) :: i,N

        dr = 1.5d0*L/dble(gdr_num_bins)
        N = size(r,dim=2,kind=i64)

        allocate(r0(3,N))
        allocate(rold(3,N))
        allocate(rnew(3,N))
        allocate(F(3,N))
        allocate(gdr(2,gdr_num_bins))


        F = 0.0_dp
        rold = r
        r0 = r  ! Saving initial configuration (for MSD)
        r = r - (L / 2.0_dp)

        call g_r(gdr, r, 1_i64, gdr_num_bins, L, cutoff)

        write(log_unit, '(A)') "time  Etot  Epot  Ekin  Tinst  Pinst  MSD Pt"

        do i=1,N_steps
            time = real(i, kind=dp)*dt
            !choose integrator depending on user?
            ! call verlet_step(rnew, r, rold, v, F, dt, L, cutoff)
            call vv_integrator(r, v, cutoff, L, dt)
            ! call euler()
            call vel_Andersen(v,nu,T)

            Epot = calc_vdw_pbc(r,cutoff,L)
            Ekin = calc_KE(v)
            Etot = Epot + Ekin
            Tinst =  calc_Tinst(Ekin,N)
            press =  calc_pressure(L, r, Tinst, cutoff)
            call compute_com_momenta(v,p_com)
            p_com_t = dsqrt(dot_product(p_com,p_com))

            rMSD = MSD(r,r0,L)
            call g_r(gdr, r, 2_i64, gdr_num_bins, L, cutoff)
            
            ! r = rnew

            call writeSystem(log_unit,lj_epsilon,lj_sigma,mass, time,Etot,Epot,Ekin,Tinst,press,rMSD,p_com_t)
            call writePositions(r, traj_unit)

        end do

        call g_r(gdr, r, 3_i64, gdr_num_bins, L, cutoff)

        call writeRDF(gdr,rdf_unit)

        deallocate(r0)
        deallocate(rold)
        deallocate(rnew)
        deallocate(F)


    end subroutine mainLoop


    subroutine verlet_step(r_new, r, r_old, v, F, dt, L, cutoff)
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
        
        ! Computing F(t)
        call calc_vdw_force(r, cutoff, L, F)
        
        ! Setting r(t+dt)
        r_new = (2.0_dp * r) - r_old + (F * dt * dt)

        ! Appplying periodic boundary conditions
        call PBC(r_new, L)

        ! Computing velocities
        v = (r_new - r_old) / (2.0_dp * dt)

        ! Setting r(t-dt) -> r(t)
        r_old(:, :) = r(:, :)

    end subroutine verlet_step


    subroutine vv_integrator(positions, velocities, cutoff, L, dt)
        !
        !  Subroutine to update the positions of all particles using the Velocity Verlet
        ! algorithm. 

        ! Args:
        !    positions  (REAL64[3,N]) : positions of all N particles, in reduced units.
        !    velocities (REAL64[3,N]) : velocities of all N partciles, in reduced units.
        !    cutoff          (REAL64) : cutoff value of the interaction.
        !    L               (REAL64) : length of the sides of the box.
        !    dt              (REAL64) : value of the integration timestep.
        
        ! Returns:
        !    positions  (REAL64[3,N]) : positions of all N particles, in reduced units.
        !    velocities (REAL64[3,N]) : velocities of all N partciles, in reduced units.

        implicit none
        double precision, dimension(:,:), intent(inout)      :: positions, velocities
        double precision, intent(in)                         :: cutoff, L, dt
        ! local variables
        double precision, dimension(3, size(positions(1,:))) :: forces

        call calc_vdw_force(positions, cutoff, L, forces)
        
        positions = positions + (dt*velocities) + (0.5d0*dt*dt*forces)
        call PBC(positions, L)

        velocities = velocities + (0.5d0*dt*forces)

        call calc_vdw_force(positions, cutoff, L, forces)
        velocities = velocities + 0.5d0*dt*forces

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


    subroutine vel_Andersen(vel,nu,temp)
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
        

        N = size(vel,dim=2,kind=i64)
        sig = dsqrt(temp) !temperature t is a parameter defined in parameters.f90
        
        do i=1,N
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


!!! NOT USED

! subroutine verlet(positions,velocities,dt,nts,cutoff,a_box,temp)
!     implicit none
!     ! _ Not used, just in case
!     ! subroutine that integrates newton's equations of motion by means of the verlet algorithm

!     ! declaration of variables
!     real(kind=dp), dimension(:,:), intent(inout)    :: positions,velocities
!     real(kind=dp), intent(in)                       :: dt,cutoff,a_box,temp
!     integer(kind=i64), intent(in)                   :: nts
!     real(kind=dp)                                   :: u_ene,kin_ene,pres,ms_dist,dt2,rij(1,1),rij_aux
!     real(kind=dp), dimension(:,:), allocatable      :: posold,forces
!     integer(kind=i64)                               :: i,j,np,i_ts

!     np = size(positions, dim=2, kind=i64)
!     allocate(forces(3,np), posold(3,np))
!     dt2 = dt*dt

!     ! opens file to write trajectories
!     !open(1, file='traj.xyz', status='new')
!     !write(1,'(g0)') np
!     !write(1,*) ''

!     ! opens file to write energies, results for analysis
!     open(2, file='energies.dat', status='new')
!     ! timestep, total energy, potential energy, kinetic energy, temperature, pressure, MSD, total momentum
!     write(2,'(8(a,1x))') 'ts','E','U','T','Temperature','P','MSD','p'
    
!     posold = positions - velocities*dt ! fabricates positions at timestep -1 using initial velocities

!     do i_ts = 1,nts

!         ! gets current forces and energies, writes energies to file
!         call calc_vdw_force(positions,cutoff,a_box,forces)
!         kin_ene = calc_KE(velocities)
!         u_ene = calc_vdw_pbc(positions,cutoff,a_box)
!         pres = calc_pressure(a_box, positions, temp, cutoff)
!         ms_dist = msd(posold,positions)
!         write(2, '(5(f12.8,1x))') dble(i_ts)*dt,kin_ene,u_ene,pres,ms_dist

!         do j=1,np
!             do i=1,3
!                 rij_aux = positions(i,j)
!                 positions(i,j) = 2.0_dp*positions(i,j) - posold(i,j) + forces(i,j)*dt2 ! computes new positions
!                 rij(1,1) = positions(i,j)
!                 call pbc(rij, a_box) ! puts single component back to simulation box
!                 positions(i,j) = rij(1,1)
!                 velocities(i,j) = 0.5_dp*(positions(i,j)-posold(i,j))/dt ! new velocity by component
!                 posold(i,j) = rij_aux
!             enddo ! xyz components loop
!             !write(1,'(a1,1x,3(f12.8,1x))') 'A', positions(1,j), positions(2,j), positions(3,j)
!         enddo ! particles loop
!         !write(1,*) ''

!     enddo ! simulation loop
!     close(2)
!     !close(1)
! end subroutine

! subroutine vel_verlet(positions,velocities,dt,nts,cutoff,a_box,temp)
!     ! _ Not used, just in case
!     implicit none
!     ! subroutine that integrates newton's equations of motion by means of the velocity verlet algorithm

!     real(kind=dp), dimension(:,:), intent(inout)    :: positions,velocities
!     real(kind=dp), intent(in)                       :: dt, cutoff,a_box, temp
!     integer(kind=i64), intent(in)                   :: nts
!     real(kind=dp)                                   :: u_ene,kin_ene,pres,dt2,rij(1,1)
!     real(kind=dp), dimension(:,:), allocatable      :: forces
!     integer(kind=i64)                               :: i,j,np,i_ts

!     np = size(positions, dim=2, kind=i64)
!     allocate(forces(3,np))
!     dt2 = dt*dt

!     !open(1, file='traj.xyz', status='new')
!     !write(1,'(g0)') np
!     !write(1,*) ''

!     open(2, file='energies.dat', status='new')
!     write(2,'(4(a,1x))') 'ts','T','U','P' ! timestep, kinetic energy, potential energy, pressure

!     do i_ts = 1,nts
!         ! gets current forces and energies, writes output
!         call calc_vdw_force(positions,cutoff,a_box,forces)
!         kin_ene = calc_KE(velocities)
!         u_ene = calc_vdw_pbc(positions,cutoff,a_box)
!         pres = calc_pressure(a_box, positions, temp, cutoff)
!         write(2, '(4(f12.8,1x))') dble(i_ts)*dt,kin_ene,u_ene,pres

!         do j=1,np
!             do i=1,3
!                 positions(i,j) = positions(i,j) + velocities(i,j) + 0.5_dp*forces(i,j)*dt2
!                 rij(1,1) = positions(i,j) 
!                 call pbc(rij, a_box) ! puts component back to simulation box
!                 positions(i,j) = rij(1,1)
!                 velocities(i,j) = velocities(i,j) + 0.5_dp*forces(i,j)*dt ! velocity by component
!             enddo ! xyz components loop
!         enddo ! particles loop

!         ! computes new forces, after moving all atoms, to get the new velocities
!         call calc_vdw_force(positions,cutoff,a_box,forces)
!         do j=1,np
!             do i=1,3
!                 velocities(i,j) = velocities(i,j) + 0.5_dp*forces(i,j)*dt ! new velocity by component
!             enddo ! xyz components loop
!         enddo ! particles loop
!         !write(1,*) ''
!     enddo ! simulation loop

!     !close(1)
!     close(2)
! end subroutine

