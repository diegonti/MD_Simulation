module integrators
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use periodic_bc, only: PBC
    use potential_m, only: calc_KE, calc_pressure, calc_vdw_force, calc_vdw_pbc, calc_Tinst, compute_com_momenta
    use simulation, only: MSD
    ! use writers_m, only: 

    implicit none
    public :: verlet, vel_verlet,vel_andersen
    private :: boxmuller

contains

    subroutine mainLoop(log_unit,traj_unit,rdf_unit,lj_epsilon,lj_sigma,mass,N_steps,dt,L,T,nu,cutoff,gdr_num_bins,r,v,F)
        ! Main Simulation Loop
        !

        implicit none
        integer(kind=i64), intent(in) :: log_unit,traj_unit,rdf_unit, N_steps, gdr_num_bins
        real(kind=dp), intent(in) :: lj_epsilon,lj_sigma,mass, L,cutoff,T,nu,dt
        real(kind=dp), dimension(:,:), intent(inout) :: r,v,F

        real(kind=dp), dimension(:,:), allocatable :: r0, rold, rnew
        real(kind=dp) :: time, Etot,Epot,Ekin,Tinst,press,rMSD,p_com_t, dr
        real(kind=dp), dimension(3) :: p_com
        real(kind=dp), dimension(gdr_num_bins) :: gdr
        integer(kind=i64) :: i,N

        dr = 1.5d0*L/dble(gdr_num_bins)
        N = size(r,dim=2,kind=i64)

        allocate(r0(3,N))
        allocate(rold(3,N))
        allocate(rnew(3,N))


        rold = r
        r0 = r  ! Saving initial configuration (for MSD)

        do i=1,N_steps
            time = i*dt
            !choose integrator depending on user?
            ! call verlet_step(rnew, r, rold, v, F, dt, L, cutoff)
            ! call vv_integrator(r, v, cutoff, L, dt)
            ! call euler()
            Epot = calc_vdw_pbc(r,cutoff,L)
            Ekin = calc_KE(v)
            Etot = Epot + Ekin
            Tinst =  calc_Tinst(Ekin,N)
            press =  calc_pressure(L, r, Tinst, cutoff)
            call compute_com_momenta(v,p_com)
            p_com_t = dsqrt(dot_product(p_com,p_com))

            rMSD = MSD(r,r0)
            call RDF(r,gdr,L,dr)

            call vel_Andersen(v,nu,T)
            r = rnew

            call writeSystem(log_unit,lj_epsilon,lj_sigma,mass, time,Etot,Epot,Ekin,T,press,rMSD,p_com_t)
            call writePositions(traj_unit,r)

        end do

        call writeRDF(dr,gdr,rdf_unit)

        deallocate(r0)
        deallocate(rold)
        deallocate(rnew)



    end subroutine mainLoop

    subroutine verlet(positions,velocities,dt,nts,cutoff,a_box,temp)
        implicit none
        ! subroutine that integrates newton's equations of motion by means of the verlet algorithm

        ! declaration of variables
        real(kind=dp), dimension(:,:), intent(inout)    :: positions,velocities
        real(kind=dp), intent(in)                       :: dt,cutoff,a_box,temp
        integer(kind=i64), intent(in)                   :: nts
        real(kind=dp)                                   :: u_ene,kin_ene,pres,ms_dist,dt2,rij(1,1),rij_aux
        real(kind=dp), dimension(:,:), allocatable      :: posold,forces
        integer(kind=i64)                               :: i,j,np,i_ts

        np = size(positions, dim=2, kind=i64)
        allocate(forces(3,np), posold(3,np))
        dt2 = dt*dt

        ! opens file to write trajectories
        !open(1, file='traj.xyz', status='new')
        !write(1,'(g0)') np
        !write(1,*) ''

        ! opens file to write energies, results for analysis
        open(2, file='energies.dat', status='new')
        ! timestep, total energy, potential energy, kinetic energy, temperature, pressure, MSD, total momentum
        write(2,'(8(a,1x))') 'ts','E','U','T','Temperature','P','MSD','p'
        
        posold = positions - velocities*dt ! fabricates positions at timestep -1 using initial velocities

        do i_ts = 1,nts

            ! gets current forces and energies, writes energies to file
            call calc_vdw_force(positions,cutoff,a_box,forces)
            kin_ene = calc_KE(velocities)
            u_ene = calc_vdw_pbc(positions,cutoff,a_box)
            pres = calc_pressure(a_box, positions, temp, cutoff)
            ms_dist = msd(posold,positions)
            write(2, '(5(f12.8,1x))') dble(i_ts)*dt,kin_ene,u_ene,pres,ms_dist

            do j=1,np
                do i=1,3
                    rij_aux = positions(i,j)
                    positions(i,j) = 2.0_dp*positions(i,j) - posold(i,j) + forces(i,j)*dt2 ! computes new positions
                    rij(1,1) = positions(i,j)
                    call pbc(rij, a_box) ! puts single component back to simulation box
                    positions(i,j) = rij(1,1)
                    velocities(i,j) = 0.5_dp*(positions(i,j)-posold(i,j))/dt ! new velocity by component
                    posold(i,j) = rij_aux
                enddo ! xyz components loop
                !write(1,'(a1,1x,3(f12.8,1x))') 'A', positions(1,j), positions(2,j), positions(3,j)
            enddo ! particles loop
            !write(1,*) ''

        enddo ! simulation loop
        close(2)
        !close(1)
    end subroutine

    subroutine vel_verlet(positions,velocities,dt,nts,cutoff,a_box,temp)
        implicit none
        ! subroutine that integrates newton's equations of motion by means of the velocity verlet algorithm

        real(kind=dp), dimension(:,:), intent(inout)    :: positions,velocities
        real(kind=dp), intent(in)                       :: dt, cutoff,a_box, temp
        integer(kind=i64), intent(in)                   :: nts
        real(kind=dp)                                   :: u_ene,kin_ene,pres,dt2,rij(1,1)
        real(kind=dp), dimension(:,:), allocatable      :: forces
        integer(kind=i64)                               :: i,j,np,i_ts

        np = size(positions, dim=2, kind=i64)
        allocate(forces(3,np))
        dt2 = dt*dt

        !open(1, file='traj.xyz', status='new')
        !write(1,'(g0)') np
        !write(1,*) ''

        open(2, file='energies.dat', status='new')
        write(2,'(4(a,1x))') 'ts','T','U','P' ! timestep, kinetic energy, potential energy, pressure

        do i_ts = 1,nts
            ! gets current forces and energies, writes output
            call calc_vdw_force(positions,cutoff,a_box,forces)
            kin_ene = calc_KE(velocities)
            u_ene = calc_vdw_pbc(positions,cutoff,a_box)
            pres = calc_pressure(a_box, positions, temp, cutoff)
            write(2, '(4(f12.8,1x))') dble(i_ts)*dt,kin_ene,u_ene,pres

            do j=1,np
                do i=1,3
                    positions(i,j) = positions(i,j) + velocities(i,j) + 0.5_dp*forces(i,j)*dt2
                    rij(1,1) = positions(i,j) 
                    call pbc(rij, a_box) ! puts component back to simulation box
                    positions(i,j) = rij(1,1)
                    velocities(i,j) = velocities(i,j) + 0.5_dp*forces(i,j)*dt ! velocity by component
                enddo ! xyz components loop
            enddo ! particles loop

            ! computes new forces, after moving all atoms, to get the new velocities
            call calc_vdw_force(positions,cutoff,a_box,forces)
            do j=1,np
                do i=1,3
                    velocities(i,j) = velocities(i,j) + 0.5_dp*forces(i,j)*dt ! new velocity by component
                enddo ! xyz components loop
            enddo ! particles loop
            !write(1,*) ''
        enddo ! simulation loop

        !close(1)
        close(2)
    end subroutine

    subroutine BoxMuller(sgm, x1, x2, xout)
    implicit none
    ! generates a random number that follows a normal distribution
        real(kind=dp), intent(in) :: sgm, x1, x2
        real(kind=dp), intent(out) :: xout
        real(kind=dp) :: pi
        pi = 4d0*datan(1d0)

        xout=sgm*dsqrt(-2d0*(dlog(1d0-x1)))*dcos(2d0*pi*x2)
    end subroutine BoxMuller

    subroutine vel_Andersen(vel,nu,temp)
    implicit none
    ! Andersen thermostat, changes the velocities in a system with a certain probability that depends on the temperature
        real(kind=dp), intent(in) :: nu,temp
        real(kind=dp), dimension(:,:), intent(inout) :: vel
        real(kind=dp) :: sig, nurand, x1, x2
        integer(kind=i64) :: i,k,N
        integer :: state_size
        integer, allocatable, dimension(:) :: state
        integer, parameter:: seed_number = 165432156

        N = size(vel,dim=2,kind=i64)
        sig = dsqrt(temp) !temperature t is a parameter defined in parameters.f90
        
        ! Setting Random Seed
        call random_seed( size=state_size )
        allocate(state(state_size))
        state = seed_number
        call random_seed( put=state )
        
        do i=1,N
        ! a random number is generated for every particle,
        ! only if this number < nu, the particle's velocity is changed

            call random_number(nurand)
            if (nurand < nu) then
                do k=1,3
                    call random_number(x1)
                    call random_number(x2)

                    call boxmuller(sig, x1, x2, vel(i,k))
                enddo
            endif
        enddo
    end subroutine

end module integrators
