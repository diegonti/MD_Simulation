module integrators
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use potential_m, only: calc_KE, calc_pressure, calc_vdw_force, calc_vdw_pbc
    use periodic_bc, only: PBC
    use simulation, only: MSD

    implicit none
    public :: verlet, vel_verlet,andersen_thermostat

contains

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
            write(2, '(5(f12.8,1x))') i_ts*dt,kin_ene,u_ene,pres,ms_dist

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
            write(2, '(4(f12.8,1x))') i_ts*dt,kin_ene,u_ene,pres

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

    ! Declare the subroutine that implements the Andersen thermostat
    subroutine andersen_thermostat(natoms, temperature, velocities, dt, collision_freq, random_number_seed)
        implicit none
        integer(kind=i64), intent(in) :: natoms, random_number_seed
        real(kind=dp), intent(in) :: temperature, dt, collision_freq
        real(kind=dp), intent(inout), dimension(:,:) :: velocities
        real(kind=dp) :: random_number
        integer(kind=i64) :: i, j

        real(kind=dp), parameter :: v_max = 4.0_dp * sqrt(temperature) ! Maximum velocity deviation
        integer(kind=i64) :: ncollisions = 0 ! Counter for the number of collisions

        ! Initialize the random number generator
        call srand(random_number_seed)
        ! Loop over all particles in the system
        do j = 1, natoms
            ! Generate a random number between 0 and 1
            random_number = rand()
            ! Check if the particle should undergo a collision
            if (random_number < collision_freq * dt) then
                ! Increment the collision counter
                ncollisions = ncollisions + 1
                ! Generate a random velocity for the particle
                do i = 1, 3
                    velocities(i,j) = v_max * (2.0_dp * rand() - 1.0_dp)
                end do
            end if
        end do
    end subroutine andersen_thermostat

end module integrators
