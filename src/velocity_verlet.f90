module velocity_verlet

    use, intrinsic :: iso_fortran_env, only : dp => real64, i64 => int16
    use            :: periodic_bc, only : PBC
    use            :: potential_m, only : calc_vsw_force

    implicit none

    public :: vv_integrator

contains
    
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
        !

        double precision, dimension(:,:), intent(inout)      :: positions, velocities
        double precision, intent(in)                         :: cutoff, L, dt
        ! local variables
        integer(kind=i64)                                    :: i, j
        double precision, dimension(3, size(positions(1,:))) :: forces

        call calc_vdw_force(positions, cutoff, L, forces)
        
        positions = positions + dt*velocities + 0.5d0*dt*dt*forces
        call PBC(positions, L)

        velocities = velocities + 0.5d0*dt*forces

        call calc_vdw_force(positions, cutoff, L, forces)
        velocities = velocities + 0.5d0*dt*forces

    end subroutine vv_integrator

end module velocity_verlet
