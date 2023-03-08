module verlet_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use            :: periodic_bc,     only: PBC
    use            :: potential_m,     only: calc_vdw_force
    implicit none

contains

    subroutine verlet_step(r_new, r, r_old, v, F, dt, a_box, cutoff)
        implicit none
        ! In/Out variables
        real(kind=dp), dimension(:,:), intent(inout) :: r_old, v, F
        real(kind=dp), dimension(:,:), intent(inout) :: r
        real(kind=dp), dimension(:,:), intent(out)   :: r_new
        real(kind=dp), intent(in)                    :: dt, a_box, cutoff
        ! Internal variables
        
        ! Computing F(t)
        call calc_vdw_force(r, cutoff, a_box, F)
        
        ! Setting r(t+dt)
        r_new = (2.0_dp * r) - r_old + (F * dt * dt)

        ! Appplying periodic boundary conditions
        call PBC(r_new, a_box)

        ! Computing velocities
        v = (r_new - r_old) / (2.0_dp * dt)

        ! Setting r(t-dt) -> r(t)
        r_old(:, :) = r(:, :)

    end subroutine verlet_step
end module verlet_mod