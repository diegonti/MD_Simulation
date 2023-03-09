module periodic_bc
        
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int16
    implicit none
    public :: PBC

contains

    subroutine PBC(positions, L)
        implicit none
        !
        ! Subroutine that applies periodic boundary conditions to a given set of particles.

        ! Args:
        !    positions (REAL64[3,N]) : The positions matrix of the system, in reduced units
        !    L              (REAL64) : Length of each side of the box.

        ! Returns:
        !    positions (REAL64[3,N]) : The positions matrix of the system, in reduced units
        !

        double precision, dimension(:,:), intent(inout) :: positions
        double precision, intent(in)                    :: L
        ! local variables
        integer(kind=i64)                               :: i, j, N, M

        N = size(positions(1,:),kind=i64)
        M = size(positions(:,1),kind=i64)

        do j = 1, M
            do i = 1, N
                if (positions(j,i) > L/2.d0) then
                    positions(j,i) = positions(j,i) - L
                else if (positions(j,i) < -L/2.d0) then
                    positions(j,i) = positions(j,i) + L
                end if
            end do
        end do
    end subroutine PBC

end module periodic_bc
