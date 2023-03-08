module simulation
        
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int16
    use            :: periodic_bc, only: PBC

    implicit none

    public  :: RDF, MSD
    private :: neighbour_distances

contains

    subroutine RDF(positions, distribution, L, dr)
        !
        !  This subrotine calculates the Radial Distribution Function of the 
        ! given configuration of the system; and, if the 'distribution' array
        ! is non-zero, the current RDF will be added to the previous one, stored
        ! in 'distribution'. This way, the final result of the RDF will give us
        ! a more precise description of the real distribution (convenient if the
        ! number of particles is not large).
        
        ! Args:
        !   positions  (REAL64[3,N]) : The positions matrix of the system (reduced units).
        !   distribution (REAL64[:]) : The current RDF of the system.
        !   L               (REAL64) : Length of the sides of the box.
        !   dr              (REAL64) : Value of the differential of the length of the
        !                               distribution

        ! Returns:
        !   distribution (REAL64[:]) : The current RDF of the system.
        !

        double precision, dimension(:,:), intent(in)         :: positions
        double precision, dimension(:), intent(inout)        :: distribution
        double precision, intent(in)                         :: L, dr
        ! local variables
        double precision, dimension(size(positions(:,1)), 1) :: r_ij
        double precision, dimension(27)                      :: neighbour_dist
        integer(kind=i64)                                    :: i, j, k, N, M

        N = size(positions(1,:),kind=i64)
        M = size(positions(:,1),kind=i64)

        do i = 1, N
            do j = i+1_i64, N
                r_ij(:,1) = positions(:,i) - positions(:,j)
                call PBC(r_ij, L)
                call neighbour_distances(r_ij(:,1), L, neighbour_dist)
                do k = 1, 27
                    if (neighbour_dist(k) <= 1.5d0*L) then
                        distribution(int(neighbour_dist(k)/dr) + 1) = &
                        distribution(int(neighbour_dist(k)/dr) + 1) + 1.0d0
                    end if
                end do
            end do
        end do

    end subroutine RDF


    subroutine neighbour_distances(r_ij_vector, L, neighbour_dist)
        !
        !  The subroutine obtains the positions of the images of particle in r_ij_vector
        ! for every neighbour box and includes them, as well as the positions of the
        ! original particle, in the local array ('neighbour_dist').
        
        ! Args:
        !   r_ij_vector    (REAL64[3]) : The position of the original particle.
        !   L                 (REAL64) : Length of the sides of the box.
        
        ! Returns:
        !   neighbour_dist (REAL64[:]) : Values of the distances of the image particles
        !

        double precision, dimension(:), intent(in)     :: r_ij_vector
        double precision, intent(in)                   :: L
        double precision, dimension(27), intent(out)   :: neighbour_dist
        ! local variables
        integer(kind=i64)                              :: i, j, k, indx
        double precision, dimension(size(r_ij_vector)) :: r_aux

        do i = -1_i64, 1_i64
            do j = -1_i64, 1_i64
                do k = -1_i64, 1_i64
                    r_aux = r_ij_vector + (/real(i,kind=dp)*L, real(j,kind=dp)*L, real(k,kind=dp)*L/)
                    indx = 9_i64*(i+1_i64) + 3_i64*(j+1_i64) + k + 2_i64
                    neighbour_dist(indx) = dsqrt(sum(r_aux * r_aux))
                end do
            end do
        end do

    end subroutine neighbour_distances


    function MSD(positions, initial_positions)
        !
        !  This function calculates the value of the Mean Squared Displacement from
        ! the values of the current positions and the initial positions of the 
        ! particles. 

        ! Args:
        !   positions         (REAL64[3,N]) : The current positions' matrix of the system.
        !   initial_positions (REAL64[3,N]) : The initial positions' matrix of the system.
        
        ! Returns:
        !   MSD (REAL64[:]) : Value of the mean squared displacement of the system.
        !

        double precision, dimension(:,:), intent(in) :: positions, initial_positions
        double precision                             :: MSD
        ! local variables
        integer(kind=i64)                            :: i, N

        N = size(positions(1,:), kind=i64)

        MSD = 0.0d0

        do i = 1, N
            MSD = MSD + sum(initial_positions(:,i) - positions(:,i))*& 
                        sum(initial_positions(:,i) - positions(:,i))
        end do
        
        MSD = MSD / dble(N)

    end function MSD

end module simulation
