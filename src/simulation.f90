module simulation
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
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

        double precision, dimension(:,:), intent(in)  :: positions
        double precision, dimension(:), intent(inout) :: distribution
        double precision, intent(in)                  :: L, dr
        ! local variables
        double precision, dimension(3, 1)             :: r_ij
        double precision, dimension(27)               :: neighbour_dist
        integer(kind=i64)                             :: i, j, k, N

        N = size(positions(1,:), kind=i64)

        do i = 1, N
            do j = i+1_i64, N
                r_ij(:,1) = positions(:,i) - positions(:,j)
                call PBC(r_ij, L)
                call neighbour_distances(r_ij(:,1), L, neighbour_dist)
                do k = 1, 27
                    if (neighbour_dist(k) <= 1.5d0*L) then
                        distribution(int(neighbour_dist(k)/dr)) = &
                        distribution(int(neighbour_dist(k)/dr)) + 1.0d0
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


    function MSD(positions, initial_positions, L)
        !
        !  This function calculates the value of the Mean Squared Displacement from
        ! the values of the current positions and the initial positions of the 
        ! particles. 

        ! Args:
        !   positions         (REAL64[3,N]) : The current positions' matrix of the system.
        !   initial_positions (REAL64[3,N]) : The initial positions' matrix of the system.
        !   L                      (REAL64) : Length of the sides of the box, in reduced units.
        
        ! Returns:
        !   MSD (REAL64[:]) : Value of the mean squared displacement of the system.
        !

        double precision, dimension(:,:), intent(in)        :: positions, initial_positions
        double precision, intent(in)                        :: L
        double precision                                    :: MSD
        ! local variables
        double precision, dimension(3,size(positions(1,:))) :: r_aux
        integer(kind=i64)                                   :: i, N

        N = size(positions(1,:), kind=i64)

        MSD = 0.0d0

        r_aux = positions - initial_positions
        call PBC(r_aux, L)

        do i = 1, N
            MSD = MSD + sum(r_aux(:,i)*r_aux(:,i))  
        end do
        
        MSD = MSD / dble(N)

    end function MSD

    subroutine g_r(gr_mat, pos, switch_case, num_bins, L, cutoff)
        ! Notes
        ! gr_mat(1,:) -> valors de distancia
        ! gr_mat(2,:) -> numero de partícules a aquesta distancia
        implicit none
        ! In/Out variables
        real(kind=dp), dimension(:,:), intent(inout) :: gr_mat
        real(kind=dp), dimension(:,:), intent(in)    :: pos
        real(kind=dp), intent(in)                    :: L, cutoff
        integer(kind=i64), intent(in)                :: switch_case, num_bins
        ! Internal variables
        integer(kind=i64), save                      :: i_ax, index_mat, j_ax, n_p, n_gdr
        real(kind=dp), save                          :: dr, dist, dv, ndg, dens
        real(kind=dp), parameter                     :: PI = 4.0_dp * atan(1.0_dp)
        real(kind=dp), dimension(3, 1)               :: rij

        select case (switch_case)
            case (1_i64)
                
                ! SWITCH = 1 => definim la memoria de la funció
                n_p = size(pos, dim=2, kind=i64)
                dens = real(n_p, kind=dp) / (L ** 3)
                dr = cutoff / real(num_bins, kind=dp)

                gr_mat(1,:) = [(real(i_ax, kind=dp)*dr, i_ax=1, num_bins)]
                gr_mat(2,:) = 0.0_dp

                n_gdr = 0_i64
            
            case (2_i64)
                
                ! SWITCH = 2 => Calculem n(r)
                n_gdr = n_gdr + 1_i64
                do j_ax = 1, n_p - 1
                    do i_ax = j_ax + 1, n_p
                        ! Calculem rij
                        rij(:, 1) = pos(:, j_ax) - pos(:, i_ax)
                        call pbc(rij, L)
                        dist = norm2(rij(:, 1))
                        
                        ! Apliquem el cutoff de maxima distancia
                        if (dist < cutoff) then  ! MA: cut-off is set as the maximum distance to account
                            index_mat = int(dist/dr, kind=i64) + 1_i64
                            gr_mat(2, index_mat) = gr_mat(2, index_mat) + 2.0_dp
                        end if
                    end do
                end do
            
            case (3_i64)
        
                ! SWITCH = 3 => Calculem g(r)
                do i_ax = 1, num_bins
                    associate(gdr => gr_mat(2, i_ax))
                        dv = (((real(i_ax, kind=dp) + 1.0_dp) ** 3) - (real(i_ax, kind=dp) ** 3)) * (dr ** 3)
                        ndg = (4.0_dp / 3.0_dp) * pi * dv * dens
                        gdr = gdr / (real(n_p, kind=dp) * ndg * real(n_gdr, kind=dp))
                    end associate
                end do
            
            end select

    end subroutine g_r

end module simulation
