module simulation
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use            :: periodic_bc, only: PBC

    implicit none

    public  :: MSD, g_r

contains

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
        ! Author: Marc Alsina <marcalsinac@gmail.com>
        !
        ! Subroutine that performs a radial distribution finction
        !
        ! Args:
        !   gr_mat (REAL64[2,:]): The allocated radial distribution function
        !   pos    (REAL64[3,N]): Atoms positions at a given time step
        !   switch_case  (INT64): The selector for the calculation. 1: allocate mmemory,
        !                         2: compute the binning, 3: normalizie the binning
        !   num_bins     (INT64): Number of bins to use to generate the rdf
        !   L           (REAL64): Simulation side box
        !   cutoff      (REAL64): The cutoff to account for interactions
        !
        ! Returns:
        !   gr_mat (REAL64[2,:]): At case 1, it doesn't return nothing, at case 2,
        !                         it is the distance array selected, and the binning 
        !                         is performed. At stage 3, it returns in the first
        !                         dimension of rank 1 the dinstances, and the second 
        !                         RDF values.

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
