module simulation
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use            :: periodic_bc, only: PBC

    implicit none

    public  :: MSD, g_r

contains

    function MSD(positions, initial_positions, L, imin, imax)
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
        integer(kind=i64), intent(in)                       :: imin, imax
        double precision                                    :: MSD
        ! local variables
        double precision, dimension(3,imax-imin+1)          :: r_aux
        integer(kind=i64)                                   :: i

!        N = size(positions(1,:), kind=i64)

        MSD = 0.0d0

        r_aux = positions(:,imin:imax) - initial_positions(:,imin:imax)
        call PBC(r_aux, L)

        do i = 1, imax-imin
            MSD = MSD + sum(r_aux(:,i)*r_aux(:,i))  
        end do
        
        !MSD = MSD / dble(N)

    end function MSD

    subroutine g_r(gr_mat, pos, switch_case, num_bins, L, cutoff, imin, imax, vlist)
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
        !   imin         (INT64): The minimum index particle
        !   imax         (INT64): The maximum index particle
        !   vlist     (INT64[:]): Verlet list for the rank who's processing
        !
        ! Returns:
        !   gr_mat (REAL64[2,:]): At case 1, it doesn't return anything, at case 2,
        !                         it is the distance array selected, and the binning 
        !                         is performed.

        implicit none
        ! In/Out variables
        real(kind=dp), dimension(:), intent(inout)   :: gr_mat
        real(kind=dp), dimension(:,:), intent(in)    :: pos
        real(kind=dp), intent(in)                    :: L, cutoff
        integer(kind=i64), intent(in)                :: switch_case, num_bins, imin, imax
        integer(kind=i64), dimension(:), intent(in)  :: vlist
        ! Internal variables
        integer(kind=i64), save                      :: index_mat
        integer(kind=i64)                            :: counter, n_neigh, i_part, j_part, neigh_index
        real(kind=dp), save                          :: dr, dist
        real(kind=dp), dimension(3, 1)               :: rij

        select case (switch_case)
            case (1_i64)
                
                ! SWITCH = 1 => definim la memoria de la funció
                dr = cutoff / real(num_bins, kind=dp)

                gr_mat(:) = 0.0_dp
            
            case (2_i64)
                
                ! SWITCH = 2 => Calculem n(r)
                counter = 1_i64

                do j_part = imin, imax
                    n_neigh = vlist(counter)  ! nombre de veïns (a distancia L/2)

                    do neigh_index = 1, n_neigh
                        ! Calculem rij
                        i_part = vlist(counter + neigh_index)
                        if (i_part < j_part) then
                            cycle
                        end if
                       
                        rij(:, 1) = pos(:, j_part) - pos(:, i_part)
                        call pbc(rij, L)
                        dist = norm2(rij(:, 1))
                        
                        ! Apliquem el cutoff de maxima distancia
                        if (dist < cutoff) then  ! MA: cut-off is set as the maximum distance to account
                            index_mat = int(dist/dr, kind=i64) + 1_i64
                            gr_mat(index_mat) = gr_mat(index_mat) + 2.0d0
                        end if
                    end do
                    counter = counter + n_neigh  +1_i64 ! index de la següent particula
                end do
            
            end select

    end subroutine g_r

end module simulation
