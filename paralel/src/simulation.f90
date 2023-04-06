module simulation
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use            :: periodic_bc, only: PBC

    implicit none

    public  :: MSD, g_r

contains

    subroutine MSD(positions, time_positions, L, i_sweep, sweeps, v_MSD, imin, imax)
        !
        !  This function calculates the value of the Mean Squared Displacement from
        ! the values of the current positions and the initial positions of the
        ! particles.

        ! Args:
        !   positions      (REAL64[3,N])       : The current positions' matrix of the system.
        !   time_positions (REAL64[sweeps,3,N]): Positions' matrix of the system for different timeframes.
        !   L              (REAL64)            : Length of the sides of the box, in reduced units.
        !   i_sweep        (INT64)             : Current time iteration.
        !   sweeps         (INT64)             : Total sweeps to compute.
        !   v_MSD          (REAL64[N_steps])   : MSD value for each time-frame
        !   imin           (INT64)             : Minimum index of particle
        !   imax           (INT64)             : Maximum index of particle

        ! Returns:
        !   v_MSD          (REAL64[N_steps])   : MSD value for each time-frame
        !

        double precision, dimension(:,:), intent(in)        :: positions
        double precision, dimension(:,:,:), intent(in)      :: time_positions
        double precision, intent(in)                        :: L
        integer(kind=i64), intent(in)                       :: i_sweep, sweeps, imin, imax
        double precision, dimension(:), intent(inout)       :: v_MSD
        ! local variables
        double precision, dimension(3,size(positions(1,:))) :: r_aux
        double precision                                    :: MSD_aux
        integer(kind=i64)                                   :: i, k, N, ind

        N = size(positions(1,:), kind=i64)

        ind = sweeps
        if (i_sweep < sweeps) then
            ind = i_sweep
        end if

        do k = 1, ind
            r_aux = positions - time_positions(k,:,:)
            call PBC(r_aux, L)
            MSD_aux = 0.0d0

            do i = imin, imax
                MSD_aux = MSD_aux + sum(r_aux(:,i)*r_aux(:,i))
            end do

            MSD_aux = MSD_aux / dble(N)
            v_MSD(i_sweep-k+1) = v_MSD(i_sweep-k+1) + MSD_aux

        end do

    end subroutine MSD


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
