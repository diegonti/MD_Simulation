module potential_m
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use            :: periodic_bc,     only: pbc
    implicit none

    public :: calc_pressure, calc_KE, calc_vdw_pbc, calc_vdw_force, update_vlist

contains

    function calc_pressure(lenth, positions, cutoff, imin, imax, vlist) result(virial)
        implicit none
        ! Author: Marc Alsina <marcalsinac@gmail.com>
        ! Subroutine to compute the pressure applying the virial theorem. 
        
        ! It combinesthe kinetick term plus the interaction term. The later depends 
        ! on the force from the interaction between pairs of molecules, which is computed
        ! using the Lennard Jones potential.

        ! Args:
        !   lenth     (REAL64)     : The side of the periodic box, in reduced units
        !   positions (REAL64[3,N]): The position's matrix of the system, in reduced units
        !   temp      (REAL64)     : The instant temperature of the system, in reduced units
        !   cutoff    (REAL64)     : The cut-off set to account por pair interactions, in reduced units
        !   imin      (INT64)      : The minimum index particle
        !   imax      (INT64)      : The maximum index particle
        !   vlist     (INT64[:])   : Verlet list for the rank who's processing

        ! Returns:
        !   press     (REAL64)     : The instant pressure at a given tempersture and atom's positions
        !                            in reduced units


        ! In/Out variables
        real(kind=dp), dimension(:,:), intent(in)   :: positions
        real(kind=dp), intent(in)                   :: lenth, cutoff
        real(kind=dp)                               :: press
        integer(kind=i64), intent(in)               :: imin, imax
        integer(kind=i64), dimension(:), intent(in) :: vlist 
        ! Internal variables
        integer(kind=i64)                           :: i_part, j_part, n_p, counter, neigh_idx, nneigh
        real(kind=dp), dimension(3)                 :: fij
        real(kind=dp), dimension(3,1)               :: rij
        real(kind=dp)                               :: dij, cutoff2, virial, vol

        press = 0.0_dp
        virial = 0.0_dp
        cutoff2 = cutoff * cutoff
        n_p = size(positions, dim=2, kind=i64)
        vol = lenth ** 3

        counter = 1_i64

        do i_part = imin, imax
            
            nneigh = vlist(counter)
            
            do neigh_idx = 1, nneigh

                j_part = vlist(counter + neigh_idx)

                ! Calculem rij
                rij(:,1) = positions(:, i_part) - positions(:, j_part)
                !rij(2,1) = positions(2, i_part) - positions(2, j_part)
                !rij(3,1) = positions(3, i_part) - positions(3, j_part)

                call pbc(rij, lenth)
                
                ! Calculem fij
                dij = (rij(1,1)**2) + (rij(2,1)**2) + (rij(3,1)**2)

                if (dij < cutoff2) then

                    fij(:) = (48.0_dp / dij**7 - 24.0_dp / dij**4) * rij(:,1)
                    !fij(2) = (48.0_dp / dij**7 - 24.0_dp / dij**4) * rij(2,1)
                    !fij(3) = (48.0_dp / dij**7 - 24.0_dp / dij**4) * rij(3,1)

                    ! Upgredagem el valor de la pressio
                    virial = virial + dot_product(rij(:,1), fij)
                end if

            end do

            counter = counter + nneigh + 1_i64

        end do
        
        ! press = (real(n_p, kind=dp)*temp / vol) + ((1.0_dp / (3.0_dp * vol)) * virial)
    end function calc_pressure

    pure function calc_KE(vel, imin, imax) result(ke)
        implicit none
        ! Author: Marc Alsina <marcalsinac@gmail.com>
        ! Subroutine to compute the kinetick energy
        !
        ! Given the atom's velocities, computes the kinetick energy associated on the system,
        ! in reduced units

        ! Args:
        !   vel (REAL64[3,N]): Atom velocities, in reduced units
        !
        ! Returns:
        !   ke  (REAL64)     : Kinetick energy, in reduced units


        ! In/Out variables
        real(kind=DP)                             :: ke
        real(kind=DP), dimension(:,:), intent(in) :: vel
        integer(kind=i64), intent(in)             :: imin, imax
        ! Internal variables
        integer(kind=I64)                         :: i

        ! Variable initialization
        ke = 0.0_DP
        
        do i = imin, imax
            ke  = ke + sum(vel(:, i)**2)
        end do

        ! Final multiplications
        ! ke = ke * 0.5_DP
    end function calc_KE

    function calc_vdw_pbc(pos, cutoff, boundary, imin, imax, vlist) result(vdw_calc)
        implicit none
        ! Author: Marc Alsina <marcalsinac@gmail.com>
        ! Subroutine to compute the potential energy by means of the Lennard-Jones potential
        !
        ! Computes the potential energy, in reduced units, for a system of identical atoms
        ! using the Lennard-Jones potential
        !
        ! Args:
        !   pos      (REAL64[3,N]): Positions of N atoms, in reduced units
        !   cutoff   (REAL64)     : Distance cut-off to consider interaction, in reduced units
        !   boundary (REAL64)     : Lenth of the periodic box, in reduced units
        !
        ! Returns
        !   vdw_calc (REAL64)     : Potential energy associated to the instant system positions,
        !                           in reduced units


        ! In/Out variables
        real(kind=DP), intent(in), dimension(:, :)  :: pos
        real(kind=DP)                               :: vdw_calc
        real(kind=DP), intent(in)                   :: cutoff, boundary
        integer(kind=i64), intent(in)               :: imin, imax
        integer(kind=i64), dimension(:), intent(in) :: vlist
        ! Internal variables
        real(kind=DP)                               :: dist, ecalc, cutoff2
        integer(kind=I64)                           :: i, j, num_particles, counter, nneigh, neigh_idx
        real(kind=DP), dimension(3,1)               :: rij

        ! Initializing parameters
        vdw_calc = 0.0_DP
        num_particles = size(pos, dim=2, kind=I64)

        cutoff2 = cutoff ** 2

        counter = 1_i64

        do j = imin, imax

            nneigh = vlist(counter)

            do neigh_idx = 1, nneigh

                i = vlist(counter + neigh_idx)

                rij(:,1) = pos(:, j) - pos(:, i)
                !rij(2,1) = pos(2, j) - pos(2, i)
                !rij(3,1) = pos(3, j) - pos(3, i)
                
                call pbc(rij, boundary)

                dist = (rij(1,1)**2) + (rij(2,1)**2) + (rij(3,1)**2)

                if (dist < cutoff2) then
                    
                    dist = dist ** 3
                    
                    ecalc = (4.0_DP * ((1.0_DP / (dist**2)) - (1.0_DP/(dist)))) - &
                            (4.0_DP * ((1.0_DP / (cutoff2**6)) - (1.0_DP/(cutoff2**3))))
                    
                    vdw_calc = vdw_calc + ecalc
                
                end if

            end do

            counter = counter + nneigh + 1_i64

        end do
    end function calc_vdw_pbc

    subroutine calc_vdw_force(pos, cutoff, boundary, forces, imin, imax, vlist)
        implicit none
        ! Author: Marc Alsina <marcalsinac@gmail.com>
        ! Subroutine to compute the force acting on each particel by means of the Lennard-Jones 
        ! potential
        !
        ! Computes the total force that acts on each atom, in reduced units, for a system 
        ! of identical atoms using the Lennard-Jones potential (i.e the spatial derivarive)
        !
        ! Args:
        !   pos      (REAL64[3,N]): Positions of N atoms, in reduced units
        !   cutoff   (REAL64)     : Distance cut-off to consider interaction, in reduced units
        !   boundary (REAL64)     : Lenth of the periodic box, in reduced units
        !
        ! Returns
        !   forces   (REAL64[3,N]): Total force acting on each atom, in reduced units.


        !In/Out variables
        real(kind=DP), intent(in), dimension(:, :)    :: pos
        real(kind=DP), intent(in)                     :: cutoff, boundary
        real(kind=DP), intent(out), dimension(:, :)   :: forces
        integer(kind=i64), intent(in)                 :: imin, imax
        integer(kind=i64), dimension(:), intent(in)   :: vlist
        ! Internal variables
        integer(kind=I64)             :: n_part, i, j, counter, neigh_idx, nneigh
        real(kind=DP)                 :: dist, cutoff2
        real(kind=DP), dimension(3,1) :: rij

        n_part = size(pos, dim=2, kind=I64)

        forces = 0.0_DP
        rij = 0.0_dp
        cutoff2 = cutoff * cutoff

        counter = 1_i64

        do i = imin, imax

            nneigh = vlist(counter)

            do neigh_idx = 1, nneigh

                j = vlist(counter + neigh_idx)
                
                rij(:,1) = pos(:, i) - pos(:, j)
                !rij(2,1) = pos(2, i) - pos(2, j)
                !rij(3,1) = pos(3, i) - pos(3, j)
                
                call pbc(rij, boundary)
                
                dist = (rij(1,1)**2) + (rij(2,1)**2) + (rij(3,1)**2)
                
                if (dist < cutoff) then
                    ! print *, i, j, dist
                    ! Calculem la forc entre particula i, j
                    forces(:, i) = forces(:, i) + (48.0_DP / dist**7 - 24.0_dp / dist**4) * rij(:,1)
                    !forces(2, i) = forces(2, i) + (48.0_DP / dist**7 - 24.0_dp / dist**4) * rij(2,1)
                    !forces(3, i) = forces(3, i) + (48.0_DP / dist**7 - 24.0_dp / dist**4) * rij(3,1)

                    !forces(1, j) = forces(1, j) - (48.0_DP / dist**7 - 24.0_dp / dist**4) * rij(1,1)
                    !forces(2, j) = forces(2, j) - (48.0_DP / dist**7 - 24.0_dp / dist**4) * rij(2,1)
                    !forces(3, j) = forces(3, j) - (48.0_DP / dist**7 - 24.0_dp / dist**4) * rij(3,1)

                end if
            end do

            counter = counter + nneigh + 1_i64

        end do
    end subroutine calc_vdw_force

    pure subroutine compute_com_momenta(vel, com_momenta, imin, imax)
        implicit none
        ! Author: Marc Alsina <marcalsinac@gmail.com>
        ! Subroutine to compute the center of mass momenta, in reduced units
        !
        ! Args:
        !   vel         (REAL64[3,N]): Velocity of the system, in reduced units
        !
        ! Returns:
        !   com_momenta (REAL64[3]): Center of mass momenta, in reduced units

        
        ! In/Out variables
        real(kind=DP), dimension(:, :), intent(in) :: vel
        real(kind=DP), dimension(3), intent(out)   :: com_momenta
        integer(kind=i64), intent(in)              :: imin, imax
        ! Internal variables
        integer(kind=I64)                          :: i_aux, n_p

        com_momenta = 0.0_DP
        n_p = size(vel, dim=2, kind=I64)

        do i_aux = imin, imax
            com_momenta(:) = com_momenta(:) + vel(:, i_aux)
        end do

    end subroutine compute_com_momenta

    pure function calc_Tinst(ke, np) result(tinst)
        implicit none
        ! Author: Marc Alsina <marcalsinac@gmail.com>
        ! Subroutine to compute the instantaneous temperature
        !
        ! Args:
        !   ke (REAL64): instant kinetic energy
        !   np (INT64) : Number of particles in the system
        !
        ! Returns:
        !   tinst (REAL64): Instant temperature

        
        ! In/Out variables
        real(kind=dp), intent(in)     :: ke
        integer(kind=i64), intent(in) :: np
        real(kind=dp)                 :: tinst

        tinst = (2.0_DP / (3.0_DP * real(np, kind=dp))) * ke

    end function calc_Tinst

    subroutine compute_vlist(L, r, cutoffv, imin, imax, vlist)
        ! Author: Marc Alsina <marcalsinac@gmail.com>
        ! Subroutine that computes the verlet list pertaining for each rank
        !
        ! Subroutine that computes the verlet list pertaining to each rank, 
        ! bounded by the particle range imin-imax.
        ! The update of the verlet list can be set at each n timesteps, or when the 
        ! outermost 
        !
        ! Args:
        !   r       (REAL64[3,N]) : Particle's positions
        !   cutoffv (REAL64)      : Verlet cuttof 
        !   imin    (INT64)       : Minimum index particle to compute
        !   imax    (INT64)       : Maximum index particle to compute
        !
        ! Returns:
        !   vlist   (INT64[:])    : Verlet list


        implicit none
        ! In/Out variables
        real(kind=dp), dimension(:,:), intent(in)    :: r
        real(kind=dp), intent(in)                    :: cutoffv, L
        integer(kind=i64), intent(in)                :: imin, imax
        integer(kind=i64), dimension(:), intent(out) :: vlist
        ! Internal variables
        integer(kind=i64)                            :: np, i, j, nneigh, pos
        real(kind=dp), dimension(3, 1)               :: rij
        real(kind=dp)                                :: dist2, cutoff2

        np = size(r, dim=2, kind=i64)
        cutoff2 = cutoffv * cutoffv

        pos = 1_I64

        do i = imin, imax

            nneigh = 0_I64

            do j = 1, np

                if (i == j) cycle

                rij(:, 1) = r(:, i) - r(:, j)
                call pbc(rij, L)
                dist2 = rij(1, 1)**2 + rij(2, 1)**2 + rij(3, 1)**2

                if (dist2 <= cutoff2) then
                    nneigh = nneigh + 1_I64
                    vlist(pos + nneigh) = j
                end if

            end do

            vlist(pos) = nneigh
            pos = pos + nneigh + 1_i64
        end do

    end subroutine compute_vlist

    function update_vlist(displacement, rmax) result(update)
        ! Author: Marc Alsina <marcalsinac@gmail.com>
        ! Subroutine that computes if the verlet list has to be updated,
        ! based on atoms displacements

        ! The actual criteria is to update the verlet list if the displacement
        ! of each particle is greater or equal half times the buffer cutoff.
        !
        ! Args:
        !   displacement (REAL64[3,N]): Displacement of each atom
        !   rmax         (REAL64)     : Verlet list cutoff
        !
        ! Returns:
        !   update       (LOGICAL)    : Update or not the verlet list
        implicit none
        ! In/Out variables
        real(kind=dp), dimension(:, :), intent(in) :: displacement
        real(kind=dp), intent(in)                  :: rmax
        logical                                    :: update


        update = any(displacement >= rmax * 0.5_DP)
    
    end function update_vlist

end module potential_m
