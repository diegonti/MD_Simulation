module initialization
    ! Module with the initialization of the position and velocity matrices.
    ! Contains two general functions that choose the desired initialization method.
    ! 3 types of initial positions are implemented -> SC, BCC, FCC.
    ! 2 types of initial velocities are implemented -> bimodal, zero.

    use, intrinsic :: iso_fortran_env, only: DP => real64, I64 => int64, input_unit, output_unit
    contains

    subroutine getInitialParams(cell,N,density,M,L,a)
        ! Returns the initial cell parameters of the simulation.
        !
        ! Args: 
        !   cell    (CHARACTER) : String with the name of the cell (sc, bcc, fcc).
        !   N       (INT64)     : Number of particles.
        !   density (REAL64)    : Density of the system.
        !
        ! Returns:
        !   M       (INT64)     : Number of cells in one direction.
        !   L       (REAL64)    : Length of the simulation box.
        !   a       (REAL64)    : Lattice parameter of the unit cell (L/M).

        character(*), intent(in) :: cell
        integer(kind=i64), intent(in) :: N
        double precision, intent(in) :: density
        integer(kind=i64), intent(out) :: M
        double precision, intent(out) :: L, a
        
        L = (real(N,kind=dp)/density)**(1.d0/3.d0)
        if (cell=="sc") then; M = nint((real(N,kind=dp))**(1.d0/3.d0),kind=i64)
        else if (cell=="fcc") then; M = nint((real(N,kind=dp)/4.d0)**(1.d0/3.d0), kind=i64)
        else if (cell=="bcc") then; M = nint((real(N,kind=dp)/2.d0)**(1.d0/3.d0), kind=i64)
        else; print*, "Select a valid initial position: sc, bcc, fcc."
        end if
    
        a = L/real(M, kind=dp)            ! Lattice parameter
        
    end subroutine getInitialParams


    subroutine changeIUnits(lj_epsilon,lj_sigma,mass,density,dt,T)
        ! Changes user Input units to reduced units.
        !
        ! Args:
        !   lj_epsilon  (REAL64) : Lennard Jones epsilon parameter for the gas (in kJ/mol).
        !   lj_sigma    (REAL64) : Lennard Jones sigma parameter for the gas (in Ang).
        !   mass        (REAL64) : Molar Mass of the gas (g/mol).
        !
        ! Inout:
        !   density     (REAL64) : Density of the system (in g/mL).
        !   dt          (REAL64) :  Time-step of the simulation (in ps).
        !   T           (REAL64) : Temperature (in K).

        double precision, intent(in) :: lj_epsilon,lj_sigma,mass
        double precision, intent(inout) :: density,dt,T
        double precision :: ru_time,ru_dens,ru_temp
        double precision, parameter :: Na = 6.0221408d23
        double precision, parameter :: kb = 1.380649d-23


        ! Conversion factors between reduced and real units
        ru_time = sqrt(mass*(lj_sigma*1d-10)**2_i64 / (lj_epsilon*1d6))*1d12    ! time in picoseconds
        ru_dens = 1d24 * mass / (Na*lj_sigma**3_i64)                            ! density in g/mL
        ru_temp = 1d3*lj_epsilon/(kb*Na)                                        ! temperature in K


        T = T / ru_temp
        density = density / ru_dens
        dt = dt / ru_time

    end subroutine changeIUnits

   
    subroutine initializePositions(M,a,r,cell)
        ! Chooses and runs the specified initial position.
        ! For the given cell, selects the unit cell matrix and initializes the positions.
        !
        ! Args:
        !   M       (INT64)       : Number of cells in one direction.
        !   a       (REAL64)      : Lattice parameter of the unit cell (L/M)
        !   cell    (CHARACTER)   : String with the name of the cell (sc, bcc, fcc).
        !
        ! Inout:
        !   r       (REAL64[3,N]) : 3xN Positions matrix.  

        implicit none
        integer(kind=i64), intent(in) :: M
        double precision, intent(in) :: a
        character(*), intent(in) :: cell
        double precision, intent(inout),dimension(:,:) :: r
        double precision, allocatable, dimension(:,:) :: ucell

        if ((index(cell,"fcc")==1) .OR. (index(cell,"FCC")==1))  then
            allocate(ucell(3,4))
            ucell = reshape( (/0.d0,0.d0,0.d0, 0.d0,0.5d0,0.5d0,  0.5d0,0.d0,0.5d0,  0.5d0,0.5d0,0.d0/),shape(ucell), order=(/1,2/))
        else if((index(cell,"sc")==1) .OR. (index(cell,"SC")==1)) then
            allocate(ucell(3,1))
            ucell = reshape( (/0.d0,0.d0,0.d0/),shape(ucell), order=(/1,2/))
        else if ((index(cell,"bcc")==1) .OR. (index(cell,"BCC")==1))  then
            allocate(ucell(3,2))
            ucell = reshape( (/0.d0,0.d0,0.d0, 0.5d0,0.5d0,0.5d0/),shape(ucell), order=(/1,2/))
        else
            print*, "Select a valid initial position: sc, bcc, fcc."
        end if

        call initializeGeneral(M,a,ucell,r)     ! Creaties position matrix with the selected unit cell

        deallocate(ucell)

    end subroutine initializePositions


    subroutine initializeVelocities(T,v,init_vel)
        ! Chooses and runs the specified initial velocities.
        !
        ! Args:
        !   T           (REAL64)      : Temperature.
        !   init_vel    (CHARACTER)   : String with the initialization method (bimodal, zero).
        !
        ! Inout:
        !   v           (REAL64[3,N]) : 3xN Velocity matrix.

        double precision,intent(in) :: T
        double precision,intent(inout),dimension(:,:) :: v
        character(*), intent(in) :: init_vel

        if ((index(init_vel,"bim")==1) .OR. (index(init_vel,"BIM")==1))  then
            call initBimodal(T,v)
        else if((index(init_vel,"zero")==1) .OR. (index(init_vel,"ZERO")==1)) then
            v = 0d0
        else
            print*, "Select a valid initial velocity: bimodal, zero."
        end if

    end subroutine initializeVelocities



    subroutine initializeGeneral(M,a,ucell,r)
        ! Initializes the position matrix for a given unit cell.
        !
        ! Args:
        !   M       (INT64)       : Number of cells in one direction.
        !   a       (REAL64)      : Lattice parameter of the unit cell (L/M).
        !   ucell   (REAL[:,:])   : Matrix with the unit cell atom's positions.
        !
        ! Inout:
        !   r       (REAL64[3,N]) : 3xN Positions matrix.   
            
        implicit none
        integer(kind=i64), intent(in) :: M
        double precision, intent(in) :: a
        double precision, dimension(:,:), intent(in) :: ucell
        double precision, intent(inout), dimension(:,:) :: r
        integer(kind=i64) :: i,j,k, p,at

        p = 1_i64
        do i=0,M-1_i64
            do j=0,M-1_i64
                do k=0,M-1_i64
                    do at=1,size(ucell, dim=2,kind=i64)
                        r(:, p) = ucell(:,at) + (/real(i,kind=dp),real(j,kind=dp),real(k,kind=dp)/)
                        !p = at + (k+(j+i*M)*M)*size(ucell, dim=2, kind=i64) 
                        p = p + 1_i64 
                    end do
                end do
            end do
        end do
        r = r*a
        
    end subroutine initializeGeneral


    subroutine initBimodal(T,v)
        ! Initial velocities as Bimodal distribution.
        !
        ! Args:
        !   T           (REAL64)      : Temperature.
        !
        ! Inout:
        !   v           (REAL64[3,N]) : 3xN Velocity matrix.

        implicit none
        double precision,intent(in) :: T
        double precision,intent(inout),dimension(:,:) :: v
        double precision :: vi
        integer(kind=i64) :: N

        N = size(v, dim=2, kind=i64)
        v = 0d0
        if (mod(N,2_i64)/=0_i64) then
            print*, "Number of particles (N) should be multiple of 2."
        end if
        vi = sqrt(T)
        v(:,1:N:2) = -vi
        v(:,2:N:2) = +vi

    end subroutine initBimodal


end module initialization