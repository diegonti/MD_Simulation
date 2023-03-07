module initialization
    ! Module with the initialization of the position and velocity matrices.
    ! Contains two general functions that choose the desired initialization method.
    ! 3 types of initial positions are implemented -> SC, BCC, FCC.
    ! 2 types of initial velocities are implemented -> bimodal, zero.

    use, intrinsic :: iso_fortran_env, only: DP => real64, I64 => int64, input_unit, output_unit
    contains

    subroutine getInitialParams(cell,N,density,M,L,a)
        ! Returns the initial cell parameters of the simulation.
        ! Inputs: 
        ! cell: String with the name of the cell (sc, bcc, fcc).
        ! N : Number of particles.
        ! density : Density of the system.
        ! Outputs:
        ! M : Number of cells in one direction.
        ! L : Length of the simulation box.
        ! a : Lattice parameter of the unit cell (L/M).

        character(*), intent(in) :: cell
        integer(kind=i64), intent(in) :: N
        double precision, intent(in) :: density
        integer(kind=i64), intent(out) :: M
        double precision, intent(out) :: L, a
        

        !# Dejo numeros real*4 porque sino hay error y la M da uno menor. Round?->INT
        L = (real(N,kind=dp)/density)**(1.d0/3.d0)
        if (cell=="sc") then; M = int((real(N,kind=dp))**(1.d0/3.d0),kind=i64)
        else if (cell=="fcc") then; M = int((real(N,kind=dp)/4.d0)**(1.d0/3.d0), kind=i64)
        else if (cell=="bcc") then; M = int((real(N,kind=dp)/2.d0)**(1.d0/3.d0), kind=i64)
        else; print*, "Select a valid initial position: sc, bcc, fcc."
        end if
    
        a = L/real(M, kind=dp)            ! Lattice parameter
        
    end subroutine getInitialParams
    

    subroutine initializePositions(M,a,r,cell)
        ! Chooses and runs the specified initial position.
        ! For the given cell, selects the unit cell matrix and initializes the positions.
        ! M : Number of cells in one direction.
        ! a : Lattice parameter of the unit cell (L/M)
        ! cell: String with the name of the cell (sc, bcc, fcc).
        ! r : 3xN Positions matrix.  

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
        ! T : Temperature.
        ! v : 3xN Velocity matrix.
        ! init_vel: String with the initialization method (bimodal, zero).

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
        ! M : Number of cells in one direction.
        ! a : Lattice parameter of the unit cell (L/M)
        ! ucell: Matrix with the unit cell atom's positions.
        ! r : 3xN Positions matrix.  
            
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
        ! T : Temperature.
        ! v : 3xN Velocities matrix.

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