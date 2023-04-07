module writers_m
    ! Module for writing the data to output files.
    ! Contains subroutines to write the system data and positions.

    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64, error_unit
    implicit none

    contains

    subroutine writeSystem(unit,lj_epsilon,lj_sigma,mass, time,E,Epot,Ekin,T,press,momentum)
        ! Writes system data to the main output file. Changes reduced units 
        ! used in simulation to real units.
        !
        ! Args:
        !   unit        (INT64)  : File unit to write on.
        !   lj_epsilon  (REAL64) : Lennard Jones epsilon parameter for the gas (in kJ/mol).
        !   lj_sigma    (REAL64) : Lennard Jones sigma parameter for the gas (in Ang).
        !   mass        (REAL64) : Molar Mass of the gas (g/mol).
        !   *args       (REAL64) : Simulation params to write. (time,E,Epot,Ekin,T,press,MSD,momentum)
        
        integer(kind=i64), intent(in) :: unit
        double precision, intent(in) :: lj_epsilon,lj_sigma,mass
        double precision, intent(in) :: time,E,Epot,Ekin,T,press,momentum
        double precision :: time_out,E_out,Epot_out,Ekin_out,T_out,press_out,momentum_out
        double precision :: ru_time,ru_dens,ru_dist,ru_temp,ru_E,ru_press,ru_vel,ru_mom
        double precision, parameter :: Na = 6.0221408d23
        double precision, parameter :: kb = 1.380649d-23

        ! Conversion factors between reduced and real units
        ru_time = sqrt(mass*(lj_sigma*1d-10)**2_i64 / (lj_epsilon*1d6))*1d12! time in ps
        ru_dist = lj_sigma                                                  ! distance in Ang
        ru_dens = 1d24 * mass / (Na*lj_sigma**3_i64)                        ! density in g/mL
        ru_E = lj_epsilon                                                   ! energy in kJ/mol
        ru_press = 1d3*lj_epsilon/((lj_sigma*1d-10)**3_i64 * Na)            ! pressure in Pa
        ru_temp = 1d3*lj_epsilon/(kb*Na)                                    ! temperature in K
        ru_vel = ru_dist/ru_time *1d-10                                     ! velocity in Ang/s
        ru_mom = 1d-3*(mass/Na)*ru_vel*1d-10                                ! linear momentum in kg*m/s

        ! Multiplying the reduced value with its correspondent conversion factor
        time_out = time * ru_time
        E_out = E * ru_E
        Epot_out = Epot * ru_E
        Ekin_out = Ekin * ru_E
        T_out = T * ru_temp
        press_out = press * ru_press
        momentum_out = momentum * ru_mom

        ! Writes to output file (coumun style)
        write(unit,'(ES18.8e4,ES18.8e4,ES18.8e4,ES18.8e4,ES18.8e4,ES18.8e4,ES18.8e4)') time_out,&
        E_out,Epot_out,Ekin_out,T_out,press_out,momentum_out

    end subroutine writeSystem


    subroutine writePositions(r,unit)
        ! Writes current position in the specified file (XYZ format).
        ! In a loop, writes trajectory.
        !
        ! Args:
        !   r       (REAL64[3,N]) : 3xN Positions matrix.  
        !   unit    (INT64)       : File unit to write on.
        implicit none
        double precision, dimension(:,:), intent(in) :: r
        integer(kind=i64), intent(in) :: unit
        integer(kind=i64) :: i, N

        N = size(r, dim=2,kind=i64)

        write(unit, '(I3)') N
        write(unit,'(A)') ''
        do i= 1,N
            write(unit,'(A,F20.8,F20.8,F20.8)') "Xe", r(1,i), r(2,i), r(3,i)
        end do


    end subroutine writePositions

    subroutine writeRdf(rdf,unit, ljsigma)
        ! Writes positions and rdf in 2 columns.
        !
        ! Args:
        !   dr      (REAL64)       : step dr
        !   rdf     (REAL64[bins]) : radial distribution function values at 0, dr, 2dr,...
        !   unit    (INT64)        : File unit to write on.
        implicit none
        ! In/Out variables
        double precision, intent(in), dimension(:,:) :: rdf
        real(kind=dp), intent(in)                    :: ljsigma
        integer(kind=i64), intent(in)                :: unit
        ! Internal variables
        integer(kind=i64)                            :: i, N

        N = size(rdf, kind=i64, dim=2)

        do i= 1,N
            write(unit,*) rdf(1,i)*ljsigma, rdf(2, i)
        end do

    end subroutine writeRdf

    subroutine writeMSD(v_MSD,unit,n_sweep,write_log,lj_sigma,lj_epsilon,dt,mass)
        ! Writes positions and rdf in 2 columns.
        !
        ! Args:
        !   v_MSD      (REAL64[N_steps] : Vector containig the MSD values in each dt
        !   unit       (INT64)          : File unit to write in.
        !   n_sweep    (INT64)          : Number of sweeps to compute for the MSD
        !   write_log  (INT64)          : Frequency to write into file
        !   lj_sigma   (REAL64)         : Reduced units of length
        !   lj_epsilon (REAL64)         : Reduced units of energy
        !   dt         (REAL64)         : Timestep
        !   mass       (REAL64)         : mass of particle

        implicit none
        ! In/Out variables
        double precision, intent(in), dimension(:)   :: v_MSD
        double precision, intent(in)                 :: lj_sigma, lj_epsilon, dt, mass
        integer(kind=i64), intent(in)                :: unit, n_sweep, write_log
        ! Internal variables
        integer(kind=i64)                            :: i, N
        double precision                             :: ru_time

        N = size(v_MSD, 1, i64)
        ru_time = 1d12 * dsqrt(mass * (1d-10*lj_sigma)**2 / (1d6*lj_epsilon))

        write(unit,*) 0.0d0, 0.0d0
        do i= 1, N, write_log
            if (i <= N-n_sweep) then
                write(unit,*) dble(i)*dt*ru_time, lj_sigma**2 *v_MSD(i)/dble(n_sweep)
            else
                write(unit,*) dble(i)*dt*ru_time, lj_sigma**2 *v_MSD(i)/dble(N-i+1)
            end if
        end do

    end subroutine writeMSD


end module writers_m
