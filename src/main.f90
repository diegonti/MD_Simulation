program main
    use, intrinsic :: iso_fortran_env, only: DP => real64, I64 => int16, input_unit, output_unit
    ! Module definitions
    use initialization
    use testing

    implicit none

    ! ~ Memory definition ~
    ! Array variables
    double precision, allocatable,dimension(:,:) :: r,v

    ! Scalar variables
    real(kind=dp) :: init_time, end_time
    double precision :: density,L,a, T
    integer(kind=i64) :: M, N

    ! String variables
    character(len=:), allocatable :: cell
    character(len=:), allocatable :: init_vel


    !!! ~ MAIN PROGRAM ~ !!!

    call cpu_time(init_time)
    write(output_unit, '(A)') 'Welcome to MDEMI!!'

    ! Test INPUTS
    N = 32_i64          ! Number of particles
    density = 0.8d0     ! Density
    T = 2d0             ! Temperature

    ! The below strings for inicialization should be gathered from the input of the user
    ! if its possible, whatever the user inputs, change it to lower case. (subroutine lowercase @diegonti)
    cell = trim("fcc")
    init_vel = trim("bimodal")

    ! Allocate Memory
    allocate(r(3,N))
    allocate(v(3,N))

    ! Initialization of the system
    call getInitialParams(cell,N,density,M,L,a)
    call initializePositions(M,a,r,cell)
    call initializeVelocities(T,v,init_vel)

    !call testMatrix(r)


    call cpu_time(end_time)
    write(output_unit, '(A,F12.8,A)') 'Execution done in: ', end_time - init_time, ' seconds.'

end program main