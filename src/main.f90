program main
    use, intrinsic :: iso_fortran_env, only: DP => real64, I64 => int16, input_unit, output_unit
    ! Module definitions

    implicit none

    ! ~ Memory definition ~
    ! Array variables

    ! Scalar variables
    real(kind=dp) :: init_time, end_time

    !String variables


    call cpu_time(init_time)
    write(output_unit, '(A)') 'Welcome to MDEMI!!'

    call cpu_time(end_time)
    write(output_unit, '(A,F12.8,A)') 'Execution done in: ', end_time - init_time, ' seconds.'
end program main