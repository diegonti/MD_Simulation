module readers_m
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64, error_unit
    implicit none

contains

subroutine read_nml(param_file, lj_epsilon, lj_sigma, mass, timestep, density, &
    andersen_nu, n_particles, n_steps, write_file, write_stats, gdr_num_bins, &
    write_frame, sim_name, cell_type, init_velocities)
    implicit none
    ! Author: Marc Alsina <marcalsinac@gmail.com>
    ! Subroutine to read the simulation parameters from a fortran nml list
    !
    ! Args:
    !   param_file (CHARACTER[*]): The path of the parameter file 

    ! In/Out variables
    real(kind=dp), intent(inout)         :: lj_epsilon, lj_sigma, mass, timestep, density, andersen_nu
    integer(kind=i64), intent(inout)     :: n_particles, n_steps, write_file, write_stats, gdr_num_bins, write_frame
    character(len=2048), intent(inout)   :: sim_name, cell_type, init_velocities
    character(len=2048), intent(in)      :: param_file
    ! Internal variables
    integer(kind=I64)                    :: unit_nr, iost
    character(len=1024)                  :: msg
    ! Namelist variables
    

    namelist /PARAMS/ lj_epsilon, lj_sigma, mass, &
    timestep, n_particles, density, n_steps, &
    write_file, write_stats, write_frame, sim_name, &
    andersen_nu, gdr_num_bins, cell_type, init_velocities

    ! Checking if file is present
    inquire(file=trim(param_file), iostat=iost)
    if (iost /= 0_I64) then
        write(unit=error_unit, fmt='(A)') "FATAL ERROR: namelist file not found"
        stop 1
    end if

    ! Reading namelist
    open(newunit=unit_nr, file=trim(param_file), access='sequential', &
    action='read', iostat=iost, iomsg=msg, form='formatted')

    read(iostat=iost, unit=unit_nr, nml=PARAMS, iomsg=msg)
    if (iost /= 0_I64) then
        write(unit=error_unit, fmt='(A,A)') "Error when loading parameters on the namelist file: ", trim(msg)
        stop 1
    end if
    close(unit_nr)

end subroutine read_nml


end module readers_m