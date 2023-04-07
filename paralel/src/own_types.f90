module own_types_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    implicit none

    type :: databloc_params_t
        ! Particle-related variables
        real(kind=dp)       :: mass
        real(kind=dp)       :: lj_epsilon
        real(kind=dp)       :: lj_sigma 
        
        ! Simulation-related variables
        real(kind=dp)       :: timestep      ! Reduced time
        real(kind=dp)       :: cutoff_set
        integer(kind=i64)   :: n_particles   
        real(kind=dp)       :: density       ! Reduced density
        integer(kind=i64)   :: n_steps       ! Equilibration steps  
        
        ! I/O variables
        integer(kind=i64)   :: write_file    ! How often do we write log information
        integer(kind=i64)   :: write_stats   ! How often do we write statistics
        integer(kind=i64)   :: write_frame   ! How often write frame
        character(len=2048) :: sim_name      ! Name of the simulation we are doing
        
        ! Thermostat variables
        real(kind=dp)       :: ref_temp      ! Reference temperature for thermostats and velocity initialization 
        real(kind=dp)       :: andersen_nu   ! Andersen parameter
        
        ! Simulation-dependent variables
        real(kind=dp)       :: box           ! Size of the box
        
        ! Analysis variables
        integer(kind=i64)   :: n_sweeps
        integer(kind=i64)   :: gdr_num_bins
        real(kind=dp)       :: gdr_max_dist
    end type 

end module
