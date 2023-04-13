#!/bin/bash
set -euo pipefail

read -rp "Minimum cores: " min_core
read -rp "Maximum cores: " max_core
read -rp "Step:          " step_core

if [[ ! -f "MDEMI.x" ]]
then
    echo "No executable has been found, submitting the compiling script. Please, relaunch once the compilation is done"
    # qsub compile_paralel.sh
    exit 1
fi

base_dir="$(pwd)"

if [[ ! -d "speedup" ]]
then
    mkdir speedup
fi

cd speedup

for ((core_num=min_core; core_num<=max_core; core_num=core_num+step_core))
do
    echo "$core_num"
    mkdir sp_${core_num}
    cd sp_${core_num}

    working_dir="$(pwd)"

# Creating the parameter file
    cat << EOF > parameters.nml
&PARAMS
! Particle related variables
  lj_epsilon = 1.77,   ! [kJ/mol]
  lj_sigma = 4.10,     ! [Ang]
  mass = 131.29,       ! [g/mol]

  ! Simulation related variables
  n_steps = 50000                ! Number of simulation steps
  timestep = 0.025,             ! Timestep [ps]
  n_particles = 64,             ! Number of particles
  density = 0.1,                ! Density [g/mL]
  cell_type = "sc",             ! Initial positions cell (sc,bcc,fcc)
  init_velocities = "bimodal",  ! Initial velocities (zero, bimodal)
                                ! The number of particles must follow the formula n_particles = f*M^3
                                ! where f is the number of atoms in the selected cell and M an integer.
                                ! f(sc)=1, f(bcc)=2, f(fcc)=4
  cutoff = 0.4                  ! Cut-off for the energy and force potential truncation
  vlcutoff = 0.6                 ! Cut-off to construct the verlet list buffer zone

  ! I/O variables
  sim_name = "simulation",    ! Simulation name
  write_stats = 10,           ! Steps between each stats writing
  write_frame = 100,          ! Steps between each trajectory writing
  write_file = 100,

  ! Thermostat variables
  andersen_nu = 0.7,        ! Thermostat probability
  temperature = 300,        ! Temperature [K]

  ! Analysis dependent variables
  gdr_num_bins = 100        ! Number of bins for RDF calculation
  n_sweeps = 100	    ! Number of sweeps to compute for the MSD

/
EOF

# Creating the submission file
    cat << EOF > submit_${core_num}.sh
#!/bin/bash

######################
# SGE options and parameters
######################
# (1) Name of the jobs
#$ -N speedup_${core_num}

# (2) Number of cores
#$ -pe smp ${core_num}

# (3) Queue
#$ -q iqtc07.q

# (4) Shell
#$ -S /bin/bash

# (5) Output files
#$ -cwd
#$ -o speedup_c${core_num}.out
#$ -e speedup_c${core_num}.err

module load openmpi/3.1.3_gcc-4.8.5

echo "Execution dir: \$(pwd)"
echo "Starting computation at \$(date)"

mpirun -np ${core_num} ${base_dir}/MDEMI.x ${working_dir}/parameters.nml

echo "Computation done at \$(date)"

exit 0

EOF

    #qsub submit_${core_num}.sh
    cd ..
    sleep 3
done

exit 0