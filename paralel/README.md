# Introduction

Here, the instructions to install and run the pararel version of this project will be explained.

## Installing

This project is designed to be compiled on Unix-type machines (Linux/macOS). The autorhs do not guarantee its correct compilation on Windows machines.

Before starting the compilation, make sure to have the following packages available in your computer:
- gfortran >= 11.3.0
- OpenMPI >= 4.1.2

For the statistics and visualization scripts the following Python modules are required: `numpy`, `scypy`, `matplotlib`, `ase`. Which can be installed with:
```
$ pip install -r scripts/requeriments.txt
or $ pip install module    (for each module)
```
Once meeting the requirements, inside this folder, just do `make`. A Makefile script will be interpreted by make, and will compile all the source files.

## Usage

In order to use the simulator, a parameter file must be supplied as a first argument list, such as `./MDEMI.x parameters.nml`. An example of parameter input file can be found at [here](./parameters.nml).

To run the program:
``` 
$ make run_parallel nproc=NP
```
Where `$NP` are the number of processors desired in order to run the program. Once the progam has finished succesfully, the output files will be created in the parent directory. Then the statistics and visualization can be generated with:
```
$ make postprocess input="name_logfile.log"
```
Which will create an `name_stats.log` in the current directory and the plots (block average and visualization) in the `plots/` (or specified) folder. 

The output folder of the plots and different flags, such as makeing a trajectory GIF or selecting the start and final frame to make the plots and statistichs can be selected adding an additional argument (optional):
```
$ make postprocess input="name_logfile.log" args="-op output_file -s start -f finish -t"
```
If wanted, the plots and visualization can be done separately with:
```
$ make stats input="name_logfile.log" args="-op output_file -s start -f finish"
$ make plots input="name_logfile.log" args="-op output_file -s start -f finish -t"
```

More info about the Python scripts and their flags in [here](scripts/README.md).

For general help about the makefile try:
```
$ make help
```
And to clean all the generated compilation files, use:
```
$ make clean
```