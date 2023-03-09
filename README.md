# Project I : MD Simulation of a VdW Gas


## Table of contents

- [Installing](#installing)
- [Usage](#usage)
- [Task Distribution](#task-distribution)
- [Contributors](#contributors)

## Installing

This project is designed to be compiled on Unix-type machines (Linux/macOS). The autorhs do not guarantee its correct compilation on Windows machines.

Before starting the compilation, make sure to be using the gfortran compiler, at least version `11.3.0` or greater. Other compilers have not been tested.

For the statistics and visualization scripts the following Python modules are required: `numpy`, `scypy`, `matplotlib`, `ase`. Which can be installed with:
```
$ pip install -r scripts/requeriments.txt
or $ pip install module    (for each module)
```

Once meeting the requirements, `git clone` this repository, and, inside the repo folder, just do `make`. A Makefile script will be interpreted by make, and will compile all the source files.

## Usage

In order to use the simulator, a parameter file must be supplied as a first argument list, such as `./MDEMI.x parameters.nml`. An example of parameter input file can be found at [here](./parameters.nml).

To compile the program use:
``` 
$ make
```
To run the program:
``` 
$ make run_serial
```
Once the progam has finished succesfully, the output files will be created in the parent directory. Then the statistics and visualization can be generated with:
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

More info about the Python scripts and their flags in [here](https://github.com/diegonti/Project_I/tree/master/scripts).

For general help about the makefile try:
```
$ make help
```
And to clean all the generated compilation files, use:
```
$ make clean
```

<br>

## Task Distribution
- Main Program &rarr; Marc Alsina
- Readers and Writers (I/O) &rarr; Emma Valdés
- Initial Conditions &rarr; Diego Ontiveros
- Interactions &rarr; Marc Alsina
- Integrators &rarr; Ignacio Prieto & Emma Valdés
- Simulation (PBC & MIC, RDF, SMD) &rarr; Maitane Fariñas
- Data Analysis (Python) &rarr; Diego Ontiveros

<br>

## Contributors

|  Marc Alsina   |  Maitane Fariñas  |  Diego Ontiveros   |  Ignacio Prieto   |  Emma Valdés  |
| -------------- | ----------------- | ------------------ | ----------------- | ------------- |
| ![malsinac](https://github.com/malsinac.png "malsinac") | ![maitanefarinas](https://github.com/maitanefarinas.png "maitanefarinas") | ![diegonti](https://github.com/diegonti.png "diegonti") | ![Ronoh97](https://github.com/Ronoh97.png "Ronoh97") | ![evaldesmartin](https://github.com/evaldesmartin.png "evaldesmartin") |
| [malsinac](https://github.com/malsinac)                                 | [maitanefarinas](https://github.com/maitanefarinas)| [diegonti](https://github.com/diegonti)                                  | [Ronoh97](https://github.com/Ronoh97)                                  | [emma](https://github.com/evaldesmartin)                                  |

For further questions please refer to the main resposible of this project: [malsinac](https://github.com/malsinac) 
