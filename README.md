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

To compile the program run:
``` 
$ make
```
To run the program:
``` 
$ make run
```
Once the progam has finished succesfully, the output files will be created in the parent directory. Then the statistics and visualization can be generated with:
```
(!subject to change)
$ make stats (-flags)
$ make plots (-flags)
```
Which will create an `name_stats.log` in the current directory and the plots (block average and visualization) in the `plots/` (or specified) folder. More info about the Python scripts in [here](https://github.com/diegonti/Project_I/tree/master/scripts).


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
