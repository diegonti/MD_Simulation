# Introduction
This is file which contains relevant information for the workflow of the following days. For the moment,
the contents of this README will be located in this file, and if it is of interest, they can be moved
to the repository README. Therefore, this file is intended to be **readed only by the developers of this repository**.

## Git workflow
Before defining the Git forkflow, it is convenient to set up the machinery:
Start by forking and cloning the [main repository](https://github.com/Eines-Informatiques-Avancades/Project_I).
Once done, set up the following remotes:
```bash
git remote add lider https://github.com/Eines-Informatiques-Avancades/Project_I.git
git remote add origin https://github.com/_your_forked_repository.git
git fetch lider
git pull lider master
```

This has to be set up only the first time.

Now, we will follow to describe the git workflow:

1. `git fetch lider && git pull lider main`
2. Do your changes. It is important that, for each relevant change that you do, do a commit, in the following way: `git add the_file_you_modified && git commit -m "The changes you did"`. This helps 
to keep track of the changes you are doing, if it is necessary to revert in some moment.
2. Once you feel confortable with the changes, **AND THE CODE COMPILES WITH THE DEFINED MAKEFILE** do a push into your github repository: `git push origin master`.
3. Once done, open a pull request into the lider's repository, and everything is done :)

_Note from the lider_: I'll try to schedule myself and accept the pull request every monday at 12:00am.

## Code
It is important to follow a familiar set of rules in order to be consistent with the code, and therefore, avoid some painful bugs:

- The matrices regarding positions and velocities will be structured as 3xnº particles. That is, each column will represent a vector of positions or a vector of velocities for each particle. This is not arbritrary, this is becuse Fortran is used to work with columns.
- To avoid extremely large function definitions, we will encapsulate all parameters in a Fortran derived type, such as:
    ```Fortran
    type :: parameters
        integer(kind=i64) :: n_steps
        .
        .
        .

    end type paramenters
    ```
    And will be defined when reading the input parameter file.
- We will work in double precission. It can be either definded by the following ways:
    ```Fortran
    use, intrinsic :: iso_fortran_env, only: DP => REAL64
    
    real(kind=DP)    :: manera_uno
    double precision :: manera_dos
    real*8           :: manera_tres 
    ```

## Task Distribution
- Main Program
- Readers and Writers (I/O) &rarr; Read input file, write outputs
- Initial Conditions &rarr; Positions (SC, BCC, FCC) and Velocities (bimodal, zero)
- Interactions &rarr; Forces, Energies (total, V, T), Pressure, Momentum
- Integrators &rarr; Velocity Verlet, Verlet, Euler // Thermostat
- Simulation &rarr; PBC, MIC, boxmuller, RDF, MSD
- Data Analysis (Python) &rarr; Stats and Visualization
