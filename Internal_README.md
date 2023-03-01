# Introduction
This file contains relevant information regarding the workflow for the following days. For the moment,
the contents of this README will be located in this file, and if it is of interest, they will be moved
to the repository README. Therefore, this file is intended to be **read only by the developers of this repository**.

## Git workflow
Before defining the Git workflow, it is convenient to set up the machinery:
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

1. `git fetch lider && git pull lider master`
2. Do your changes. It is important that, for each relevant change that you do, do a commit, in the following way: `git add the_file_you_modified` (as many times as files you have modified) `git commit -m "The changes you did"`. This helps 
to keep track of the changes you are doing, if it is necessary to revert in some moment.
2. Once you feel comfortable with the changes, **AND THE CODE COMPILES WITH THE DEFINED MAKEFILE** do a push into your github repository: `git push origin master`.
3. Once done, open a pull request into the lider's repository, and everything is done :)

_Note from the lider_: I'll try to schedule myself and accept the pull request every Monday at 12:00am.

## Code
It is important to follow a familiar set of rules in order to be consistent with the code, and therefore, avoid some painful bugs:

- We will stick to the Fortran 2018 standard. Therefore, all compilations must include the `-std=f2018` flag
- Matrices regarding positions and velocities will be structured as 3xnº particles. That is, each column will represent a vector of positions or a vector of velocities for each particle. This is not arbitrary, this is because Fortran is used to work with columns.
- To avoid extremely large function definitions, we will encapsulate all parameters in a Fortran derived type, such as:
    ```Fortran
    type :: parameters
        integer(kind=i64) :: n_steps
        .
        .
        .

    end type parameters
    ```
    And will be defined when reading the input parameter file.
- We will work in double precission. It can be either definded by the following ways:
    ```Fortran
    use, intrinsic   :: iso_fortran_env, only: DP => REAL64
    
    real(kind=DP)    :: manera_uno
    double precision :: manera_dos
    real*8           :: manera_tres 
    ```
    As well, when defining a numerical value, define its numerical type with the `_DP` or `d` delimitiers at the end of the numerical definition.
### Code style

The indentation is mandatory in order to force readability. Remember, a piece of code is more often read than written, so try to apply this fact when you are writing your code. 

Variables names must have sense in the context of its application, and don't use cryptic names.

Documentation is also mandatory. Since this is a collaboration project, we have to make sure that our code can be interpreted by other humans. Therefore, we will make use of the [Google documentation style](https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html). The language of implementation must be done in English.

## Task Distribution
- Main Program
- Readers and Writers (I/O) &rarr; Read input file, write outputs
- Initial Conditions &rarr; Positions (SC, BCC, FCC) and Velocities (bimodal, zero)
- Interactions &rarr; Forces, Energies (total, V, T), Pressure, Momentum
- Integrators &rarr; Velocity Verlet, Verlet, Euler // Thermostat
- Simulation &rarr; PBC, MIC, boxmuller, RDF, MSD
- Data Analysis (Python) &rarr; Stats and Visualization

## Builds
A build is simply a compilation of the whole program. We will work on two builds:

- **Debug**: Focused on the program to work, and with the minimum number of bugs as possible. Therefore, we will make use of debugging compiler options, and we will try to solve those warnings as much as possible.
- **Release**: Focused on speed, we will use this build to produce results. Note that this code will start once we validate or finish the Debugging release.

Recall that some smart people in an office probably located in California has spent his/her time working out on smart checking that trigger warnings when parts of our code doesn't make sense or are susceptible to doesn't work properly. Therefore, the purpose of a warning is not to flood your screen making you feel like you're hacking the CNI, but to give you some smart advises in order to avoid bugs. Of course, there are warnings, like -Wtabs, that aren't quite relevant, but others can save us from hours of painful debugging.