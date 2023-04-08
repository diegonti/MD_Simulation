# Project I : MD Simulation of a VdW Gas

This is the official reposotory for the Advanced Computational Tools master subject final project, from the Atomistic and Multiscale Computational Modelling in Physics, Chemistry and Biochemistry master.

The puropose of this project is to develop a serial and pararel version of a Molecular Dynamics simulator of a Lennard-Jones fluid. Therefore, the serial code can be found at the [serial directory](./serial/), and the paralel code at the [paralel directory](./paralel/)

## Table of contents

- [Notes for developers](#notes-for-developers)
- [Task Distribution](#task-distribution)
- [Contributors](#contributors)

## Notes for developers
This section contains relevant information regarding the workflow for the following days. Therefore, this section is indicated for everybody who wants to contribute to this project.

### Git workflow
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

### Code
It is important to follow a familiar set of rules in order to be consistent with the code, and therefore, avoid some painful bugs:

- We will stick to the Fortran 2018 standard. Therefore, all compilations must include the `-std=f2018` flag
- Matrices regarding positions and velocities will be structured as 3xnº particles. That is, each column will represent a vector of positions or a vector of velocities for each particle. This is not arbitrary, this is because Fortran is used to work with columns.
- We will work in double precission. It can be either definded by the following ways:
    ```Fortran
    use, intrinsic   :: iso_fortran_env, only: DP => REAL64
    
    real(kind=DP)    :: manera_uno
    double precision :: manera_dos
    real*8           :: manera_tres 
    ```
    As well, when defining a numerical value, define its numerical type with the `_DP` or `d` delimitiers at the end of the numerical definition.
#### Code style

The indentation is mandatory in order to force readability. Remember, a piece of code is more often read than written, so try to apply this fact when you are writing your code. 

Variables names must have sense in the context of its application, and don't use cryptic names.

Documentation is also mandatory. Since this is a collaboration project, we have to make sure that our code can be interpreted by other humans. Therefore, we will make use of the [Google documentation style](https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html). The language of implementation must be done in English. It is also important in the documentation to specify the type of the data that is handeled. When matrices are involved, they will be defined as `data_type[x,y,z,...]` where x, y, z are the dimensions along the specified rank.

### Builds
A build is simply a compilation of the whole program. We will work on two builds:

- **Debug**: Focused on the program to work, and with the minimum number of bugs as possible. Therefore, we will make use of debugging compiler options, and we will try to solve those warnings as much as possible.
- **Release**: Focused on speed, we will use this build to produce results. Note that this code will start once we validate or finish the Debugging release.

Recall that some smart people in an office probably located in California has spent his/her time working out on smart checking that trigger warnings when parts of our code doesn't make sense or are susceptible to doesn't work properly. Therefore, the purpose of a warning is not to flood your screen making you feel like you're hacking the CNI, but to give you some smart advises in order to avoid bugs. Of course, there are warnings, like -Wtabs, that aren't quite relevant, but others can save us from hours of painful debugging.

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
