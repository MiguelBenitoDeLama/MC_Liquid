# Monte Carlo Liquid Simulation program

Program to perform a liquid simulation in order to obtain the radial distribution function using a Monte Carlo approach. The program can use two different potential types, Lennard-Jones and Stillinger potentials. This program is part of the assesment of "Computational Chemistry Programming Project" of the TCCM and it was written following the manual provided by Jeremy Harvey.    

Installation
------------

To compile the program, run `make` from the directory that contains the f90 files. The Makefile is prepared for gfortran, if your compiler is not gfortran modify FC and FFLAGS. The compilation will create an executable "program.exe" which is used to run the program.

Input
-----

The input of the program is introduced using the file "input.in", with the format shown in the input examples that are contained in the "input_examples" directory.

The input contains the parameters of the simulation and representation of the radial distribution function, as well as some options of the program. The potential type of the simulation can be selected between a Lennard Jones potential "LJ" and a Stillinger potential "ST", the program can use or not a neighbour list method in order to speed up the calculation and the initial structure can be a given structure through a ".xyz" file or it can be randomly generated.

Output
------

The program generates three output files:

"MC_simulation.log": Contains the parameters used in the simulation and the results of the simulation (success ratio of the MC simulation and cpu time).

"initial_structure.xyz": Contains the number of particles and initial structure used in the simulation.

"radial.dat": Contains two columns, the distance for which the radial distribution function is calculated and the value of the radial distribution function for each distance.
