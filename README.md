# SPH Fluid Simulator (C++ / OpenGL)
**An interactive fluid program based on Smoothed Particle Hydrodynamics.**

![Github Stars](https://img.shields.io/github/stars/lijenicol/SPH-Fluid-Simulator?style=flat)
![GitHub forks](https://img.shields.io/github/forks/lijenicol/SPH-Fluid-Simulator?style=flat)
![GitHub issues](https://img.shields.io/github/issues/lijenicol/SPH-Fluid-Simulator?style=flat)
![GitHub](https://img.shields.io/github/license/lijenicol/SPH-Fluid-Simulator)

Smoothed Particle Hydrodynamics (SPH) is a method that can be used to approximate solutions to the Navier-Stokes equations which describe the motion of fluid substances. It is a powerful technique which has its applications in fields such as fluid dynamics and astrophysics. This is my C++/OpenGL implementation of it.

This project started out as a college project, however, I have since added to it so that people can use it as a reference for their own projects. I found that there wasn't too much documentation on implementing this in C++, which caused much of a headache when I couldn't get it working the day the project was due. I hope in making all of this code accessible on GitHub, that it will be able to help someone out who is in a similar situation.

## Video demo:

[![SPH Demo](https://img.youtube.com/vi/FRoIgCHV93U/0.jpg)](https://www.youtube.com/watch?v=FRoIgCHV93U)

## Referencing The Code:

The source code for all of the interesting SPH maths/physics can be found in `src/SPHSystem.cpp`. That file includes all of the code which calculates the densities/forces/positions of particles and that is where the main `update()` method is.

## Neighborhood Search

A core part of the SPH algorithm is calculating properties of particles based off properties of other close particles. This means that we need to determine which particles are "close" to another, which can be computationally expensive. 

The naive solution is the "brute force" approach, which is: `for every particle, loop through every other particle, and if both particles are within a certain distance, mark them as "close"`. This is a classic case of an `O(n^2)` algorithm, an algorithm which doesn't scale well when more particles are added.

Because the distance at which we mark particles "close" is a constant (called `h`), we can use another technique called **Spatial Hashing**. This technique divides the 3D environment into a grid, where the length of each cubic cell is equal to `h`. For each particle, we determine the cell that it is in (on the grid we created) and then hash the cell - we use that hash to place the particle inside a "particle table" whose key is the hash (to handle hash collisions, there is a linked list for each table value). Once the particle table is built, we can then use it to determine which particles are "close" to another. For more info, check out the articles listed in the section **More information on SPH**.

## Multithreading

If you have seen my video on YouTube, you may notice that this program was rather slow - this was due to the fact that it was running only one thread. Since the release of the video, I have improved the performance dramatically by allowing this program to run on more than one thread, which means that the forces/positions of particles can be calculated concurrently, reducing the execution time of the `update()` method. To further increase performance, this is something that a GPU would crush as GPUs have a lot more cores which can increase the amount of calculations which can be done in parallel.

## More Information on SPH ##

These two papers helped out a lot in my journey:

- **SPH Fluids in Computer Graphics**, *State of The Art Report* <br>https://cg.informatik.uni-freiburg.de/publications/2014_EG_SPH_STAR.pdf
- **Particle-Based Fluid Simulation for Interactive Applications** <br>https://matthias-research.github.io/pages/publications/sca03.pdf

## Dependencies:

To build the `sph` executable, there are a few dependencies which are required. These dependencies are:

- GLUT / freeglut
- GLEW
- GLM
- CMake >= 3.13.0

## Building:

This project relies on CMake to generate build files. Once the dependencies are installed, running CMake at the root of this repository will generate the build files necessary to craft the `sph` executable.

Example of building in Linux:

```
cmake -Bbuild
cd build
make
```

Once built, simply run the command below to start the program:

```
./sph
```