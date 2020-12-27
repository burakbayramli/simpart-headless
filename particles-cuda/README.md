# Particles Cuda

This example comes from CUDA development kit `samples/5_Simulations`
folder. It simulates 16384 particles through the GPU. Running and
compiling it requires CUDA + GPU.

I took out as much OpenGL dependencies as possible, but there are
still some left.

Contributions are welcome.

Sample: particles
Minimum spec: SM 3.0

This sample uses CUDA to simulate and visualize a large set of particles and their physical interaction.  Adding "-particles=<N>" to the command line will allow users to set # of particles for simulation.  This example implements a uniform grid data structure using either atomic operations or a fast radix sort from the Thrust library

Key concepts:
Graphics Interop
Data Parallel Algorithms
Physically-Based Simulation
Performance Strategies
