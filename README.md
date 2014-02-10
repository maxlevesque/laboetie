# laboetie

laboetie is an elektrokinetic lattice boltzmann code with moment propagation and electrokinetics.

## Authors

Maximilien Levesque  
Benjamin Rotenberg

## Acknowledgments

laboetie is written on the basis of a code by Capuani, Frenkel, Rotenberg et al.

## Strategy

laboetie is written in fortran. Some Fortran 2003 or more advanced functions require not-too-old versions of gcc-gfortran.
You are expected to find allocatable arrays, modules, object oriented programming, do concurrent, among others.  

## License

This program, its sources, manual, etc., MUST NOT BE NOR DISTRIBUTED NOR MODIFIED NOR ANYTHING WITHOUT ASKING BENJAMIN ROTENBERG AND MAXIMILIEN LEVESQUE.

## Github

Github repo originaly created on 2013/02/01.

## Installation instructions

You need to have scons installed. SCons is a modern GNU make with much easier syntax and above all very good dependence tracking.
Check your linux distribution tutorials to install SCons. A simple `sudo make install scons` is sufficient in Ubuntu, or `sudo yum install scons` in Fedora derivatives.

Once installed, just type `scons` in the folder you extracted laboetie in.

## Inputs (lb.in)

* `lx` Number of nodes in x direction
* `ly` Number of nodes in y direction
* `lz` Number of nodes in z direction
* `lbmodel` Lattice Boltzmann geometric model for velocities, e.g., D3Q19
* `timestepmax_for_PoissonNernstPlanck` Maximum number of timesteps in trying to find the equilibrium distribution of charged solutes. Should be 1 salt-free systems.
