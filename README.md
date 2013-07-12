# laboetie

laboetie is an elektrokinetic lattice boltzmann code with moment propagation and electrokinetics.

## Authors

Maximilien Levesque  
Benjamin Rotenberg

## Acknowledgments

laboetie is written on the basis of a code by Capuani, Frenkel, Rotenberg et al.(?).

## Strategy

laboetie is written in modern fortran (2003+).  
You are expected to find allocatable arrays, modules, object oriented programming, do concurrent, among others.  
It should be coded with test-driven development in mind.  

## License

There is a license associated to this program. Read the LICENSE file.

## Github

Github repo originaly created on 2013/02/01.

## Installation instructions

You need to have scons installed. SCons is a modern GNU make with much easier syntax and above all very good dependence tracking.
Check your linux distribution tutorials to install SCons. A simple `sudo make install scons` is sufficient in Ubuntu, or `sudo yum install scons` in Fedora derivatives.

Once installed, just type `scons` in the folder you extracted laboetie in.

