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

## Parallelism

The moment propagation is parallelized. It uses the OPENMP API. It is enabled by default (see `-fopenmp` in the Makefile).  
To disable openmp parallelism, remove `-fopenmp` from line 13 of Makefile.

By default, laboetie will use all the threads of your computer.

How do you control the number of threads used by OPENMP?  
You should export OMP_NUM_THREADS=3 in your terminal before executing laboetie if you want laboetie to use 3 threads:  
Thus, to compile and execute laboetie limited to 8 threads on my computer, I type in my terminal:
```bash
$ make
$ export OMP_NUM_THREADS=8
$ ./laboetie
```

Please note that the run time, user time, system time and cpu time printed to you by laboetie are not usable anymore when OMP_NUM_THREADS > 1. This is an issue that needs to be corrected.

Parallelism is implemented by dividing the system along almost independent slices along the z direction: if your system has {100,100,1} nodes along the x, y and z directions, respectively, then it is absolutely useless to multithread your job: use `export OMP_NUM_THREADS=1`.

For now, OPENMP in laboetie is memory bound. Preliminary speed-ups for systems of few hundred of nodes are as follow:
```
# number of threads, speed-up
1, 1
2, 1.9
3, 2.6
4, 3.1
6, 2.7
8, 2.7
```

I would recommand to use 2 or 4 threads only.

## Inputs (lb.in)

* `lx` Number of nodes in x direction
* `ly` Number of nodes in y direction
* `lz` Number of nodes in z direction
* `lbmodel` Lattice Boltzmann geometric model for velocities, e.g., D3Q19
* `timestepmax_for_PoissonNernstPlanck` Maximum number of timesteps in trying to find the equilibrium distribution of charged solutes. Should be 1 salt-free systems.
* `D_equil` number of steps for equilibrating charges (finding PB solution). 1 if charge (salt) free fluid.
* `t_equil` number of steps for equilibrating the flux without constraints
* `tmom` number of steps for equilibrating the flux without constraints
* `tmax` between tmom and tmax, moment propagation is done
* `geometryLabel`
                  # 0 for a custom cell (written in geom.in)
                  # 1 for a slit, i.e. two walls at z=zmin and z=zmax
                  # 2 for a cylinder along Z. lx have to be equal to ly.
                  # 3 for a body centered cubic cell with solid spheres in contact
                  # 4 for disc benichou with exists at 0, 3, 6 and 9 oclock
                  # 5 for corrugated wall (sinusoidal and mirored wrt to y=ly/2.
                  # 6 for a slit with various diameter (2D)
                  # 7 for a tube with various diameter (3D)
                  # 8 for a spherical cavity (3D)
                  # 9 for sphere with 6 exits as asked by benichou
                  # 10 for the supercell given by Xudong, Vincent, Marie and Benjamin from COMSOL
* `stripes` ?
* `initialSolventDensity` fluid density in LB units
* `f_ext` external force
* `charge_distrib` sol) charge distributed in the whole solid. int) charge distributed on interfacial nodes only
* `sigma` = 0.0 # charge distributed in solid
* `bjl` = 0.4 # bjerum length
* `lambda_D` = -2.0 # old debye_l
* `D_plus` = 0.05
* `D_minus` = 0.05
* `D_iter` = 1
* `tracer_Db` = 0.01   # bulk diffusion coefficient of the tracer
* `tracer_Ds` = 0.0   # surface diffusion coefficient of the tracer
* `tracer_z` = 0.0     # charge of the tracer
* `tracer_ka` = 0.1    # adsorption coefficient of the tracer
* `tracer_kd` = 0.01    # desorption coefficient of the tracer
* `c_tracer_0` = 0.1
* `left wall vel` = 0.0
* `right wall vel` = 0.0
* `Tprint_eq` = 1000
* `Tprint_run` = 10000
* `elec_slope` = 0.0 0.0 0.0 # external electric field
* `lncb_slope` = 0.0 0.0 0.0 # external gradient of salt concentration
* `f_gen` = 0.0 0.0 0.0



## Outputs

All outputs files are found in `output/`.

### supercell.xsf

A 3-dimensional representation of the supercell in [xsf format](http://www.xcrysden.org/doc/XSF.html).
*supercelf.xsf* can be opened with [VMD](http://www.ks.uiuc.edu/Research/vmd/): ```vmd -xsf output/supercell.xsf```.
The color code is:
* `pink` Solid nodes
* `green` Interfacial fluid nodes, i.e., fluid nodes close to a solid node
* `white` Non-interfacial fluid nodes
