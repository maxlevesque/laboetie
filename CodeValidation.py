from matplotlib.pyplot import figure, show
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
from numpy import arange, sin, pi, linspace, array
import math as m
from scipy.optimize import fsolve
from StringIO import StringIO

# You should change xi_initial_guess, whose value is read graphically from the first plot.

PathLocationInput = './lb.in'
PathLocRho = 'output/c_plus_alongZ.dat'
PathLocV = 'output/mass-flux_profile_along_z.dat'

####################################
# Read Input file and create a tag #
####################################

i=0
dico={}
with open(PathLocationInput) as var:
        for line in var:                           # read all the lines
            ligne=line.split()
            if len(ligne)>2 and  ligne[1]=='=' :   # Here we do not take into account the column with '=' signs
                                                   # We also record only lines that have two columns at least 
                dico[ligne[0]]=0
            if len(ligne)>2 and '#' not in line:   # We skip commented lines
                ligne = line.split()
                dico[ligne[0]]=ligne[2:len(ligne)] # we point the variable name to its correponding numerical value
                                                   # to skip the equal sign, ligne starts at 2 instead of 1

#################################
# Finding the Alpha coefficient #
#################################

PI = m.pi                         # PI constant
ENNE = float(dico['lz'][0])-2;    # channel width (only liquid nodes)
Nx = float(dico['lx'][0])         # Box Dimension in x-direction
Ny = float(dico['ly'][0])         # Box Dimension in y-direction
Nz = float(dico['lz'][0])         # Box Dimension in z-direction
bjl = float(dico['bjl'][0]);      # bjerrum length
sigma = float(dico['sigma'][0]);  # charge distribution


if sigma<0.1 :
    Alpha = 2/ENNE * ( PI*m.fabs(sigma/(2*Nx*Ny))*ENNE*bjl )**(0.5)    
else:
    AAA = PI*ENNE*bjl*m.fabs(sigma/(2*Nx*Ny))
    # Equation we wish to solve
    func = lambda xi : AAA - xi* np.tan(xi) 

    # We plot the equation in order to find the zeroes graphically.
    # Choose one zero as initial guess.
    xi = np.linspace(-2.5, 2.5, 201)

    '''
    plt.plot(xi, func(xi))
    plt.xlabel("xi")
    plt.ylabel("func")
    plt.grid()
    plt.show()
    '''


    xi_initial_guess = 1.5706 # this value is read graphically 

    # Use the numerical solver to find the roots
    xi_solution = fsolve(func, xi_initial_guess)

    print "The solution is xi = %f" % xi_solution
    print "at which the value of the expression is = %f" % func(xi_solution)
    Alpha = 2*xi_solution/ENNE

# We now have a value for Alpha
print "Alpha = ", Alpha

############################
# Theoretical Calculations #
############################

# Physical Values and constants

Rho = 1000     # Fluid density kg/m^3 (rho=1000 for water) 
T = 298.15     # Temperature in Kelvin
Kb = 1.38E-23  # Boltzmann constant
Eta = 0.001    # Dynamic viscosity kg/(m.s)
l_b = 0.7E-9   # Bjerrum length (m)
Xi = 2.8837297

# Lattice Boltzmann units of the simulation
f_ext = 0.0#float(dico['f_ext'][0]);
delta_x = l_b/bjl 
delta_t = delta_x**2/(6*Eta/Rho)        

# x - data
Box = linspace(1.5, Nz-0.5, 1000)     # First node (i.e. 1) and last node (i.e. Nz) are solid
                                      # The walls are between the solid node and the fluid node, i.e. 1.5 and
                                      # Nz-0.5
SF = float(Nz)/2 + 0.5                # The theoretical function is centred in zero. We shift the function
                                      # by SF in order to match the numerical results

# y - data
Length = len(Box)
RhoPlus= [None]*Length
v = [None]*Length
for i in range(Length):
    j = Box[i]
    RhoPlus[i] = Alpha**2/( 2*PI*bjl * ( m.cos(Alpha*(j-SF)) )**2 )          # Density profile - Theoretical Sol.
                                                                             # for Poisson-Boltzmann (PB) eq.
    v[i] = m.pow(ENNE,2) * 3./4 *(1-4*(j-SF)*(j-SF)/(m.pow(ENNE,2))) * f_ext # Velocity profile - Poiseuille flow
    #print j

###################
# LB Calculations #
###################

num_linesRho = sum(1 for line in open(PathLocRho))
num_linesV = 142#sum(1 for line in open(PathLocV))
SpaceRho = [None]*num_linesRho
DensityPlus = [None]*num_linesRho
SpaceV = [None]*num_linesV
Jflux = [None]*num_linesV

i=0
with open(PathLocRho) as var:
        for line in var: # read all the lines
                if '#' not in line:
                    line = line.split()           # split all the columns in groups of characters/numbers
                    SpaceRho[i] = float(line[0])  # reads the first column i.e. 0 since we count from zero
                    DensityPlus[i] = float(line[1])
                    i +=1
i=0
with open(PathLocV) as var:
        for line in var:
                line = line.strip()
                if not line or line.startswith('#'): # line is either blank or commented,
                    continue                         # thus ignore it
                    
                line = line.split()        
                SpaceV[i] = float(line[0])  
                Jflux[i] = float(line[2])
                i +=1
                    
#####################
##    Plots        ##
#####################

# If the pressure gradient (i.e. f_ext) is equal to zero, the Poiseuille Flow is not plot
if f_ext!=0.0:
    plot((Box),v,'r-', label='Theoretical')
    plot(SpaceV[int(Nz):],Jflux[int(Nz):],'ko',label='Numerical')
    legend(loc='upper center')
    plt.xlabel('$z-direction$', fontsize = 18)
    plt.ylabel('$Flux$ $Profile$ - $\\rho$ $v(z)$', fontsize = 18)
    plt.xlim(0,Nz+1)
    plt.grid()
    show()
    
#if sigma!=0.0:
plot((Box),RhoPlus,'r-', label='Theoretical')         # Theoretical 
plot(SpaceRho,DensityPlus,'ko', label = 'Numerical')  # LB calculations
legend(loc='lower center')
plt.xlabel('$z-direction$', fontsize = 18)
plt.ylabel('$\\rho_+(z)$', fontsize = 18)
plt.xlim(0,Nz+1)
plt.grid()
plt.savefig('DensityProfile.eps', format='eps', dpi=1000)
show()

#####################
##    For info     ##
#####################

# Physical units (PU)
sigmaPU = sigma/(2*Nz*Nx) * 1/(delta_x**2) # Charge distribution
LPU = Nz*delta_x                           # Channel width

print '----------------------------------------------------------'
print 'your channel width is L = %e m' % LPU
print 'your charge distribution is sigma = %e m^(-2)' % sigmaPU
