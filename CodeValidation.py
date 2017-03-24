from matplotlib.pyplot import figure, show
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
from numpy import arange, sin, pi, linspace, array
import math as m
from scipy.optimize import fsolve
from StringIO import StringIO
from pathlib import Path


# You should change xi_initial_guess, whose value is read graphically from the first plot.

PathLocationInput = './lb.in'
#PathLocRho = './output/c_plus_alongZ.dat'
PathLocRho = './output/c_plus_alongZCHARGESINIT.dat'
PathLocV = 'output/mass-flux_profile_along_z.dat'
#PathLocF = 'output/soluteForceEqZ.dat'
PathLocF = 'output/SFZ.dat'
#PathLocPHI = 'output/phi.dat'
PathLocPHI = 'output/PHI_PNP.dat'


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

# Let us check the geometry of our simulation. If it is not a slit pore skip all the rest                                                   
# 1 => slit pore; 0 => geom.in; 11 => geom.pbm; -1 => bulk
geometryLabel = float(dico['geometryLabel'][0]);
if geometryLabel!=1:
    print 'You are not modelling a slit. Analytical solutions cannot be compared.'
    print 'Change geometryLabel in your input file and set it equal to 1.'
else:
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
    elec_slopex = float(dico['elec_slope'][0]) # E-field in x-direction
    elec_slopey = float(dico['elec_slope'][1]) # E-field in y-direction
    elec_slopez = float(dico['elec_slope'][2]) # E-field in z-direction

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
    # Use the numerical solver to find the roots
    xi_initial_guess = 1.5 # this value is read graphically 
    xi_solution = fsolve(func, xi_initial_guess)

    print "The solution is xi = %f" % xi_solution
    print "at which the value of the expression is = %f" % func(xi_solution)
    Alpha = 2*xi_solution/ENNE
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
    Cs = (1./3)**(0.5)               # Speed of sound in LB units
    f_ext = float(dico['f_ext'][1]); # TODO: f_ext should have different indices according to its
                                     # non-zero value
    delta_x = l_b/bjl 
    delta_t = delta_x**2/(6*Eta/Rho)        

    # -------------
    # x - data
    # -------------
    Box = linspace(1.5, Nz-0.5, 1000)     # First node (i.e. 1) and last node (i.e. Nz) are solid
                                          # The walls are between the solid node and the fluid node, i.e. 1.5 and
                                          # Nz-0.5
    SF = float(Nz)/2 + 0.5                # The theoretical function is centred in zero. We shift the function
                                          # by SF in order to match the numerical results
    # -------------
    # y - data
    # -------------
    SurfArea = (Nx*Ny) # Surface area of walls
                       # Laboetie does a SUM(*) where * = c_plus, phi....
                       # Thus we have to multiply our analytical results by
                       # the surface area = Nx*Ny

    lambda_D = float(dico['lambda_D'][0]); # Debye length
    if lambda_D !=0.0:
        # Low Potential Prefactor
        PrefactLP1 = 1./(8*PI*bjl*lambda_D**2) 
        PrefactLP2 = 4*PI*sigma*bjl*lambda_D/(2*Nx*Ny)
        PrefactLP3 = m.sinh(ENNE / ( 2*lambda_D) )
        # EOF prefactor
        EOFpf = elec_slopey*sigma*lambda_D/(Ny*ENNE)

                       
    Length = len(Box)
    RhoPlus= [None]*Length
    RhoPlusLowPot = [None]*Length
    v = [None]*Length
    F = [None]*Length
    phi_a = [None]*Length
    EOFv = [None]*Length
    for i in range(Length):
        j = Box[i]
        #------------------------------------------------------------------------------------------------
        if lambda_D ==0: # no salt added
            RhoPlus[i] = (Alpha**2/( 2*PI*bjl * ( m.cos(Alpha*(j-SF)) )**2 ))* SurfArea  # Density profile - Theoretical Sol.
                                                                                         # for Poisson-Boltzmann (PB) eq.
            # !!!!!!!!!!!!!! Attention !!!!!!!!!!!!!!!!!!!!
            # The solute_force here below has a factor 1./3, i.e. beta = 1/KbT.
            # We need to check that the force is defined as such in the code.
            F[i] = ( 1./3*(- RhoPlus[i] * 2*Alpha*m.tan(Alpha*(j-SF))) )                 # Solute Force
            #phi_a[i] = ( 2*m.log(m.cos(Alpha*(j-SF))) ) * SurfArea                       # Potential Phi
            phi_a[i] = ( 2*m.log(m.cos(Alpha*(j-SF))/(m.cos(Alpha*ENNE/2))) ) * SurfArea                       # Potential Phi
            EOFv[i] =  1./(PI*bjl) * m.log(m.cos(Alpha*(j-SF))/(m.cos(Alpha*ENNE/2))) * elec_slopey * SurfArea #/1.9 # Electro-osmotic-flow
                                                                                                                # No salt is considered
        else: # salt added
            RhoPlusLowPot[i] = PrefactLP1 * ( 1 - PrefactLP2 * m.cosh( (j-SF)/lambda_D )/(PrefactLP3) ) * SurfArea
            EOFv[i] = (EOFpf * ( m.cosh((j-SF)/lambda_D) - m.cosh(ENNE / ( 2*lambda_D) ) ) / m.sinh(ENNE / ( 2*lambda_D) )) * SurfArea
        #------------------------------------------------------------------------------------------------
        if f_ext !=0:
            v[i] = ( m.pow(ENNE,2) * 3./4 *(1-4*(j-SF)*(j-SF)/(m.pow(ENNE,2))) * f_ext ) * SurfArea    # Velocity profile - Poiseuille flow
        


    ###################
    # LB Calculations #
    ###################


    my_file_RHO = Path(PathLocRho)
    if my_file_RHO.is_file(): # Then file exists
        num_linesRho = sum(1 for line in open(PathLocRho))
        SpaceRho = [None]*num_linesRho
        DensityPlus = [None]*num_linesRho
        i=0
        with open(PathLocRho) as var:
                for line in var: # read all the lines
                        if '#' not in line:
                            line = line.split()           # split all the columns in groups of characters/numbers
                            SpaceRho[i] = float(line[0])  # reads the first column i.e. 0 since we count from zero
                            DensityPlus[i] = float(line[1])
                            i +=1
                            
    my_file_V = Path(PathLocV)
    if my_file_V.is_file(): # Then file exists                        
        num_linesV = sum(1 for line in open(PathLocV)) 
        SpaceV = [None]*num_linesV
        Jflux = [None]*num_linesV
        LowMachV = [None]*num_linesV
        i=0
        with open(PathLocV) as var:
                for line in var:
                        line = line.strip()
                        if not line or line.startswith('#'): # line is either blank or commented,
                            continue                         # thus ignore it
                        line = line.split()        
                        SpaceV[i] = float(line[0])  
                        Jflux[i] = float(line[2])            # TODO: the value here is set to 2, but can be different
                                                             # according to which direction f_ext is set to.
                                                             # Make an automatic change!
                        LowMachV[i] = float(Jflux[i]/Cs)
                        i +=1

    my_file_F = Path(PathLocF)
    if my_file_F.is_file(): # Then file exists
        num_linesF = sum(1 for line in open(PathLocF))
        SpaceF = [None]*num_linesF
        SF = [None]*num_linesF
        i=0
        with open(PathLocF) as var:
                for line in var:
                        line = line.strip()
                        if not line or line.startswith('#'): # line is either blank or commented,
                            continue                         # thus ignore it
                        line = line.split()        
                        SpaceF[i] = float(line[0]) # Spacing
                        SF[i] = float(line[1])     # Solute Force
                        i +=1

    my_file_PHI = Path(PathLocPHI)
    if my_file_PHI.is_file(): # Then file exists
        num_linesPHI = sum(1 for line in open(PathLocPHI))
        SpacePHI = [None]*num_linesPHI
        PHI = [None]*num_linesPHI
        i=0
        with open(PathLocPHI) as var:
                for line in var:
                        line = line.strip()
                        if not line or line.startswith('#'): # line is either blank or commented,
                            continue                         # thus ignore it    
                        line = line.split()        
                        SpacePHI[i] = float(line[0]) # Spacing
                        PHI[i] = float(line[1])      # Potential Phi
                        i +=1
        # When phi is computed a constant is added to the integrand. Therefore the function needs to be
        # shifted by this same constant, which happens to be equal to cst = PHI[int(Nz/2)]
        PHIs = np.array(PHI) #- PHI[int(Nz/2)] # PHI shifted


                        
    #####################
    ##    Plots        ##
    #####################

    # If the pressure gradient (i.e. f_ext) is equal to zero, the Poiseuille Flow is not plot
    if f_ext!=0.0:
        if my_file_V.is_file():
            # Poiseuille Flow plot
            plot((Box),v,'r-', label='Theoretical')
            plot(SpaceV[int(Nz):],Jflux[int(Nz):],'ko',label='Numerical')
            legend(loc='lower center')
            plt.xlabel('$z-direction$', fontsize = 18)
            plt.ylabel('$Flux$ $Profile$ - $\\rho$ $v(z)$', fontsize = 18)
            plt.xlim(0,Nz+1)
            plt.grid()
            show()

            #if lambda_D ==0:

    # If the walls are not charged, then none of the following graphs are plot    
    if sigma!=0.0:
        if my_file_RHO.is_file(): # Then file exists
            # RhoPlus plot
            if lambda_D ==0:
                plot((Box),RhoPlus,'r-', label='Theoretical')        # Theoretical
            else:
                plot((Box),RhoPlusLowPot,'r-', label='Theoretical')
            plot(SpaceRho,DensityPlus,'ko', label = 'Numerical') # LB calculations
            legend(loc='lower center')
            plt.xlabel('$z-direction$', fontsize = 18)
            if lambda_D ==0:
                plt.ylabel('$\\rho_+(z)$', fontsize = 18)
            else:
                plt.ylabel('$\\rho_+(z) - Low$ $Potential$', fontsize = 18)
            plt.xlim(0,Nz+1)
            plt.grid()
            #plt.savefig('DensityProfile.eps', format='eps', dpi=1000)
            show()
            
        if lambda_D ==0: # i.e. no salt
            if my_file_F.is_file(): # Then file exists
                # Solute Force Plot
                plot((Box),F,'r-', label='Theoretical')             # Theoretical 
                plot(SpaceF,SF,'ko', label = 'Numerical')           # LB calculations
                legend(loc='lower center')
                plt.xlabel('$z-direction$', fontsize = 18)
                plt.ylabel('$Solute$ $Force$', fontsize = 18)
                plt.xlim(0,Nz+1)
                plt.grid()
                show()
            if my_file_PHI.is_file(): # Then file exists
                # Potential PHI plot
                plot((Box),phi_a,'r-', label='Theoretical')         # Theoretical 
                plot(SpacePHI, PHIs, 'ko', label = 'Numerical')     # LB calculations
                legend(loc='lower center')
                plt.xlabel('$z-direction$', fontsize = 18)
                plt.ylabel('$Potential$ $\phi$', fontsize = 18)
                plt.xlim(0,Nz+1)
                plt.grid()
                show()

    if elec_slopey!=0.0:
        if my_file_V.is_file():
            # EOF - plot - electro-osmotic-flow with no salt
            plot((Box),EOFv,'r-', label='Theoretical')
            plot(SpaceV[int(Nz):],Jflux[int(Nz):],'ko',label='Numerical')
            legend(loc='lower center')
            plt.xlabel('$z-direction$', fontsize = 18)
            plt.ylabel('$EOF - Flux$ $Profile$ - $\\rho$ $v(z)$', fontsize = 18)
            plt.xlim(0,Nz+1)
            plt.grid()
            show()
            
    #####################
    ##    For info     ##
    #####################

    # Physical units (PU)
    sigmaPU = sigma/(2*Nz*Nx) * 1./(delta_x**2) # Charge distribution
    LPU = Nz*delta_x                           # Channel width
    print '----------------------------------------------------------'
    print 'your channel width is L = %e m' % LPU
    print 'your charge distribution is sigma = %e m^(-2)' % sigmaPU

    # Check Low Mach Number Condition
    '''
    if my_file_F.is_file(): # Then file exists
        LAST = len(F)       # solid node
        FN1 = 2             # first fluid node
        FN2 = LAST-1        # last fluid node
        Fmax1 = F[FN1]      # Force at FN1
        Fmax2 = F[FN2]      # Force at FN2
        LowMachCond1 = (Fmax1*delta_t/2)/Cs
        LowMachCond2 = (Fmax1*delta_t/2)/Cs
        print '--------------- Low Mach Number Condition ----------------'
        print '---------------  In Solute Force profile  ----------------'
        print ' The Condition is the following : Ratio << 1'
        print 'At first fluid node your LMN Ratio = %e m' % LowMachCond1
        print 'At last fluid node your LMN Ratio = %e m'  % LowMachCond2

    # TODO: Change the following checks
        if LowMachCond1>1.e-1 or LowMachCond2>1.e-1:
            print 'The low Mach Number condition is not fullfilled'
            print ' Check output/soluteForceEqZ*.dat'
    '''
    if f_ext!=0.0:
        if my_file_V.is_file():
            if any(t>1.e-1 for t in LowMachV):
                print '--------------- Low Mach Number Condition ----------------'
                print '---------------   In Mass Flux Profile    ----------------'
                print ' The Condition is the following : Ratio << 1'
                print 'The low Mach Number condition is not fullfilled'
                print ' Check output/mass-flux_profile_along_z.dat'

