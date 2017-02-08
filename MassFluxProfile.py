#!/usr/bin/python
import string, re, struct, sys, math, os, time
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
from pylab import *


num_lines = sum(1 for line in open('output/mass-flux_profile_along_z.dat'))
Space = [None]*num_lines
Jflux = [None]*num_lines
i=0
PathLocation = 'output/mass-flux_profile_along_z.dat'
with open(PathLocation) as var:
        for line in var: # read all the lines
                #if '#' not in line:
                line = line.split()        # split all the columns in groups of characters/numbers
                Space[i] = float(line[0])  # reads the first column i.e. 0 since we count from zero
                Jflux[i] = float(line[2])
                i +=1

plot(Space,Jflux,'k--')
legend(loc='lower right')
plt.xlabel('$z-direction$', fontsize = 18)
plt.ylabel('$j$-flux', fontsize = 18)
#savefig('LogLogGraph.eps', format='eps', dpi=1000)
show()
