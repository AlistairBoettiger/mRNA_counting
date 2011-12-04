#~/Documents/Berkeley/Levine_Lab/Projects/Enhancer_Modeling/Code/ -tt


# Adapted from Matlab comp_hbdata_binomial.m
# General changes: 	Python indexes from zero, not 1
# 			loops begin with: and end with indentation change, 

""" Simulate Python data Help scripts.
Call as python sim_hbSD.py & to run 
Solves for the expected distribution of fraction of activated nuclei among embryos.  
"""

import sys 	# used to count system inputs, 
import scipy as sp  # needed for sp.special.gammainc function
import math    # has .sqrt, .pi etc
import numpy   # needed for .linspace etc
import matplotlib.pyplot as plt # needed for plotting options
from scipy.special import gamma 
from numpy.random import binomial


from numpy import * # This lets us call numpy functions natively without having to type numpy.
#from pylab import *

# Define a main() function that prints a little greeting.
def main():
  # Get the name from the command line, using 'World' as a fallback.
  if len(sys.argv) >= 2:
    input = sys.argv[1]
  else:
    V = 200 # number of concentration points to check
    time = 10  # min in cell cylce
    Ts = time*60/20 # number of time points to check
    D = 4.5 # diffusion rate of bcd (according to Dostatni)
    cG = 4.8 # c Gregor
    N = 50

    C = numpy.linspace(2,7,V)  # range of concentrations to explore; 
    Th = numpy.linspace(1,1000,1000) # range of thresholds number of molecules to explore
    low = cG-cG*.1 # 10% less than threshold bcd concentration
    high = cG + cG*.1 #  10% more than threshold bcd concentration
    miss_rate = numpy.zeros(Ts)
    bf = numpy.zeros([Ts,N])  
   
    a_effs = numpy.array([.03, .05, .03 + .05 ])/1.9  # Effective enhancer size
    ploop = [.22, .13, .22*.13]   # Looping Probability
    col = ['c','g','r'] # color for plotting
   
    Ps = numpy.zeros(Ts)
    bfD = {}
    
    for e in range(3):
     	a = a_effs[e] 
	    # Main loop    
        for t in range(Ts): # note python indexes from zero
	    T = t*20; # % 7*60;
	    phalf =   sp.special.gammainc(a*cG*D*T,Th)   # 
	    thresh =  numpy.argmin(abs(phalf - .5)) 
	    pint =  1-sp.special.gammainc(a*C*D*T,Th[thresh]);
	    cut = numpy.argmin( abs(C - low))  # find where lines intersect
	    miss_rate[t] = pint[cut];    
	    
#        fig1 = figure()
#	plot(C,pint)
#	# matplotlib.pyplot.axvline(x=low,color = 'c')
#	hold(True) 
#	plot([low,low],[0,1],'c')
#	plot([high,high],[0,1],'r')
#	draw()
 #       show()
	
            p = 1-(1-ploop[e])*(1-miss_rate[t])
	
	   # print ploop[e]
	   # print miss_rate[t]
	   #  Ps[t] = p
	
	    binodist = numpy.zeros(N)
	    for k in range(N):
	       binodist[k] =  choose(N,k)*math.pow(p,k)*math.pow((1-p),(N-k))
	
	   # figure(2)
	   # plot(binodist)
	   # draw()
	
	    bf[t,:] = binodist # N binomial draws 
    
   #     print Ps 
   #     print miss_rate
        bfD[e] = sum(bf[i,:] for i in range(Ts))
        
    xbf = numpy.linspace(0,1,N)
    fig2 = plt.figure(4) 
    plt.plot(xbf,bfD[0],'b-',xbf,bfD[1],'g-',xbf,bfD[2],'r-')
    plt.draw()
    plt.title(['simulation data, N=', N,  'T= ',T]);
    plt.show()

	
def choose(n, k):
  ntok = 1
  for t in xrange(min(k, n-k)):
     ntok = ntok*(n-t)//(t+1)
  return ntok
	

# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
  main()
