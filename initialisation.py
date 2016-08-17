import numpy as np
import math
from modeSolver_STEP import *
#---------------------------------------------------------
#Initialisation
#---------------------------------------------------------
lamda = 0.6361                  #wavelength
a = 14.07                       #core_radius 
n1 = 1.45205                      #core_refractive_index
n2 = 1.44681                       #cladding_refractive_index
l,m = modes(1,4)                #modes define

#---------------------------------------------------------
#waveguide parameters
#---------------------------------------------------------
v,u,w = wave_parameters(a,lamda,n1,n2)
n = v/np.pi - l/2               #number of modes
n = round(n)
z = find_asymptotes(l,n)

#---------------------------------------------------------
#solve for charateristic equation
#---------------------------------------------------------
core = core_char_eqn(u,l)
clad = clad_char_eqn(w,l)
char = char_eqn(u,w,l)

#---------------------------------------------------------
#find roots
#---------------------------------------------------------
solved_u,new_w = find_roots(l,v,w,z,n)

#-------------------------------------------------------
#Mode plots
#------------------------------------------------------- 
radius_range = np.arange(0, 2*a ,0.01)
theta = np.arange(0, 2*np.pi ,0.01)
#theta = np.linspace(0, 2 * np.pi ,1000)
R, Theta = np.meshgrid(radius_range,theta)
X,Y = R*np.cos(Theta),R*np.sin(Theta)

Z1 = mode_distribution_polar(a,l,m,solved_u,new_w)
Z2 = mode_distribution_cartesian(a,l,m,solved_u,new_w,R,Theta)
