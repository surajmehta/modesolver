from scipy.special import jv
import matplotlib.pyplot as plt

from initialisation import *

#-------------------------------------------------------
#plot characteristic equation
#-------------------------------------------------------
plt.figure(1)
plt.plot(u,core,label = "core")                     #plot core
plt.plot(u,clad,label = "clad")                     #plot cladding
#plt.plot(u,char,label = "char")                     #plot char = core - cladding
plt.plot(solved_u, core_char_eqn(solved_u,l), 'o')  #plot intersect points
plt.legend()
plt.grid()  
plt.xlabel("u_parameter")
plt.title("characteristic Equation")
plt.xlim(0,round(v))
plt.ylim(0,clad_char_eqn(v,l)+1)

#-------------------------------------------------------
#plot mode distribution
#-------------------------------------------------------
plt.figure(2)
plt.plot(radius_range,Z1)

#-------------------------------------------------------
#plot Intensity
#-------------------------------------------------------
plt.figure(3)
plt.pcolormesh(X, Y, Z2)

#show figures
plt.show()