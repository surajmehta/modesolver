from scipy.special import jv, kv, jn_zeros
from scipy.optimize import fsolve
import numpy as np

#---------------------------------------------------------
#Utility Functions
#---------------------------------------------------------
def sqrt_calculate(a,b):
    """
    
    """
    return np.sqrt(((a)**2-(b)**2))  

def find_roots(l,v,w,z,n):
    """
    """
    nn = [None]*n
    nn[0]=(0 + z[0])/ 2
    for i in range(1,int(n)):
        if (i == 0):
            continue
        nn[i]= (z[i-1]+z[i])/2
    
    solved_u = fsolve(char_eqn,nn,args=(v,l))
    new_w = sqrt_calculate(v,solved_u)
    return solved_u,new_w

def modes(l,m):
    """
    modes
    """
    if(m < 1):
        return 0
    l = l
    m=m-1
    return l,m

def wave_parameters(a,lamda,n1,n2):
    """
    """
    #Numerical Aperture
    NA = sqrt_calculate(n1,n2)
    
    v = 2 * np.pi * a * NA / lamda      #v = NA * Ko * a
    u = np.arange(0,v,0.01) #u = Kt * a int(round(v))
    w = sqrt_calculate(v,u)             #w = gamma * a    
    return v,u,w

#---------------------------------------------------------
#Characteristic Equations
#---------------------------------------------------------
def core_char_eqn(u,l):
    """
    """
    core_eqn = u * jv(l+1,u) / jv(l,u)
    return core_eqn

def find_asymptotes(l,n):
    """
    """
    z = jn_zeros(l,n)    
    return z

def clad_char_eqn(w,l):
    """
    """
    clad_eqn = w * kv(l+1,w) / kv(l,w)
    return clad_eqn

def char_eqn(u,v,l):
    """
    """
    w = sqrt_calculate(v,u)
    core_eqn = core_char_eqn(u,l)
    clad_eqn = clad_char_eqn(w,l)
    char = core_eqn - clad_eqn
    return char

#---------------------------------------------------------
#Mode Equations
#---------------------------------------------------------
def mode_distribution_polar(a,l,m,u,w):
    """
    """
    radius_range = np.arange(0, 2*a ,0.01)
    Z=np.zeros(radius_range.shape,np.float64)
    Z[np.where(radius_range<=a)] = jv(l, (u[m] * radius_range[np.where(radius_range<=a)]/ a))
    Z[np.where(radius_range>a)] = jv(l,u[m]) / kv(l, w[m]) * kv(l,(w[m] * radius_range[np.where(radius_range>a)]/ a))
    return Z

def mode_distribution_cartesian(a,l,m,u,w,R,Theta):
    """
    """
    Z=np.zeros(R.shape,np.float64)
    Z[np.where(R<=a)] = jv(l, (u[m] * R[np.where(R<=a)]/ a)) * np.cos(l*Theta[np.where(R<=a)])
    Z[np.where(R>a)] = jv(l,u[m]) / kv(l, w[m]) * kv(l,(w[m] * R[np.where(R>a)]/ a)) * np.cos(l*Theta[np.where(R>a)])
    return Z