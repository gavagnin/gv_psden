import sys
import numpy as np
import pylab as pl
from binaryio import load_binary
M=2
eps2=0.05

def load(fname):
    """
this function load the binary file and save particles' position and velocity. Length is needed to know the number of particles.
    """
	    D=load_binary(fname)
	    X = D['particles']['x'][:, 0]
    	V = D['particles']['v'][:, 0]
    	length=len(X)
	return dict(X=X,V=V, length=length)
def energy(r, v):
"""
this function calculate the ernergy of each particle from its position and velocity. 
"""
    	U=-M/np.sqrt(r**2+eps2)
    	T=0.5 * v**2
    	return T+U
def r_recover(E, v):
"""
this function derive particle position from energy and velocity
"""
    	r=np.sqrt((M/(0.5*v**2-E))**2-eps2)
    	return r
def func(fname):

"""
this function 
"""
	ic=load(fname)
	x_min=np.min(abs((ic['X'])))
	x_max=np.max(abs((ic['X'])))
	v_min= np.min(abs((ic['V'])))
	v_max= np.max(abs((ic['V'])))
	E_min=energy(x_min, v_min)
	E_max=energy(x_max, v_max)
	step_E=(E_max-E_min)/(0.5*ic['length'])
	E=np.arange(E_min, E_max, step_E, dtype=np.float)
	step_v=abs(v_max-v_min)/(0.5*ic['length'])
	if step_v==0:
		v=v_min
	else:
		v=np.arange(v_min, v_max, step_v, dtype=np.float)
	x=r_recover(E,v)
	return dict(x=x, v=v, E=E, ic=ic)

zerovel_init=func('output/3_vel/zv00000.binary')
a= zerovel_init['ic']['X'][11341]
w=abs((zerovel_init['x']-a)) < 10**-5.8

print zerovel_init['x'][w], a, zerovel_init['E'][w]
zerovel_fin=func('output/3_vel/zv00096.binary')
b=abs((zerovel_fin['E']-zerovel_init['E'][w])) < 10**-4.85
print zerovel_fin['E'][b], zerovel_fin['x'][b], zerovel_fin['v'][b]
e=abs(zerovel_fin['ic']['X']-zerovel_fin['x'][b]) < 10**-5.5
print zerovel_fin['ic']['V'][e]
plusvel_init=func('output/3_vel/pv00000.binary')
q=abs((plusvel_init['E']-zerovel_init['E'][w])) < 10**-5.2

print plusvel_init['x'][q], a, plusvel_init['E'][q]
plusvel_fin=func('output/3_vel/pv00096.binary')
f=abs((plusvel_fin['E']-plusvel_init['E'][q])) < 10**-4.9
print plusvel_fin['E'][f], plusvel_fin['x'][f], plusvel_fin['v'][f]









