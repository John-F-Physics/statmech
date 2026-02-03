import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

fAxis = np.linspace(-0.5,0.5,10)
fracX, fracY = np.meshgrid(fAxis,fAxis)
#Box length L
L = 5
X, Y = 5*fracX, 5*fracY
r = np.sqrt(X**2+Y**2)
#Sphere radius R
R = 2

def V(r):
    condition = (r <= R)
    print(condition)
    return condition*10000

pot = V(r)
print(np.roll(pot,2,axis=0))
