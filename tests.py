import numpy as np
import matplotlib.pyplot as plt

#Lennard-Jones potential
def lj(sigma,epsilon,r):
    largeR = r >= 1
    return largeR*4*epsilon*(np.power(sigma/r,12)-np.power(sigma/r,6))

sigma = 1
epsilon = 1

x = np.linspace(-3,+3,1000)
y = np.linspace(-3,+3,1000)

X, Y = np.meshgrid(x,y)

R = np.sqrt(X**2+Y**2)
V = lj(sigma,epsilon,R)
plt.pcolormesh(X,Y,V,cmap="plasma")
plt.show()


