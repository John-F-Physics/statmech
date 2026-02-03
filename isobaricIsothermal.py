import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

#Number of particles to generate (rough estimate, this is held constant throughout the simulation but is stochastic)
numParticlesGenerate = 1
#Define mesh size (must be odd)
meshSize = 3
#Find central index of meshgrid
centralIndex = int((meshSize-1)/2)

#Meshgrid representing fractional coordinates from -1/2 to +1/2. Make sure there is an odd number of points so that the potential mesh can be centred
fracCoordAxis = np.linspace(-0.5,0.5,meshSize)
fracX, fracY, fracZ = np.meshgrid(fracCoordAxis,fracCoordAxis,fracCoordAxis)

#3D boolean array representing wherever particles are present in a space or not
P = numParticlesGenerate/(meshSize**3)
particleBoolArray = np.random.choice([0,1],p=[1-P,P],size=np.shape(fracX))
#Find true total number of particles
N = np.sum(particleBoolArray)
#List of indices of particles, representing positions of particles in the boolean array
particleIndices = []
for i in range(meshSize):
    for j in range(meshSize):
        for k in range(meshSize):
            if particleBoolArray[i][j][k] == 1:
                particleIndices.append([i,j,k])
particleIndices = np.array(particleIndices)

#Produce the potential mesh using the fractional coordinate meshes and the box size L
def potMesh(fracX,fracY,fracZ,L,epsilon,sigma):
    X, Y, Z = fracX*L, fracY*L, fracZ*L
    R = np.sqrt(X**2+Y**2+Z**2)
    #Muck about with the central point so that there is no division by zero
    R[centralIndex][centralIndex][centralIndex] = 1
    #Potential has been halved due to double counting pair potentials
    V = 2*epsilon*((sigma**12)/(R**12)-(sigma**6)/(R**6))
    #Change central point in potential to zero so that the particle doesn't interact with itself, this ain't no Higgs field
    V[centralIndex][centralIndex][centralIndex] = 0
    return V

#Function that rolls the potential mesh so that the central point lines up with the required coordinates
#The rolling is required so that periodic boundary conditions are satisfied
def lineup(V,indexX,indexY,indexZ):
    #Find out how far to roll
    rollX = indexX-centralIndex
    rollY = indexY-centralIndex
    rollZ = indexZ-centralIndex
    #X-axis is axis 1, Y-axis is axis 0, Z-axis is axis 2
    rolledV = np.roll(V,rollX,axis=1)
    rolledV = np.roll(rolledV,rollY,axis=0)
    rolledV = np.roll(rolledV,rollZ,axis=2)
    return rolledV

epsilon = 1
sigma = 1
L = 5
V = potMesh(fracX,fracY,fracZ,L,epsilon,sigma)
print(V)
    
