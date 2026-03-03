import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

#Number of particles to generate (rough estimate, this is held constant throughout the simulation but is stochastic)
numParticlesGenerate = 10
#Define mesh size (must be odd)
meshSize = 51
#Define pressure of system (Pa)
pressure = 1.01e+5
#Define temperature of system (K)
temperature = 298
#Define physical constants
#Lennard-Jones constants
epsilon = (1.77/(6.0221408e+23))*1000
sigma = 4.10e-10
L = 50e-10
#Boltzmann constant
k = 1.38e-23

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

#Find the full potential energy of the system
def getFullEnergy(fracX,fracY,fracZ,L,epsilon,sigma,particleBoolArray,particleIndices):
    #Create the meshgrid of potential values
    V = potMesh(fracX,fracY,fracZ,L,epsilon,sigma)
    #Initialise total energy
    totEnergy = 0
    #Cycle through each particle
    for i in range(0,len(particleIndices)):
        #Get index of particle i
        index = particleIndices[i]
        shiftedV = lineup(V,*index)
        totEnergy += np.sum(particleBoolArray*shiftedV)
    return totEnergy

#Main Monte Carlo simulation


#Get the initial system energy
currentEnergy = getFullEnergy(fracX,fracY,fracZ,L,epsilon,sigma,particleBoolArray,particleIndices)
#System initial volume
vol = L**3


#Make small updates to the system
for j in range(100):
    #Make a small fractional change to L
    newL = L+np.random.uniform(-0.1,0.1)*L
    #Choose a particle from the particleIndices array
    chosenParticle = np.random.randint(0,len(particleIndices))
    #Modify the index, but make sure the indices are within bounds
    newIndex = particleIndices[chosenParticle]+np.random.randint(-5,6,size=3)
    newIndex += (newIndex > meshSize-1)*(-meshSize)
    newIndex += (newIndex < 0)*(meshSize-1)
    #Create new particle arrays
    newBoolArray = np.copy(particleBoolArray)
    newBoolArray[*newIndex] = 1
    newBoolArray[*particleIndices[chosenParticle]] = 0
    newIndices = np.copy(particleIndices)
    newIndices[chosenParticle] = newIndex
    #Find energy of new system configuration
    newEnergy = getFullEnergy(fracX,fracY,fracZ,newL,epsilon,sigma,newBoolArray,newIndices)
    #Find the Metropolis acceptance criterion for the new state in the Markov chain
    #Define new volume
    newVol = newL**3
    acceptanceCriterion = ((newVol/vol)**N)*np.exp(-(newEnergy-currentEnergy+pressure*(newVol-vol))/(k*temperature))
    print(acceptanceCriterion)
    #Accept the new state if the acceptance criterion is greater than 1
    if acceptanceCriterion >= 1:
        L = newL
        currentEnergy = newEnergy
        particleBoolArray = newBoolArray
        vol = newVol
        particleIndices = newIndices
    #Accept the new state with a finite probability
    else:
        zetta = np.random.uniform(0,1)
        if acceptanceCriterion > zetta:
            #Accept the new state
            L = newL
            currentEnergy = newEnergy
            particleBoolArray = newBoolArray
            vol = newVol
            particleIndices = newIndices
