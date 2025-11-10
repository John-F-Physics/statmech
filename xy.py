#The classical XY vector model
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

#Hamiltonian of the XY model
def hamil(J,lattice):
    #Sum of external field term and neighbour interactions
    #+x neighbour term
    rNeighbour = -J*np.sum(np.cos(lattice-np.roll(lattice,-1,axis=1)))
    #-x neighbour term
    lNeighbour = -J*np.sum(np.cos(lattice-np.roll(lattice,+1,axis=1)))
    #+y neighbour term
    tNeighbour = -J*np.sum(np.cos(lattice-np.roll(lattice,-1,axis=0)))
    #-y neighbour term
    bNeighbour = -J*np.sum(np.cos(lattice-np.roll(lattice,+1,axis=0)))
    return lNeighbour+rNeighbour+tNeighbour+bNeighbour

class Lattice:

    def __init__(self,size,J,T):
        """
        Classical XY model with a lattice of (size×size), interaction strength J and temperature T.
        """
        #Size of lattice (N by N)
        self.size = size
        #Stength of interaction between neighbours
        self.J = J
        #Generate angles for vectors at each lattice point
        self.lattice = np.random.uniform(0,2*np.pi,size=[size,size])
        self.T = T
        self.hamilVal = hamil(self.J,self.lattice)
        self.scale = (8.617333262e-5*self.T)/self.J

    #Update the system
    def run(self):
        #Initial Hamiltonian
        initHamil = self.hamilVal
        #Create a new lattice by changing a fraction of of the vectors by a small amount
        newLattice = np.copy(self.lattice)
        #Select random indices from the lattice to modify
        indices = np.random.choice([0,1],p=[0.95,0.05],size=[self.size,self.size])
        for i in range(self.size):
            for j in range(self.size):
                if indices[i][j] == 1:
                    newLattice[i][j] += np.random.uniform(-0.1,0.1)
        #Get the Hamiltonian of the new lattice
        newLatticeHamil = hamil(self.J,newLattice)
        #Check if new lattice will be accepted using Boltzmann distribution
        #Boltzmann constant measured in electronvolts
        h = np.exp(-(newLatticeHamil-initHamil)/(8.617333262e-5*self.T))
        if h >= 1:
            self.lattice = newLattice
            self.hamilVal = newLatticeHamil
        else:
            #Acceptance criterion with probability h if h is less than 1
            if np.random.choice([True,False],p=[h,1-h]):
                self.lattice = newLattice
                self.hamilVal = newLatticeHamil

def main():
    l = Lattice(30,1,10)
    print(l.scale)
    #The characteristic parameter is k_B*T/J
    numRuns = 100000
    runPoints = np.arange(0,numRuns)
    hamilPoints = np.zeros(numRuns)
    for i in range(numRuns):
        #Progress scale
        if i%1000 == 0:
            print(i)
        l.run()
        hamilPoints[i] = l.hamilVal
    plt.plot(runPoints,hamilPoints)
    plt.show()
    #Plot the final lattice
    x = np.arange(0,30)
    y = np.arange(0,30)
    X, Y = np.meshgrid(x,y)
    X = X.flatten()
    Y = Y.flatten()
    angles = l.lattice%(2*np.pi)
    U = np.cos(angles).flatten()
    V = np.sin(angles).flatten()
    #Fetch colour array for vectors C that depends on the angle of the vectors
    plt.quiver(X,Y,U,V,angles,cmap="hsv")
    plt.show()

if __name__ == "__main__":
    main()
        
        
        
