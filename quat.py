import numpy as np

class quat:

    def __init__(self,a,b,c,d):
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.array = np.array([a,b,c,d])
        self.abs = np.sqrt(a**2+b**2+c**2+d**2)
        self.vect = np.array([b,c,d])
        self.re = a

    def helloWorld(self):
        print("helloWorld")

    def __add__(self,other):
        return quat(*(self.array+other.array))

#Test code
h = quat(1,2,3,4)
print(h.vect)
print(h.re)
l = quat(1,5,3,4)
z = h+l
