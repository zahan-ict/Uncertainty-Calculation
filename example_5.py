import numpy as np
from unc import Unc
import math

atext = """ ----------- Operations with Matrix ----------------- """
print(atext)

# ################ MATRIX MULTIPLICATION EXAMPLE 1  ################################
n = 3
size = np.ones((n, 3))
print('Lets perform Matrix Calculation.')
print("\nMatrix 1:\n")
A = Unc(np.matrix('11,12,13;21,22,23;31,32,44'), size) # Size variable assume as uncertainty values of the matrix

print("\nMatrix 2:\n")
B = Unc(np.matrix('10,20,30;40,50,60;70,80,90'), size) # Matrix B

print('\nNow perform Matrix  multiplication calculation.')
input("\nEnter any key to view result ...........\n")
C = A * B  # Sum Up
print("\nResults:\n", C)
print("NOTE:By default the value inside parenthesis represents the absolute uncertainty\n")

###########################################################################
