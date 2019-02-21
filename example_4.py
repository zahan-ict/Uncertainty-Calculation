import numpy as np
from unc import Unc
import math

atext = """ ----------- Operations with Matrix ----------------- """
print(atext)
# ################ MATRIX ADDITION EXAMPLE 1  ################################
n = 3
size = np.ones((n, 3))

print("\nMatrix 1:\n")
A = Unc(np.matrix('11,12,13;21,22,23;31,32,44'), size,'') # Size variable assume as uncertainty values of the matrix
# A = Unc(np.matrix('11,12,13;21,22,23;31,32,44'), size, 'R')  to see relative uncertainty
print("\nMatrix 2:\n")
B = Unc(np.matrix('10,20,30;40,50,60;70,80,90'), size, '') # Matrix B

print("NOTE:By default the value inside parenthesis represents the absolute uncertainty\n")
print('\nNow perform Matrix  addition calculation.')
input("Enter any key to view result ...........\n")
C = A + B  # Sum Up
print("Result:\n", C)

###########################################################################

print('\nNow try with different uncertainty value')

# ################ MATRIX ADDITION EXAMPLE 2 ################################
input("Enter any key to continue ...........")
n = 3
size = np.ones((n, 3))
print("\nMatrix 3:\n")
A = Unc(np.matrix('11,12,13;21,22,23;31,32,44'), np.matrix('0.2,0.5,0.6;0.3,0.01,0.08;0.9,0.1,0.30')) # Size variable assume as uncertainty values of the matrix

print("\nMatrix 4:\n")
B = Unc(np.matrix('10,20,30;40,50,60;70,80,90'), size) # Matrix B
input("Enter any key to view result ...........\n")
C = A + B  # Sum Up
print("\n Result:\n", C)
print("NOTE:By default the value inside parenthesis represents the absolute uncertainty\n")
###########################################################################


print('\nNow try Matrix subtraction')

# ################ MATRIX ADDITION EXAMPLE 2 ################################
input("Enter any key to continue ...........")
n = 3
size = np.ones((n, 3))
print("\nMatrix 3:\n")
A = Unc(np.matrix('11,12,13;21,22,23;31,32,44'), np.matrix('0.2,0.5,0.6;0.3,0.01,0.08;0.9,0.1,0.30')) # Size variable assume as uncertainty values of the matrix
print("\nMatrix 4:\n")
B = Unc(np.matrix('10,20,30;40,50,60;70,80,90'), size) # Matrix B
input("Enter any key to view result ...........\n")
C = A - B  # Sum Up
print("Result:\n", C)
print("NOTE:By default the value inside parenthesis represents the absolute uncertainty\n")
###########################################################################

