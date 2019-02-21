
import numpy as np
from unc import Unc
import math

atext = """ ----------- Operations with Arrays ----------------- """
print(atext)

# ################ Statical Uncertainty   Example   ################################
print('\nLet’s assume that length of an object has measured five times')
a = np.array([72, 77, 82, 85, 88, 90, 86, 70, 91, 95])
print('\n', a)
result = Unc(a, '', '', 'Math', 'Gp1') # Size variable assume as uncertainty values of array

input("\nEnter any key to see results  ...........")

print(result.statUnc(),'\n')
# ############################################################################################


# ################ Sum up Uncertainty   Example   ################################
input("\nEnter any key to continue  ...........")
print('Lets sum up multiple uncertainty')

a = np.array([1, 2, 3, 4, 5, 6, 7, 7, 9, 10])
b = np.array([0.005, 0.02, 0.003, 0.04, 0.005, 0.06, 0.07, 0.07, 0.09, 0.02])
print("\nValues and Uncertainties  are:\n", a, '\n', b)
A = Unc(a, b) # Size variable assume as uncertainty values of array
A = Unc(a, b, '') # Size variable assume as uncertainty values of array

input("\nEnter any key to see results  ...........")
print(' Results is:', A.sum())

# ############################################################################################

input("\nEnter any key to continue ...........")

# ################ Adding Multiple  Categories of Uncertainty   Example   ################################
print("\n Let's perform group of array addition calculation:")
a = np.array([1, 2, 3, 4, 5, 6, 7, 7, 9, 10])
b = np.array([0.5, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.7, 0.9, 0.2])
print("\n(a) Values and Uncertainties  are:\n", a, '\n', b)
A = Unc(a, b) # Size variable assume as uncertainty values of array

x = np.array([11, 12, 13, 14, 15, 16, 17, 18, 19, 20])
y = np.array([0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.3, 0.1, 0.2, 0.3])
print("\n(b) Values and Uncertainties  are:\n", x, '\n', y)

B = Unc(x, y) # Array B

input("\nEnter any key to see results  ...........")
C = A+B   # Sum Up
print("\nResults:", C)
print("\nNOTE: By default the value in parenthesis represents the absolute uncertainty\n")

#####################################################################################

input("\nEnter any key to continue ...........")

# ################ Array Multiplication(x) Example  ################################
print("\nNow let's try basic multiplication C = Ax")

x = np.array([11, 12, 13, 14, 15, 16, 17, 18, 19, 20])
y = np.array([0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.3, 0.1, 0.2, 0.3])

print("\n(c) Values and Uncertainties  are:\n", x, '\n', y)

input("\nEnter any key to continue ...........")


A = Unc(x, y) # Matrix B
x = Unc(2)

C = A * x  # Sum Up
print("Results:", C)

print("\nNOTE: By default the value in parenthesis represents the absolute uncertainty\n")

#####################################################################################

input("\nEnter any key to continue ...........")

print("\nNow let's try basic division C = A/x")
# ################ Array Multiplication(x) Example  ################################
x = np.array([11, 12, 13, 14, 15, 16, 17, 18, 19, 20])
y = np.array([0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.3, 0.1, 0.2, 0.3])

print("\n(d) Values and Uncertainties  are:\n", x, '\n', y)
A = Unc(a, b)

input("\nEnter any key to continue ...........")
x = Unc(2)
C = A / x  # Sum Up
print("\nResults:", C)
print("\nNOTE: By default the value in parenthesis represents the absolute uncertainty\n")

#####################################################################################
