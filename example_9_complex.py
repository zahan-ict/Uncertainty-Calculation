from ucom import Ucom
import numpy as np
import cmath

atext = """ -------------------------Complex Array Operations----------------------------------- """
print(atext)

# ################ Complex Statical Uncertainty   Example   ################################
print("\n=======Let’s assume that complex length of an object has measured six times=========")

a = np.array([1+3j, 2+5j, 3+5j, 2+1j, 0+1j, 10+3j])
print("\nComplex Measurements are:\n", a)
H = Ucom(a)
input("\nEnter any key to see results  ...........")
print('Results is:', H.comStatisticUnc())
############################################################################################################


# ######## Perform Complex Array Addition Calculation ###########################################################
print("\n==============Let's sum up multiple complex uncertainties==================")

a = np.array([1+3j, 2+5j, 3+5j, 2+1j, 0+1j, 10+3j])
a_unc = np.array([0+0.1j, 0+0.002j, 0+0.001j, 0+0.301j, 0+0.003j, 0+0.004j])
print("\nComplex Values and Uncertainties  are:\n", a, '\n', a_unc)
H = Ucom(a, a_unc)

input("\nEnter any key to see results  ...........")
print('Results is:', H.Csum())
input("\nEnter any key to continue ...........")
############################################################################################################


# ################ Adding Multiple  Categories of Uncertainty   Example   ################################
print("\n==============Adding Multiple  Categories of complex Uncertainty ==================")

print("\n Let's perform group of complex array addition calculation:")
a = np.array([1+3j, 2+5j, 3+5j, 2+1j, 0+1j, 10+3j])
b = np.array([0+0.1j, 0+0.002j, 0+0.001j, 0+0.301j, 0+0.003j, 0+0.004j])
print("\n(a) Complex Values and Uncertainties  are:\n", a, '\n', b)
A = Ucom(a, b) # Size variable assume as uncertainty values of array

x = np.array([0+1j, 1+2j, 0+3j, 5+1j, 0+1j, 0+8j])
y = np.array([0+0.3j, 0+0.004j, 0+0.007j, 0+0.80j, 0+0.20j, 0+0.300j])
print("\n(b) Complex Values and Uncertainties  are:\n", x, '\n', y)

B = Ucom(x, y) # Array B

input("\nEnter any key to see results  ...........")
C = A+B   # Sum Up
print("\nResults:", C)
print("\nNOTE: By default the value in parenthesis represents the absolute uncertainty\n")


