# #######################################################
# ================ Perform basic Calculation =======
# #######################################################

import numpy as np
from unc import Unc
import math
x = 3  # value 1
y = 2  # value 2
x_u = 0.5  # uncertainties 1
x_y = 0.7  # uncertainties 2
a = Unc(x, x_u)  # a = Unc(x, x_u,'R','category name','departments') to see relative uncertainty
b = Unc(y, x_y)

print('\nPerform Basic Calculation')
print("NOTE: By default the value in parenthesis represents the absolute uncertainty\n")
input("Enter any key to continue ...........")
# ########################### Addition Calculation ##################
print('First value is:', x, 'And Uncertainty is:', x_u,'Second value is:', y, 'And Uncertainty is:', x_y)
print('\nPerform Addition calculation:')
Addition = a+b
print(Addition)
#################################################

# ######## Subtraction  Calculation ##############
print('\nPerform Subtraction calculation:')
input("Enter any key to continue ...........")
Subtraction = a-a
print(Subtraction)
#################################################


# ######## Multiplication Calculation ##############
print('\nPerform Multiplication calculation:')
input("Enter any key to continue ...........")
Multiplication = a*b
print(Multiplication)
#################################################

# ######## Division Calculation ##############
print('\nPerform Division calculation:')
input("Enter any key to continue ...........")
Multiplication = a/b
print(Multiplication)
#################################################

# ######## Power Calculation ##############
print('\nPerform Power calculation:')
input("Enter any key to continue ...........")
Power = a**2
print(Power)
#################################################

# ######## Square Root Calculation ##############
print('\nPerform Square Root calculation:')
input("Enter any key to continue ...........")
sqrt = a.sqrt()
print(sqrt)
#################################################

# ######## Natural logarithm - ln(x) Calculation ##############
print('\nPerform Natural logarithm calculation:')
input("Enter any key to continue ...........")
ln = a.ln()
print(ln)
#################################################

# ######## Normal log (base-10 logarithm) -  Calculation ##############
print('\nPerform base-10 logarithm calculation:')
input("Enter any key to continue ...........")
log =a.ulog()
print(log)
#################################################

# ######## e^x -  Calculation ##############
print('\nPerform e^x calculation:')
input("Enter any key to continue ...........")
uexp =a.uexp()
print(uexp)
#################################################

# ######## sin - Calculation ##############
print('\nPerform sin calculation:')
input("Enter any key to continue ...........")
sin = np.sin(a)
print(sin)
#################################################

# ######## cos - Calculation ##############
print('\nPerform cos calculation:')
input("Enter any key to continue ...........")
cos = np.cos(a)
print(cos)
#################################################

# ######## tan - Calculation ##############
print('\nPerform tan calculation:')
input("Enter any key to continue ...........")
tan = np.tan(a)
print(tan)
#################################################

# ######## cot - Calculation ##############
print('\nPerform cot calculation:')
input("Enter any key to continue ...........")
cot = a.cot()
print(cot)
#################################################


#################################################

# ######## asin - Calculation ##############
print('\nPerform asin calculation:')
input("Enter any key to continue ...........")
c = Unc(0.6, 0.2)
asin = np.arcsin(c)
print(asin)
#################################################


# ######## acos - Calculation ##############
print('\nPerform acos calculation:')
input("Enter any key to continue ...........")
acos = np.arccos(c)
print(acos)
################################################

# ######## atan - Calculation ##############
print('Perform atan calculation:')
input("Enter any key to continue ...........")
atan = np.arctan(c)
print(atan)
################################################

# ######## sinhx - Calculation ##############
print('\nPerform sinhx calculation:')
input("Enter any key to continue ...........")
sinh = np.sinh(c)
print(sinh)
################################################

# ######## coshx - Calculation ##############
print('\nPerform coshx calculation:')
input("Enter any key to continue ...........")
cosh = np.cosh(c)
print(cosh)
################################################

# ######## tanhx - Calculation ##############
print('\nPerform tanhx calculation:')
input("Enter any key to continue ...........")
tanh = np.tanh(c)
print(tanh)
################################################

# ######## asinhx - Calculation ##############
print('\nPerform asinhx calculation:')
input("Enter any key to continue ...........")
asinh = np.arcsinh(c)
print(asinh)
################################################

# ######## acoshx - Calculation ##############
print('\nPerform acoshx calculation:')
input("Enter any key to continue ...........")
d = Unc(1.5, 0.2)
acosh = np.arccosh(d)
print(acosh)
################################################

# ######## atanhx - Calculation ##############
print('\nPerform atanhx calculation:')
input("Enter any key to continue ...........")
atanh = np.arctanh(c)
print(atanh,"end")
################################################







