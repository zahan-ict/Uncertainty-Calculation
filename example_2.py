import numpy as np
from unc import Unc
import math
a = Unc(3, 0.5)
b = Unc(2, 0.7)

print('\nUncertainty value with constant calculation:')
print("NOTE: By default the value in parenthesis represents the absolute uncertainty\n")
input("Enter any key to continue ...........")
# ######## Constant calculation with addition ##################
# ########################### Addition Calculation ##################
print('Value and uncertainty is:', a)
print('Now if we Add(+) pi=3.14 Constant right side of the calculation, example: Unc(3, 0.5) + Unc(math.pi)')
c = a + Unc(math.pi)
print(c)


input("Enter any key to continue ...........")
print('\nNow if we Add(+) pi=3.14 Constant left side of the calculation, example: Unc(math.pi) + Unc(3, 0.5)')
d = Unc(math.pi) + a
print(d)

#################################################################
input("Enter any key to continue ...........\n")
print('Now if we Subtract (-) pi=3.14 Constant right side of the calculation, example: Unc(3, 0.5) - Unc(math.pi)')
# ######## Constant calculation with subtraction ##################
c = a - Unc(math.pi)
print(c)

print('Now if we Subtract (-) pi=3.14 Constant left side of the calculation, example: Unc(math.pi) - Unc(3, 0.5)')
d = Unc(math.pi) - a
print(d)
input("Enter any key to continue ...........\n")
#################################################################


# ######## Constant calculation with multiplication ##################
print('\nNow if we Multiply (*) pi=3.14 Constant right side of the calculation, example: Unc(3, 0.5) * Unc(math.pi)')
c = a * Unc(math.pi)
print(c)

print('\nNow if we Multiply (*) pi=3.14 Constant left side of the calculation, example: Unc(math.pi) * Unc(3, 0.5)')
d = Unc(math.pi) * a
print(d)
input("Enter any key to continue ...........")
#####################################################################
# ######## Constant calculation with division ##################
a = Unc(3, 0.5)
print('pi is ',math.pi)
print('\nNow if we divided (/) pi=3.14 Constant as  example: Unc(3, 0.5) / Unc(math.pi)sss')
c = a / Unc(math.pi)
print(c)

print('\nNow if we divided (/) pi=3.14 Constant as  example: Unc(math.pi) / Unc(3, 0.5)')
d = Unc(math.pi) / a
print(d)
#################################################################

input("Enter any key to continue ...........\n")
# ######## Natural logarithm - ln(x) calculation with pi constant ##############
print('\nNow if add (+) pi=3.14 Constant  with log as  example: a.ulog() + Unc(math.pi)')
k = a.ulog() + Unc(math.pi)
print(k)

# #################################################################

