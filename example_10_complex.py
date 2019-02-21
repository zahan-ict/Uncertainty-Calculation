from ucom import Ucom
import numpy as np
import cmath

x = 0+3j  # value 1
x_u = 0+0.5j  # uncertainties 1
y = 0+4j
y_u = 0+.5
a = Ucom(x, x_u)
b = Ucom(x, x_u)

# ######## Complex Power Calculation ###################################################
print('\nPerform complex power calculation where value = 0+3j and uncertainty = 0+0.5j ')
input("\nEnter any key to continue ...........\n")
Power = a**2
print('Result is:', Power)
print("NOTE: The value in  brackets represents the absolute uncertainty")
input("\nEnter any key to continue ...........")

# ######## Perform complex Square Root Calculation ###########################################
print('\nPerform complex square root calculation where value = 0+4j and uncertainty = 0+0.5j ')
sqrt = b.sqrt()
print('Result is:', sqrt)
print("NOTE: The value in  brackets represents the absolute uncertainty")
input("\nEnter any key to continue ...........")
print('\nNow perform multiplication of this two value  ')
H = Power*sqrt
print('Result is:', H)
print("NOTE: The value in  brackets represents the absolute uncertainty")
###################################################################
