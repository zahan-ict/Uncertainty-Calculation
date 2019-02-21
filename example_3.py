import numpy as np
from unc import Unc
import math
a = Unc(30, 0.2)

# ######## trigonometric calculation on degree ##################
print('By default trigonometric calculations  are performed in "Radian". If we want to get result in "Degree" than we need to call degree function\n')

print('Let Value and uncertainty is :', a)
print("NOTE:The value in parenthesis represents the absolute uncertainty\n")

input("Enter any key to continue ...........\n")
print('Lets perform sin calculation')
sin = Unc.degree(a, 'sin')
print(sin)

input("\nEnter any key to continue ...........\n")
print('Lets perform cos calculation')
cos = Unc.degree(a, 'cos')
print(cos)
input("\nEnter any key to continue ...........\n")
print('Lets perform tan calculation')
tan = Unc.degree(a, 'tan')
print(tan)
input("\nEnter any key to continue ...........\n")
print('Lets perform cot calculation')
cot = Unc.degree(a, 'cot')
print(cot)


#################################################################



