from ucom import Ucom
import numpy as np
import cmath
# atext = """ ----------- Complex Operations  ----------------- """
# print(atext)

# ######## Complex Addition Calculation ##################
print("\n==============Perform complex Addition calculation==================")
print("\nFirst value = 1+2j and  uncertainty = 0+.2j")
print("Second value = 5+3j  and  uncertainty = 0+.1j ")

input("\nEnter any key to get result ...........")
H = Ucom(2j, 0.2j)
K = Ucom(3j, 0.1j)
C = H + K
print(C)
print("\nNOTE: By default  value inside the bracket represents absolute uncertainty\n")
###########################################################################

input("\nEnter any key to continue ...........")

# ######## Complex Subtraction  Calculation #################################
print("\n==============Lets perform complex  subtraction calculation==================")
print("\nFirst value = 1+2j and  uncertainty = 0+.2j")
print("Second value = 5+3j  and  uncertainty = 0+.1j ")

input("\nEnter any key to get result ...........\n")

H = Ucom(1+2j, 0+.2j)
K = Ucom(5+3j, 0+.1j)
C = H - K
print(C)

print("\nNOTE: By default  value inside the bracket represents absolute uncertainty\n")
###########################################################################


input("\nEnter any key to continue ...........\n")

# ######## Complex Multiplication  Calculation #################################
print("\n==============Lets perform  complex  multiplication calculation==================")
print("\nFirst value = 1+2j and  uncertainty = 0+.2j")
print("Second value = 5+3j  and  uncertainty = 0+.1j ")

input("\nEnter any key to get result ...........\n")

H = Ucom(1+2j, 0+.2j)
K = Ucom(5+3j, 0+.1j)
C = H * K
print(C)
print("\nNOTE: By default the value in parenthesis represents the absolute uncertainty\n")
###########################################################################

input("\nEnter any key to continue ...........\n")

# ######## Complex division  Calculation #################################
print("\n==============Lets perform  complex  division calculation==================")
print("\nFirst value = 1+2j and  uncertainty = 0+.2j")
print("Second value = 5+3j  and  uncertainty = 0+.1j ")

input("\nEnter any key to get result ...........\n")

H = Ucom(1+2j, 0+.2j)
K = Ucom(5+3j, 0+.1j)
C = H / K
print(C)
print("\nNOTE: The value in  brackets represents the absolute uncertainty\n")
###########################################################################
