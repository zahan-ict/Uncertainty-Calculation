from ucom import Ucom
import numpy as np
import cmath

# ######## complex  sin calculation ##################################################
print('\nPerform complex sin calculation where value = 0+3j and uncertainty = 0+0.1j ')
a = Ucom(0+3j, 0+0.1j)
sinx = np.sin(a)
input("\nEnter any key to continue ...........\n")
print(sinx)
print("NOTE: The value in  brackets represents the absolute uncertainty\n")
###################################################################################


# ######## complex  cos calculation ###################################################
print('\nPerform complex cos calculation where value = 0+2j and uncertainty = 0+0.3j ')

b = Ucom(0+2j, 0+0.3j)
cosx = np.cos(b)
input("\nEnter any key to see continue ...........\n")
print(cosx)
print("NOTE: The value in  brackets represents the absolute uncertainty\n")
print('\nNow add this values into another complex number 0+3j and uncertainty 0+2j \n')
input("Enter any key to see result ...........\n")
print('\nResult is:')
D=Ucom(0+3j, 0+2j)
H = sinx + cosx
HK = H + a
print(HK)
######################################################################################

