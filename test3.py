from unc import Unc
Ug=0
R2=1e3
R3= 1e3
R4= 1.001e3
U0= 5

u_R2 = 5
u_R3 = 5
u_R4 = 1

u_U0 = 1.00
u_Ug = 0.01

R2 = Unc(R2, u_R2)
R3 = Unc(R3, u_R3)
R4 = Unc(R4, u_R4)
U0 = Unc(U0, u_U0)
Ug = Unc(Ug, u_Ug)

I2=U0/(R3+R4)
print(I2)
I1= (I2*R4+Ug)/R2
R1 = U0/I1-R2
print(R1)
