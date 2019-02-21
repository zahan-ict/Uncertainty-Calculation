from unc import Unc

R1 = Unc(3, 0.5)
R2 = Unc(2, 0.2)
R3 = Unc(1, 0.2)
R4 = Unc(1, 0.1)

print(id(R1), id(R2), id(R3), id(R4))

K3 = R2+R1+R4


print(K3)




#print(R1+R2+R3)

