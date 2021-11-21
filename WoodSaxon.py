from decimal import Decimal
from Num_TISE import shooting

import matplotlib.pyplot as plt
import math

#Inicializacao

#All values found that reprodcuce Krane: 
#[V0,A] = [95, 67], [95, 68], [96, 66], [96, 67], [97, 66], [97, 67], [98, 65], [98, 66], [99, 64], [99, 65], [99, 66], [100, 64], [100, 65], [101, 63], [101, 64], [102, 63], [102, 64], [103, 62], [103, 63], [104, 62], [105,62], [106,61]

#parameters of potential (CHANGE)
V0 = 96 #in MeV
A = 67 #number of nucleons

#constants (DON'T CHANGE)
hc = 197.327 #in Mev.fm
m_proton = 938.272 #in Mev
R0 = 1.25 #in fm
a = 0.524 #in fm

#useful relations
R = R0*(A**(1/3))
E_ref = ((hc)**2)/(2*m_proton*(R**2))
V0_adimensional = V0/E_ref
a_adimensional = a/R



def V(x):#Woods-Saxon potential, also includes centrifugal barrier
        return -((V0_adimensional)/(1+math.exp((x-1)/a_adimensional)))

cond_fronteira = [0,0]
Tamanho = [10**(-10),3.5]

Energy_Values = []
orbital = ['s','p','d','f','g','h','i','j']
degeneracy_sum = []
subshell = []

#threshold to generate the magic numbers
threshold = 7.5


#cycle through different values of angular momenta
for l in range(8): #it goes through s,p,d,f,g,h,i,j orbitals
        Function = shooting(Tamanho,cond_fronteira,V,l)
        Energies,Dif,Zeros_Interval = Function.Energy_Function(-V0_adimensional,0) 
        #Function.Energies_plot(Energies,Dif)
        Zeros = Function.Find_Zeros(Zeros_Interval)
        #Function.Print_Psi(Zeros)
        for j in range(len(Zeros)):
                Energy_Values.append(Zeros[j])
                #subshell has energy,degeneracy(float),name and degeneracy(string)
                subshell.append((Zeros[j],2*(2*l+1),str(j+1)+orbital[l],str(2*(2*l+1))))
        #print(Zeros)
#print(Energy_Values)

def takefirst(elem):
        return elem[0]
subshell.sort(key=takefirst)
#print("subshell's")
#print(subshell)


#print("diferenças de energias entre shell")
dif = [0]*(len(subshell)-1)
for k in range (len(subshell)-1):
        dif[k]=subshell[k+1][0]-subshell[k][0]
#print(dif)

for j in range(len(subshell)):
        plt.text(1.55,0.05+subshell[j][0],subshell[j][2])
        plt.text(0.4,0.05+subshell[j][0],subshell[j][3])

degeneracy_sum.append(subshell[0][1])
for n in range(len(subshell)-1):
        degeneracy_sum.append(subshell[n+1][1]+degeneracy_sum[n])
#print("Degenerescência acumulativa")
#print(degeneracy_sum)
for m in range(len(dif)):
        if(dif[m]>threshold): #energy difference to change shell is 9
                plt.text(0.9,0.05+subshell[m][0],degeneracy_sum[m])
plt.eventplot(Energy_Values,orientation='vertical',colors='k')
plt.axis("off")
plt.title("A="+str(A))
plt.show()
'''for i in range(len(Energy_Values)):
    plt.text(1.55,0.05+Energy_Values[i],text[i])
    plt.text(0.4,0.05+Energy_Values[i],degeneracy[i])
plt.eventplot(Energy_Values,orientation='vertical',colors='k')
plt.axis("off")
plt.show()
#Function.Print_Psi(Zeros)
#print(math.sqrt(Zeros[0]))'''

