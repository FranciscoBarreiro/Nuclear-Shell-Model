from decimal import Decimal
from Num_TISE import shooting

import matplotlib.pyplot as plt
import math


#Inicializacao

#parameters of potential (CHANGE)
V0 = 50 #in MeV
A = 175 #number of nucleons

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

cond_fronteira = [0,0]
Tamanho = [10**(-10),3.5]

def takefirst(elem):
    return elem[0]


Energy_Values = []
orbital = ['s','p','d','f','g','h','i','j']
degeneracy_sum = []
subshell = []
magic_numbers = [2,8,20,28,50,82,126]


lamda_0 = -25 #initial lambda
lamda_f = -50 #final lambda
lamda_step = -1 #increment of lambda

#cycle for lambda 
for lamda in range(lamda_0,lamda_f,lamda_step):
    print("Value of lamda is",lamda)
    Energy_total = []
    Energy_Values = []
    #cycle through different values of angular momenta
    for l in range(8): #it goes through s,p,d,f,g,h,i,j orbitals
        if l!=0:
            j_possibilities = [l-(1/2),l+(1/2)]
        else:
            j_possibilities = [1/2]

        for j in j_possibilities:
            def V(x):#Woods-Saxon potential, also includes centrifugal barrier
                return -((V0_adimensional)/(1+math.exp((x-1)/a_adimensional))) + lamda*(1/2)*(1/x)*(E_ref/(m_proton))*(j*(j+1)-l*(l+1)-(3/4))*(V0_adimensional/a_adimensional)*((math.exp((x-1)/a_adimensional))/((1+math.exp((x-1)/a_adimensional))**2))
            Function = shooting(Tamanho,cond_fronteira,V,l)
            Energies,Dif,Zeros_Interval = Function.Energy_Function(-V0_adimensional,0)
            #Function.Energies_plot(Energies,Dif)
            Zeros = Function.Find_Zeros(Zeros_Interval)
            #Function.Print_Psi(Zeros)
            for i in range(len(Zeros)):
                Energy_total.append([Zeros[i],str(i+1)+orbital[l]+str(int(2*j))+'/'+str(2),int(2*j+1)])
                #subshell has energy,degeneracy(float),name and degeneracy(string)
                #subshell.append((Zeros[i],2*j+1,str(i+1)+orbital[l]+str(int(2*j))+'/'+str(2),str(2*j+1)))
            
    Energy_total.sort(key=takefirst)
    #print(Energy_total)

    #calculate degeneracy
    Deg_acumulated = []
    temp = 0
    for i in range(len(Energy_total)):
        temp = temp + Energy_total[i][2]
        Deg_acumulated.append(temp)
    #print(Deg_acumulated)

    dif = []
    for i in range(len(Energy_total)-1):
        if Deg_acumulated[i]<184:
            dif.append([Energy_total[i+1][0] - Energy_total[i][0],Deg_acumulated[i]])
    dif.sort(key=takefirst,reverse=True)
    #print(dif)

    #print("Interval :", dif[len(magic_numbers)][0], "-",dif[len(magic_numbers)-1][0])

    magic_prediction = []
    for k in range(len(magic_numbers)):
        magic_prediction.append(dif[k][1])

    magic_prediction.sort()
    #print(magic_prediction)

    if magic_prediction == magic_numbers:
        print("\nGenerates the correct magic numbers!!\n")
        print("Interval of Difference to generate magic numbers:", dif[len(magic_numbers)][0], "-",dif[len(magic_numbers)-1][0]+'\n')
    
    



'''subshell.sort(key=takefirst)
print("subshell's")
print(subshell)


print("diferenças de energias entre shell")
dif = [0]*(len(subshell)-1)
for k in range (len(subshell)-1):
        dif[k]=subshell[k+1][0]-subshell[k][0]
print(dif)

for j in range(len(subshell)):
        plt.text(1.55,0.05+subshell[j][0],subshell[j][2])
        plt.text(0.4,0.05+subshell[j][0],subshell[j][3])

degeneracy_sum.append(subshell[0][1])
for n in range(len(subshell)-1):
        degeneracy_sum.append(subshell[n+1][1]+degeneracy_sum[n])
print("Degenerescência acumulativa")
print(degeneracy_sum)
for m in range(len(dif)):
        if(dif[m]>9): #energy difference to change shell is 9
                plt.text(0.9,0.05+subshell[m][0],degeneracy_sum[m])
plt.eventplot(Energy_Values,orientation='vertical',colors='k')
plt.axis("off")
plt.title("A="+str(A))
plt.show()'''