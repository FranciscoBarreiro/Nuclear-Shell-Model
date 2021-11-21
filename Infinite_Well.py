from Num_TISE import shooting

import matplotlib.pyplot as plt
import math

#Inicializacao
def V(x): #infinite well
        return 0
cond_fronteira = [0,0]
Tamanho = [0.00001,1]

Energy_Values = []
orbital = ['s','p','d','f','g','h','i','j']
theoretical_zeros= [[3.142,6.283,9.425],[4.493,7.725,10.904],[5.763,9.095],[6.988,10.417],[8.183,11.705],[9.35581],[10.5128],[11.657]]

#threshold to generate the magic numbers
threshold = 9

degeneracy_sum = []
subshell = []
#cycle through different values of angular momenta
for l in range(8): #it goes through s,p,d,f,g,h,i,j orbitals
        Function = shooting(Tamanho,cond_fronteira,V,l)
        #Plot de Energias entre 0 e 80
        Energies,Dif,Zeros_Interval = Function.Energy_Function(0,140)
        #Function.Energies_plot(Energies,Dif,"Infinite potential well for l="+str(l))
        Zeros = Function.Find_Zeros(Zeros_Interval)
        #Function.Print_Psi(Zeros)
        for j in range(len(Zeros)):
                Energy_Values.append(Zeros[j])
                subshell.append((Zeros[j],2*(2*l+1),str(j+1)+orbital[l],str(2*(2*l+1))))
                error = (abs(Zeros[j]-(theoretical_zeros[l][j]**2))/(theoretical_zeros[l][j]**2))*100
                print("The error for l=",l," and n=", j+1," is")
                print(error,"%")

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

for i in range(len(subshell)):
    plt.text(1.55,0.05+subshell[i][0],subshell[i][2])
    plt.text(0.4,0.05+subshell[i][0],subshell[i][3])

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
plt.show()
#Function.Print_Psi(Zeros)
#print(math.sqrt(Zeros[0]))
