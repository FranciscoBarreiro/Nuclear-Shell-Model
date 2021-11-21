from Num_TISE import shooting

import matplotlib.pyplot as plt
import math


#Inicializacao


#constants (DON'T CHANGE)
hc = 197.327 #in Mev.fm
m_proton = 938.272 #in Mev
R0 = 1.25 #in fm
a = 0.524 #in fm


cond_fronteira = [0,0]
Tamanho = [10**(-10),3.5]


orbital = ['s','p','d','f','g','h','i','j']
#degeneracy_sum = []
#subshell = []
A_vals = []
Energy_total = []

#orbital order
orbital_order = ['1s','1p','1d','2s','1f','2p','1g','2d','3s','1h','2f','3p','1i','2g','3d','4s']
p = [4,3,3,2,2,1,1]
magic_indices = [1,2,4,6,7,10,12]

good_values = []
#print("We have V0 = ",V0)

A0 = 60 #initial A
Af = 70 #final A (not inclusive)
inc_A = 1 #increment between A's


V0_ref = 95 #initial V0
V0_final = 110 #final V0 (not inclusive)
inc_V0 = 1 #increment between V0's

def takefirst(elem):
    return elem[0]

for V0 in range(V0_ref,V0_final,inc_V0):
    print("V0 value is", V0)
    for A in range(A0,Af,inc_A):
        print("A value is", A)
        break_happened = 0
        Energy_total = []
        ordered_names = []
        dif = []
        magic_measured = [0]*(len(magic_indices))

        #useful relations
        R = R0*(A**(1/3))
        E_ref = ((hc)**2)/(2*m_proton*(R**2))
        V0_adimensional = V0/E_ref
        a_adimensional = a/R
        #print(A)

        #Woods-Saxon potential, also includes centrifugal barrier
        def V(x,V0_adimensional=V0_adimensional,a_adimensional=a_adimensional):
            return -((V0_adimensional)/(1+math.exp((x-1)/a_adimensional)))
        
        
        for l in range(7): #it goes through s,p,d,f,g,h,i orbitals
            Energy_Values = []
            Function = shooting(Tamanho,cond_fronteira,V,l)
            Energies,Dif,Zeros_Interval = Function.Energy_Function(V(0),0)
            #Function.Energies_plot(Energies,Dif)
            Zeros = Function.Find_Zeros(Zeros_Interval)
            #Function.Print_Psi(Zeros)
            #print(Zeros)
            if len(Zeros)>=p[l]:
                for j in range(p[l]):
                    Energy_Values.append([Zeros[j],str(j+1)+orbital[l]])
            else:
                break_happened = 1
                print("CAN'T GENERATE ALL ORBITALS -- POTENTIAL NOT DEEP ENOUGH!")
                break

            #print(Energy_Values)
            for j in range(len(Energy_Values)):
                Energy_total.append(Energy_Values[j])
            #print(Energy_total)
        
        if break_happened == 0:
            Energy_total.sort(key=takefirst) #in ascending order
            #print(Energy_total)
            #print(len(Energy_total))
            #Only get names
            for k in range(len(Energy_total)):
                ordered_names.append(Energy_total[k][1])

            #print(ordered_names)
            #Get differences
            for k in range (len(Energy_total)-1):
                dif.append([Energy_total[k+1][0]-Energy_total[k][0],k+1])
            
            dif.sort(reverse=True, key=takefirst) #in descending order
            for k in range(len(magic_indices)):
                magic_measured[k] = dif[k][1]

            magic_measured.sort()
            #print(magic_measured)
            if magic_measured==magic_indices and ordered_names==orbital_order:
                good_values.append([V0,A])
                print("\nFOUND!\n")
                #print(magic_measured)
                #print(ordered_names)

print("The values of [V0,A] that generate Krane's spectrum are ")
print(good_values)




        
                

                        


'''
A0 = 40
for i in range(9):
        break_happened = 0
        A = A0 + i*5
        #useful relations
        R = R0*(A**(1/3))
        E_ref = ((hc)**2)/(2*m_proton*(R**2))
        V0_adimensional = V0/E_ref
        a_adimensional = a/R
        print(A)
        def V(x,V0_adimensional=V0_adimensional,a_adimensional=a_adimensional):#Woods-Saxon potential, also includes centrifugal barrier
                return -((V0_adimensional)/(1+math.exp((x-1)/a_adimensional)))
        Energy_Values = []
        for l in range(7): #it goes through s,p,d,f,g,h,i,j orbitals
                Function = shooting(Tamanho,cond_fronteira,V,l)
                Energies,Dif,Zeros_Interval = Function.Energy_Function(V(0),0)
                #Function.Energies_plot(Energies,Dif)
                Zeros = Function.Find_Zeros(Zeros_Interval)
                #Function.Print_Psi(Zeros)
                #for j in range(len(Zeros)):
                #        Energy_Values.append(Zeros[j])
                #        subshell.append((Zeros[j],2*(2*l+1),str(j+1)+orbital[l],str(2*(2*l+1))))
                #print(Zeros)
                if len(Zeros)>=p[l]:
                        Energy_Values.append(Zeros[p[l]-1])
                        #Energy_total[i].append(Zeros[p[l]])
                else:
                        break_happened = 1
                        print("DIDN'T FIND")
                        break
        #print(Energy_Values)
        if break_happened==0:
                print("FOUND!")
                A_vals.append(A)
                #print(Energy_Values)
                Energy_total.append(Energy_Values)
                #print(Energy_total)
                '''

'''
transpose_Energy_total = [*zip(*Energy_total)]
#print(transpose_Energy_total)

clr = ["b","r","g","y","c","m","orange"]
for k in range(len(clr)):
        plt.plot(A_vals,transpose_Energy_total[k], color=clr[k])
plt.show()
'''
'''
def takefirst(elem):
        return elem[0]
subshell.sort(key=takefirst)
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
plt.show()
'''



'''for i in range(len(Energy_Values)):
    plt.text(1.55,0.05+Energy_Values[i],text[i])
    plt.text(0.4,0.05+Energy_Values[i],degeneracy[i])
plt.eventplot(Energy_Values,orientation='vertical',colors='k')
plt.axis("off")
plt.show()
#Function.Print_Psi(Zeros)
#print(math.sqrt(Zeros[0]))'''