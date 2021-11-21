from decimal import Decimal
from Num_TISE import shooting

import matplotlib.pyplot as plt
import math


#Inicializacao

#parameters of potential (CHANGE)
V0 = 50 #in MeV
#A = 65 #number of nucleons
lamda = -33 #strength of spin-orbit
A_range = [*range(20,201, 20)]

#constants (DON'T CHANGE)
hc = 197.327 #in Mev.fm
m_proton = 938.272 #in Mev
R0 = 1.25 #in fm
a = 0.524 #in fm
orbital = ['s','p','d','f','g','h','i','j']

#Energy_total = []
Subshell_states = []
orbitals = []
Energy_name=[]
#Dif_allA = []
magic_indices = [1,2,4,6,7,10,12]



def takefirst(elem):
    return elem[0]

for A in A_range:
    #useful relations
    R = R0*(A**(1/3))
    E_ref = ((hc)**2)/(2*m_proton*(R**2))
    V0_adimensional = V0/E_ref
    a_adimensional = a/R

    cond_fronteira = [0,0]
    Tamanho = [10**(-10),3.5]

    Energy_Values = []
    Energy_total = []
    dif = []
    magic_measured = [0]*(len(magic_indices))
    
    degeneracy_sum = []
    subshell = []
    '''Energy_Values.clear()
    degeneracy_sum.clear()
    subshell.clear()'''

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
                    Energy_Values.append(Zeros[i])
                    #subshell has energy,degeneracy(float),name and degeneracy(string)
                    tp = str(int(2*j))+'/'+str(2)
                    subshell.append((Zeros[i],A,2*j+1,str(i+1)+(r'${%s}_{%s}$' % (orbital[l],tp)),str(2*j+1)))
                    if(A==A_range[-1]):
                        orbitals.append(str(i+1)+(r'${%s}_{%s}$' % (orbital[l],tp)))
                        Energy_name.append((Zeros[i]*E_ref,str(i+1)+(r'${%s}_{%s}$' % (orbital[l],tp))))
                    #print(subshell[-1])
                    Subshell_states.append(subshell[-1])
            #print(Zeros)
    #print("Energy_name com A=", A)
    #print(Energy_Values)

    subshell.sort(key=takefirst) #in ascending order
    #print(Energy_total)
    #print(len(subshell))


#print(orbitals)
#print(transpose_Subshell_states)
all_states = []
for orb in orbitals:
    temp = []
    for state in Subshell_states:
        if(state[3] == orb):
            temp.append((state[1],state[0]*E_ref))
    all_states.append(temp)
    #print(temp)
#print(all_states)


def takesecond(elem):
    return elem[-1][1]
#all_states.sort(key=takesecond)

#plt.figure(figsize=(4,30))
#plt.axes().set_aspect(6)

itr = 0
for k in all_states:
    itr = itr + 1
    #print(k)
    state = [*zip(*k)]
    plt.plot(state[0], state[1])
    '''if (int(itr) % 2) == 0:
        plt.text(210, state[1][-1],orbitals[all_states.index(k)],fontsize=10)
    else:
        plt.text(221, state[1][-1],orbitals[all_states.index(k)],fontsize=10)'''

Energy_name.sort(key=takefirst)
for a in range(len(Energy_name)):
    if int(a) % 2 == 0:
        plt.text(210, Energy_name[a][0], Energy_name[a][1],fontsize=10)
    else:
        plt.text(221, Energy_name[a][0], Energy_name[a][1],fontsize=10)


plt.axvline(x=65, color= 'gray', linestyle='--')
plt.axvline(x=110, color= 'gray', linestyle='--')
plt.axvline(x=155, color= 'gray', linestyle='--')
#plt.text(35, -115, "[2,8,20,28,50,82,126]", color='red')

#plt.legend(orbitals)
#plt.legend(orbitals, bbox_to_anchor=(1.02,1), loc='upper left', borderaxespad=0)
plt.xlabel("A")
plt.ylabel("E(MeV)")
plt.savefig('fig1.png', dpi=300)
plt.show()
