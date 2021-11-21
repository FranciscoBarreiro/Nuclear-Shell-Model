import matplotlib.pyplot as plt
import math

class shooting:
    #Class variables
    epsilon = 10**-3 #for evaluating the first point after 0
    number_of_points = 10**3 #for creation of psi function
    error = 10**-3 #bissection method (may be set to 10**-7)
    Step_number = 100 #Number of points for energy graph (may be set to 500)

    #Construtor
    def __init__(self,cond_iniciais,cond_fronteira,potential,l=0): #,potential)
        self.x_min = cond_iniciais[0]
        self.x_max = cond_iniciais[1]
        self.psi0 = cond_fronteira[0]
        self.psin = cond_fronteira[1]
        self.psi1 = cond_fronteira[0] + self.epsilon
        self.V = potential
        self.l = l
        #self.potential = potential

    def k_2(self,x,E): 
        return E - self.V(x) - (self.l*(self.l+1))/(x**2)
    
    #Gera a funcao de onda
    def wave_function(self,E):
        step = (self.x_max-self.x_min)/self.number_of_points
        psi = []
        psinorm = 0
        psi.append(self.psi0)
        psinorm = psinorm + ((self.psi0)**2)*step
        psi.append(self.psi1)
        psinorm = psinorm + ((self.psi1)**2)*step
        for i in range(1,self.number_of_points):
            x_i = i*step + self.x_min
            temp = (2*(1-(5/12)*self.k_2(x_i,E)*step**2)*psi[i] - (1+(1/12)*self.k_2(x_i-step,E)*step**2)*psi[i-1])/(1+(1/12)*self.k_2(x_i+step,E)*step**2)
            psi.append(temp)
            psinorm = psinorm + (temp**2)*step

        #normalize
        for i in range(len(psi)):
            psi[i] = psi[i]/(math.sqrt(psinorm))
        
        return psi

    def Energy_Function(self,E_0,E_max):
        Energy_step = (E_max - E_0)/self.Step_number
        Energies = []
        Dif = []
        Zeros_Interval = []  
        for i in range(self.Step_number): 
            Energies.append(E_0 + Energy_step*i)
            psi = self.wave_function(E_0 + Energy_step*i)
            Dif.append(psi[self.number_of_points]-self.psin)
            if Dif[i]*Dif[i-1] < 0 and i!=0 :
                Zeros_Interval.append(Energies[i-1])
                Zeros_Interval.append(Energies[i])

        return Energies,Dif,Zeros_Interval

    def Energies_plot(self,Energies,Dif,title=""):
        n = len(Energies)
        Zeros = [0]*n
        plt.plot(Energies,Dif,'r',Energies,Zeros,'g--')
        if title != "" :
            plt.suptitle(title)
        plt.show()
        
    def bissection_method(self,a,b):
        precision = 10**(-12)
        while True :
            x = (a+b)/2
            psi = self.wave_function(x)
            if abs(self.psin-psi[self.number_of_points])<self.error or b-precision<=a<=b+precision:
                return x 
            psia = self.wave_function(a)
            #print(self.psin-psi[self.number_of_points])
            if (self.psin-psi[self.number_of_points])*(self.psin-psia[self.number_of_points]) < 0:
                b = x
            else :
                a = x
    
    def Find_Zeros(self,Zeros_Interval):
        Number_of_Zeros2x = len(Zeros_Interval)
        Zeros = []
        for i in range(0,Number_of_Zeros2x,2) :
            Zeros.append(self.bissection_method(Zeros_Interval[i],Zeros_Interval[i+1]))
            #print(i)
        return Zeros
    
    def Print_Psi(self,Zeros):
        xs = []
        step = (self.x_max-self.x_min)/self.number_of_points
        for i in range(self.number_of_points+1):
            xs.append(self.x_min + step*i)
        for i in range(len(Zeros)):
            psi = self.wave_function(Zeros[i])
            plt.plot(xs,psi)
            plt.suptitle("State"+str(i)) 
            plt.show()