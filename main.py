import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import moduleExNum
import time
import sys


toolbar_width = 40

# setup toolbar
sys.stdout.write(__file__)
sys.stdout.write("[%s]" % (" " * toolbar_width))
sys.stdout.flush()
sys.stdout.write("\b" * (toolbar_width+1)) # return to start of line, after '['

"""""
2D Ising Model using Metropolis–Hastings algorithm

Numeric Exercise for Course Statistical Physics
written by Nir Goldfriend 
December 2017
"""""

class SpinGrid(object):
    def __init__(self, J,H, T,N=40,binary=True):
        if not isinstance(N,int): #not all(isinstance(i, float) for i in [J, H,T]) or
            raise ValueError("N must be an integer and j,H be a real number")
        self.J = J # exchange energy
        self.H = H # strength of magnetic field
        self.N = N #grid_size
        self.T = float(T) # temperature
        if binary:
            self.spin_grid = (2 * np.ceil(np.random.rand(N, N) - 0.5) - 1) / 2  # generating N*N 1/2's spin grid
        else:
            self.spin_grid = np.random.rand(N, N) - 0.5  # generating N*N  spin grid that can get each value from -1/2-1/2
        self.probMetropolis() #Metropolis–Hastings algorithm
        self.dEnergy()
        self.getParam()
        self.M = 0  # total magnetization
        self.E = 0  # total energy
        self.chi=0 # # susceptibility
        self.C_V=0 # heat_capacity
        #        self.Mlist = np.array([])  # list that holds magnetization values
#        self.Elist = np.array([])  # list that holds energy values

    def dEnergy(self,i=1,j=1,spin_state=True): # spin_state gets True if we're looking its real state of the option to flip it
        i=int(i) # casting to integers
        j=int(j) # casting to integers
        if spin_state:
            spin = self.spin_grid[i][j]
        else:
            spin = -1 * (self.spin_grid[i][j])
        energy_H = self.H * spin  # The energy from the H component
        energy_exchng = 0  # # The energy from the exchnage, we now gonna build it
        if i == 0:  # spin_grid[0,j]
            top = self.spin_grid[self.N - 1][j]
        else:
            top = self.spin_grid[i - 1][j]
        if i == self.N - 1:  # spin_grid[boundry,j]
            bottom = self.spin_grid[0][j]
        else:
            bottom = self.spin_grid[i + 1][j]
        if j == 0:  # spin_grid[i,0]
            left = self.spin_grid[i][self.N - 1]
        else: 
            left = self.spin_grid[i][j - 1]
        if j == self.N - 1:  # spin_grid[i,boundry]
            right = self.spin_grid[i][0]
        else:
            right = self.spin_grid[i][j + 1]
        energy_exchng = energy_exchng + self.J * spin * \
                                        (top + bottom + left + right)
        dE = energy_H + energy_exchng
        return dE

    def probMetropolis(self,steps=10000):
        beta = np.power(self.T, -1)
        for k in range(steps):
            i,j=np.random.randint(0,self.N,2)
            #for i in range(self.N):
          #  for j in range(self.N):
            deltaE=self.dEnergy(i,j,True)-self.dEnergy(i,j,False)
            if deltaE<0:
                self.spin_grid[i][j] = -self.spin_grid[i][j]
            elif np.random.rand()<np.exp(-1*deltaE*beta):
                self.spin_grid[i][j] = -self.spin_grid[i][j]
#            plt.imshow(self.spin_grid)
#           plt.show()
#            print(deltaE)



    def getParam(self):
        Mlist = np.array([])  # list that holds magnetization values
        Elist = np.array([])
        M=0
        for i in range(self.N):
            for j in range(self.N):
                Elist=np.append([Elist],self.dEnergy(i,j))
                Mlist = np.append([Mlist], self.spin_grid[i][j])
#                E +=self.dEnergy(self,i,j,True)
                #M +=self.spin_grid[i][j]
        E = np.mean(Elist)
        M = np.mean(Mlist)
        self.E = E
        self.M =M
#        self.Mlist=Mlist
#        self.Elist = Elist
        self.chi=np.mean((np.diff([Mlist,M*np.ones(np.size(Mlist))],axis=0)**2)) #susceptibility
        self.C_V=np.mean((np.diff([Elist,E*np.ones(np.size(Elist))],axis=0)**2)) #heat_capacity

 #       return E,M,chi,C_V

class Plots:
    def __init__(self,Tmax,J,H,Tmin=1e-10,N=40,steps=5000,binary=True):
        self.N = N  # size of lattice
        self.H = H  # strength of magnetic field
        self.J = J
        self.Tmax=Tmax
        self.Tmin=Tmin
        self.steps = steps
        self.title = 'Number of Spins: {}, J: {}, H  {},\n Number of Steps {}'.format(self.N, self.J,self.H, self.steps)
        self.T=np.linspace(Tmin,Tmax,200)
        self.E = np.array([])
        self.M = np.array([])
        self.chi = np.array([])
        self.C_V = np.array([])
        self.binary=binary
        self.x=np.array([])



    def AqcuireData(self): # building the lattices for different temperatures
        Lattice = SpinGrid(self.J, self.H, self.Tmin, self.N, self.binary)
        Lattice.probMetropolis(self.steps)
        for t in self.T:
            Lattice.T=t
            progress = int(((199 / self.Tmax)*t))
            print('Total Progress: ',progress*100/199,'%')
            sys.stdout.write(round((40//199)*progress) * "#")
            sys.stdout.flush()
            Lattice.probMetropolis(self.steps)
            Lattice.getParam()
            self.chi = np.append([self.chi], [Lattice.chi])
            self.E = np.append([self.E], [Lattice.E])
            self.C_V = np.append([self.C_V], [Lattice.C_V])
            self.M=np.append([self.M], [Lattice.M])


    def AnalitycM(self):
        beta = np.power(self.T, -1)
        g=5e-1
        x=beta*self.H
        j=0.5
        M=self.N*g*moduleExNum.brillouin(j,x)
        self.x=x
        return M


def part1():
    NumEx2017=Plots(4,0,0.6,1e-6,20,5000)
    NumEx2017.AqcuireData()
    meanM=NumEx2017.M
    meanE=NumEx2017.E
    analitycM=NumEx2017.AnalitycM()
    x=NumEx2017.x
    T=NumEx2017.T
    C_V=np.power(T,-1)*NumEx2017.C_V
    chi=np.power(T,-1)*NumEx2017.chi
    analityc_chi=np.power(T,-1)*moduleExNum.chi_anal(x)
    #analitycC_V=(1/NumEx2017.N)*NumEx2017.H*np.append(0,analitycM[1:len(T)]/np.diff(np.power(T,-1)))
    fig, ax = plt.subplots()
    ax.plot(T,C_V,linestyle='dashed', marker='.',markersize=2,label=r'$C_{V}$',linewidth=0.1)
    ax.plot(T,chi,linestyle='dashed', marker='h',markersize=2 ,label=r'$\chi$',linewidth=0.1)
    ax.plot(T,meanM,linestyle='dashed', marker='*',markersize=2, label=r'$<M>$',linewidth=0.1)
    ax.plot(T,meanE,linestyle='dashed', marker='o',markersize=2,label=r'$<E>$', linewidth=0.1)
    ax.plot(T,(1/(NumEx2017.N))*analitycM,label='Analityc <M>', linewidth=1.5)
    ax.plot(T,analityc_chi,label=r'$Analityc\chi$', linewidth=1.5)
    ax.plot(T,(1/(NumEx2017.N))*NumEx2017.H*analitycM,label=r'$Analityc<E>$', linewidth=1.5)
    ax.plot(T,np.log(2)*NumEx2017.H*analityc_chi,label=r'$AnalitycC_{V}$', linewidth=1.5)
    ax.grid(alpha=0.5)
    ax.legend(fontsize = 10)
    ax.set_title(NumEx2017.title)
    ax.set_xlabel('Temperature[J/kB]')
    #plt.ylim([0 ,max([max(C_V),max(chi),max(meanM),max(meanE)])])
    plt.plot()
    plt.show()


def part2():
    NumEx2017=Plots(2,0.9,0.00002,1e-6,40,5000,binary=True) #Tmax,J,H,Tmin=0,N=40,steps=5000,binary=True
    NumEx2017.AqcuireData()
    meanM=NumEx2017.M
    meanE=NumEx2017.E
    T=NumEx2017.T
    C_V=(np.power(T,-1)*NumEx2017.C_V)
    chi=np.power(T,-1)*NumEx2017.chi
    #analitycC_V=(1/NumEx2017.N)*NumEx2017.H*np.append(0,analitycM[1:len(T)]/np.diff(np.power(T,-1)))
  



def showing_the_lattice():
    #Part Extra for seeing how the grids look like under different parameters
    ne=SpinGrid(1,0.02,2e7,10,binary=True) #J,H, T,N=40,binary=True)
    ne.probMetropolis(500)
    ne.getParam()
    print(ne.C_V,ne.E)
    print(ne.M)
    print(ne.chi)
    fig, ax = plt.subplots()
    data = ne.spin_grid
    cax = ax.imshow(data,interpolation='spline36',cmap='YlOrBr',vmin=-0.5,vmax=0.5) # The colormap specifcly were chosen to mimick shophomre sencond year phsics labs in magnets
    ax.set_title('Lattice for T: {}; N: {};\n H: {} ; J: {}\n MetropolisIterations: {}'.format(ne.T,ne.N,ne.H,ne.J,500000000))
    # Add colorbar, make sure to specify tick locations to match desired ticklabels
    cbar = fig.colorbar(cax, ticks=[-0.5, 0.5])
    cbar.ax.set_yticklabels([r'$\frac{-1}{2}$', r'$\frac{1}{2}$'])  # vertically oriented colorbar
    plt.show()


#main progarm
def main():
    part1()
    part2()
    showing_the_lattice()
    pass
if __name__ == '__main__':
    main()