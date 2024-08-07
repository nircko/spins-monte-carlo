from scipy.integrate import  quad
import numpy as np
from scipy import constants
import matplotlib.pyplot as plt

#Ex 8 numeric integration
#working in SI units
sample_rate=1e3
eta=1
A_D=100**2*(1e-6)**2 # in meter squared
A_L=1*(1e-2)**2 # in meter squared
f=2*1e-2 # in meters
R=100e6 # resistance
wave_len=np.array([[3.5,5],[8,14]])*1e-6
T_noise=1200 # Temperature noise
T=[300.1,300]
#T_target=301
#T_background=300
I_shot=0
c2=lambda T,wav: constants.h*constants.speed_of_light/(constants.k*T*wav)
I_jhonson=np.sqrt(4*constants.k*T_noise/R)
print("I_jhonson=",I_jhonson)
for i in range(len(wave_len)):
    for t in T:
        print('window=',wave_len[i,0],":\t:",wave_len[i,1])
        print('T=',t)
        a=c2(t,wave_len[i,1])
        print("a=",a)
        radiant_emittance = lambda x: x**-4*(np.exp(a/x)-1)**-1
        W_at= quad(radiant_emittance,wave_len[i,0]/wave_len[i,1],1)
        print("W_at_T=",t,"is=\t", W_at,'\n')
        I_at= W_at[0]*2*constants.speed_of_light*constants.elementary_charge*A_L*A_D*eta/(np.pi*f**2*wave_len[i,1]**3)
        print("For t=",t,"K\t", "\tI_at=",I_at)
        I_shot = np.sqrt(2 * constants.elementary_charge * I_at * sample_rate)
        print("Shot noise=",I_shot,"\n\n\n\n\n#########################")
    temperature=T[1]
    NETD=f**2*(I_shot+I_jhonson)*np.pi*constants.k*temperature**2/(2*A_D*A_L*eta*constants.elementary_charge*constants.h*constants.speed_of_light**2)
    c1=constants.h*constants.speed_of_light/(constants.k*temperature)
    dQdT= lambda x:np.exp(c1/x)/((np.exp(c1/x)-1)*x**5)
    integrate_dQdT=quad(dQdT,wave_len[i,0],wave_len[i,1])
    print("NETD=",NETD/integrate_dQdT)


