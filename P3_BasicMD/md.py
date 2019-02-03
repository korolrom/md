import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
import num_int as num
import SHO as sho
#input parameters:
m=1.
omega=1.
dt=1e-1
t=10.

#future sampler
def sample_q():
    return np.random.rand()
def sample_p():
    return np.random.rand()+1
    #p is shifted so there are no particles in the classically forbidden region
p0=sample_p()
q0=sample_q()
e=sho.get_ham(p0,q0,omega,m)

#Analytical solution for harmonic oscillator
amp=np.sqrt(2*e/m)/omega
phi=np.arcsin(q0/amp)
da=-omega**2
AnalPQ=np.array([sho.anal_p(amp,omega,m,phi,t,dt), sho.anal_q(amp,omega,m,phi,t,dt)])
AnalE=np.array([sho.get_ke(AnalPQ[0,:],m),sho.get_pe(AnalPQ[1,:],omega,m)])

#Numerical integration of EoM to get trajectories
PQ=np.array([AnalPQ, num.Euler(p0,q0,t,dt,m,da),num.Verlet(p0,q0,t,dt,m,da),
            num.vVerlet(p0,q0,t,dt,m,da),num.Beeman(p0,q0,t,dt,m,da)])
#print(np.shape(PQ))

#get the energies on the trajectory
E=np.moveaxis(np.array([sho.get_ke(PQ[:,0,:],m),sho.get_pe(PQ[:,1,:],omega,m)]),1,0)
#print(np.shape(E))
#print(np.shape(PQ))

#Compute global energy conservation and the standard deviation in E
Ebar=np.cumsum(np.sum(E,axis=1),axis=-1)/(np.arange(int(t/dt))+1)
print(np.shape(Ebar))
print(np.shape((np.sum(E,axis=1)-Ebar)**2))
Esigma=np.cumsum((np.sum(E,axis=1)-Ebar)**2,axis=-1)/(np.arange(int(t/dt))+1)

# Plots of p and q with time, energies with time
n=np.size(E,axis=0)
Style=[['r-','r--','r:','r-.'],['g-','g--','g:','g-.'],['b-','b--','b:','b-.'],['k-','k--','k:','k-.'],['m-','m--','m:','m-.'],['c-','c--','c:','c-.']]
fig = plt.figure()
subj=1
for item in [PQ,E]:
    axes_1 = fig.add_subplot(1, 2, subj)
    for i in range(n):
        plt.plot(np.arange(0,t,dt),item[i,0,:],Style[i][0])
        plt.plot(np.arange(0,t,dt),item[i,1,:],Style[i][1])
        if subj==2:
            plt.plot(np.arange(0,t,dt),np.sum(item[i,:,:],axis=0),Style[i][2])
            plt.plot(np.arange(0,t,dt),Ebar[i,:],Style[i][3])
            plt.plot(np.arange(0,t,dt),np.sqrt(Esigma[i,:]),Style[i][3])
    subj=subj+1
plt.show()
#print(p0,q0)


# x=sp.symbols('x m omega')
# strPot=["m / 2 * omega ** 2 * x ** 2"]
# symPot=[sp.parsing.sympy_parser.parse_expr(item) for item in strPot]
# funcs = [sp.lambdify(x + m + omega, f) for f in symPot]
# print(funcs)

