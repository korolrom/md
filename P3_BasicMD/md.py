"""
#!/usr/bin/env python #use in Linux
This is an MD simulation that integrates EoM for the SHO and Morse potentials
and plots the resulting trajectories, together with a check for the global and local
energy conservation
"""
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
import num_int as num
import analytic as ana
#input parameters:
m=1.
dt=1e-2
t=10.
morseTag=0
if morseTag: #Morse
    alpha=1.
    de=1.
    omega=alpha**2*de
else:
    omega=1.

#future sampler
def sample_q():
    """This will one day be a great sampler"""
    return np.random.rand()
def sample_p():
    """This will also one day be a great sampler"""
    return np.random.rand()+1
    #p is shifted so there are no particles in the classically forbidden region
def new_a_SHO(position):
    """Gives an acceleration kick (i.e. force kick divided by mass for the SHO potential)"""
    return -omega**2*position

def new_a_Morse(position):
    """Gives an acceleration kick (i.e. force kick divided by mass for the Morse potential)"""
    return 2*alpha*de*(np.exp(-alpha*position)-1)

p0=sample_p()
q0=sample_q()
e=ana.pe_SHO(q0,omega,m)+ana.ke(p0,m)

#Analytical solution for harmonic oscillator
amp=np.sqrt(2*e/m)/omega
phi=np.arcsin(q0/amp)
AnalPQ=np.array([ana.p_SHO(amp,omega,m,phi,t,dt), ana.q_SHO(amp,omega,m,phi,t,dt)])
AnalE=np.array([ana.ke(AnalPQ[0,:],m),ana.pe_SHO(AnalPQ[1,:],omega,m)])

#Numerical integration of EoM to get trajectories
if morseTag: #Morse
    PQ=np.array([AnalPQ, num.Euler(p0,q0,t,dt,m,new_a_Morse),num.Verlet(p0,q0,t,dt,m,new_a_Morse),
            num.vVerlet(p0,q0,t,dt,m,new_a_Morse),num.Beeman(p0,q0,t,dt,m,new_a_Morse)])
else:
    PQ=np.array([AnalPQ, num.Euler(p0,q0,t,dt,m,new_a_SHO),num.Verlet(p0,q0,t,dt,m,new_a_SHO),
            num.vVerlet(p0,q0,t,dt,m,new_a_SHO),num.Beeman(p0,q0,t,dt,m,new_a_SHO)])

#Get the energies on the trajectory
if morseTag: #Morse
    E=np.swapaxes(np.array([ana.ke(PQ[:,0,:],m),ana.pe_Morse(PQ[:,1,:],alpha,de)]),1,0)
else:
    E=np.swapaxes(np.array([ana.ke(PQ[:,0,:],m),ana.pe_SHO(PQ[:,1,:],omega,m)]),1,0)

#Compute global energy conservation and the standard deviation in E
Ebar=np.cumsum(np.sum(E,axis=1),axis=-1)/(np.arange(int(t/dt))+1)
Esigma=np.cumsum((np.sum(E,axis=1)-Ebar)**2,axis=-1)/(np.arange(int(t/dt))+1)

# Plots of p and q with time, energies with time
n=np.size(E,axis=0)
ColorList=[char for char in "rgbkmcy"]
MarkerList=[char for char in ".o^s+xD"]
LineList=['-','--',':','-.']
LineStyle=[[line+color for line in LineList] for color in ColorList]
#LineStyle=[['r-','r--','r:','r-.'],['g-','g--','g:','g-.'],['b-','b--','b:','b-.'],['k-','k--','k:','k-.'],['m-','m--','m:','m-.'],['c-','c--','c:','c-.'],['y-','y--','y:','y-.']]
ScatterStyle=[[marker+color for marker in MarkerList] for color in ColorList]
#ScatterStyle=[['r.','ro','r^','rs','r+','rx','rD'],['g.','go','g^','gs','g+','gx','gD'],['b-.','b.','bo','b^','bs','b+','bx','bD'],['k.','ko','k^','ks','k+','kx','kD'],['m.','mo','m^','ms','m+','mx','mD'],['c.','co','c^','cs','c+','cx','cD'],['y.','yo','y^','ys','y+','yx','yD']]
Curves=['SHO Analytic','Euler','Verlet','vVerlet','Beeman']
font=12
fig = plt.figure()
subj=1
for item in [PQ,E]:
    axes_1 = fig.add_subplot(1, 2, subj)
    for i in range(n):#[0,4]:
        plt.plot(np.arange(0,t,dt),item[i,0,:],LineStyle[i][0],label=Curves[i])
        plt.plot(np.arange(0,t,dt),item[i,1,:],LineStyle[i][1])
        if subj==1:
            plt.legend()
            plt.text(0.01, 0.01, 'p (solid line)\nq (dashed line)', fontsize=font,
            horizontalalignment='left', verticalalignment='bottom', transform = axes_1.transAxes)
            plt.ylabel("Phase space coordinates", fontsize=font)
        else:
            plt.plot(np.arange(0,t,dt),np.sum(item[i,:,:],axis=0),LineStyle[i][2])
            plt.plot(np.arange(0,t,dt),Ebar[i,:],LineStyle[i][3])
            plt.plot(np.arange(0,t,dt),np.sqrt(Esigma[i,:]),LineStyle[i][3])
            plt.ylabel("Energy", fontsize=font)
            plt.text(0.01, 0.01, "KE (solid), PE (dashed)\nlocal E (dotted), global E (dash-dotted)\nstandard deviation of E (dots)",
            fontsize=font, horizontalalignment='left', verticalalignment='bottom', transform = axes_1.transAxes)

    plt.xlabel("time, t", fontsize=font)
    subj=subj+1
plt.tight_layout()
#plt.savefig("fig_wo_cpcp.pdf")
plt.show()

#print(np.shape(E))
#print(np.shape(PQ))
#print(p0,q0)


# x=sp.symbols('x m omega')
# strPot=["m / 2 * omega ** 2 * x ** 2"]
# symPot=[sp.parsing.sympy_parser.parse_expr(item) for item in strPot]
# funcs = [sp.lambdify(x + m + omega, f) for f in symPot]
# print(funcs)

