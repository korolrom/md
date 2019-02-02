import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
#input parameters:
m=1
omega=1
dt=1e-1
t=10
#compute hamiltonian
def pot(q):#H-dep
    return m/2*omega**2*q**2
def ke(p):
    return p**2/2/m
def ham(p,q):
    return pot(q)+ke(p)
#future sampler
def sampleq():
    return np.random.rand()
def samplep():
    return np.random.rand()+1
    #p is shifted so there are no particles in the classically forbidden region
p0=samplep()
q0=sampleq()
E=ham(p0,q0)
x=sp.Symbol('x')
print(pot(x))
print(-sp.diff(pot(x)))
#Analytical solution for harmonic oscillator
A=np.sqrt(2*E/m)/omega
phi=np.arcsin(q0/A)

def analq(A,omega,phi,t,dt):
    return A*np.sin(omega*np.arange(0,t,dt)+phi)

def analp(A,omega,phi,t,dt):
    return m*omega*A*np.cos(omega*np.arange(0,t,dt)+phi)

#things to plot for analytics: analPos, analVel, pot(analPos), ke(analVel), ham(analPos, analVel)

#Numerics
def Euler(p,q,t,dt):
    eulerPQ=np.zeros((2,int(t/dt)))
    for i in np.arange(0,int(t/dt)):
        eulerPQ[:,i]=[p, q]
        #print(eulerPQ[:,i])
        a=-omega**2*q #H-dep
        q=q+p/m*dt+a/2*dt**2 # position from Newton's eqtns
        p=p+m*a*dt # momentum from Newton's equations
    return eulerPQ

def Verlet(p,q,t,dt):
    verletPQ=np.zeros((2,int(t/dt)))
    a=-omega**2*q #H-dep
    q_1=q-p/m*dt-a/2*dt**2#get position from previous step
    for i in np.arange(0,int(t/dt)):
        verletPQ[:,i]=[p, q]
        # print(verletPQ[:,i])
        a=-omega**2*q #H-dep, force is evaluated at current timestep
        q1=2*q-q_1+a*dt**2 #new position at the next timestep
        p=(q1-q_1)*m/dt/2 #compute the momentum across this timestep
        q_1=q #saving the old position for propagation in the next timestep 
        q=q1 #saving the new position at the end of this timespep
    return verletPQ


def Verlet(p,q,t,dt):
    verletPQ=np.zeros((2,int(t/dt)))
    a=-omega**2*q #H-dep
    q_1=q-p/m*dt-a/2*dt**2#get position from previous step
    for i in np.arange(0,int(t/dt)):
        verletPQ[:,i]=[p, q]
        # print(verletPQ[:,i])
        a=-omega**2*q #H-dep, force is evaluated at current timestep
        q1=2*q-q_1+a*dt**2 #new position at the next timestep
        p=(q1-q_1)*m/dt/2 #compute the momentum across this timestep
        q_1=q #saving the old position for propagation in the next timestep 
        q=q1 #saving the new position at the end of this timespep
    return verletPQ
#print(p0,q0)

plt.plot(np.arange(0,t,dt),analq(A,omega,phi,t,dt),'r-')
plt.plot(np.arange(0,t,dt),analp(A,omega,phi,t,dt),'r--')
# print(Euler(p0,q0,t,dt))
plt.plot(np.arange(0,t,dt),Euler(p0,q0,t,dt)[0,:],'k--')
plt.plot(np.arange(0,t,dt),Euler(p0,q0,t,dt)[1,:],'k')
#print(Verlet(p0,q0,t,dt))
plt.plot(np.arange(0,t,dt),Verlet(p0,q0,t,dt)[0,:],'b--')
plt.plot(np.arange(0,t,dt),Verlet(p0,q0,t,dt)[1,:],'b')
plt.show()
