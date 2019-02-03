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
# x=sp.Symbol('x')
# print(pot(x))
# print(-sp.diff(pot(x)))
#Analytical solution for harmonic oscillator
A=np.sqrt(2*E/m)/omega
phi=np.arcsin(q0/A)

def analq(A,omega,phi,t,dt):
    return A*np.sin(omega*np.arange(0,t,dt)+phi)

def analp(A,omega,phi,t,dt):
    return m*omega*A*np.cos(omega*np.arange(0,t,dt)+phi)

#things to plot for analytics: analPos, analVel, pot(analPos), ke(analVel), ham(analPos, analVel)

print(p0,q0)

plt.plot(np.arange(0,t,dt),analq(A,omega,phi,t,dt),'r-')
plt.plot(np.arange(0,t,dt),analp(A,omega,phi,t,dt),'r--')
print(Euler(p0,q0,t,dt))[:2,:2]
plt.plot(np.arange(0,t,dt),Euler(p0,q0,t,dt)[0,:],'k--')
plt.plot(np.arange(0,t,dt),Euler(p0,q0,t,dt)[1,:],'k')
#print(Verlet(p0,q0,t,dt))
plt.plot(np.arange(0,t,dt),Verlet(p0,q0,t,dt)[0,:],'b--')
plt.plot(np.arange(0,t,dt),Verlet(p0,q0,t,dt)[1,:],'b')
plt.plot(np.arange(0,t,dt),vVerlet(p0,q0,t,dt)[0,:],'c--')
plt.plot(np.arange(0,t,dt),vVerlet(p0,q0,t,dt)[1,:],'c')
plt.plot(np.arange(0,t,dt),Beeman(p0,q0,t,dt)[0,:],'m--')
plt.plot(np.arange(0,t,dt),Beeman(p0,q0,t,dt)[1,:],'m')
plt.show()


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
    VerletPQ=np.zeros((2,int(t/dt)))
    a=-omega**2*q #H-dep
    q_1=q-p/m*dt-a/2*dt**2#get position from previous step
    for i in np.arange(0,int(t/dt)):
        VerletPQ[:,i]=[p, q]
        # print(VerletPQ[:,i])
        a=-omega**2*q #H-dep, force is evaluated at current timestep
        q1=2*q-q_1+a*dt**2 #new position at the next timestep
        p=(q1-q_1)*m/dt/2 #compute the momentum across this timestep
        q_1=q #saving the old position for propagation in the next timestep 
        q=q1 #saving the new position at the end of this timespep
    return VerletPQ


def vVerlet(p,q,t,dt):#velocity Verlet, same result as Verlet
    vVerletPQ=np.zeros((2,int(t/dt)))
    a=-omega**2*q
    q_1=q-p/m*dt-a/2*dt**2#get position from previous step
    a_1=-omega**2*q_1
    for i in np.arange(0,int(t/dt)):
        vVerletPQ[:,i]=[p, q]
        # print(vVerletPQ[:,i])
        a=-omega**2*q #H-dep, force is evaluated at current timestep
        q=q+p/m*dt+a*dt**2 #naive position
        p=p+(a+a_1)*m*dt/2 #not naive momentum
        a_1=a #saving the old acceleration(force) for propagation in the next timestep 
    return vVerletPQ

def Beeman(p,q,t,dt):#Beeman: q's are the same as Verlet, v's are more accurate but with no t-reversal symmetry;
    BeemanPQ=np.zeros((2,int(t/dt)))
    a=-omega**2*q
    q_1=q-p/m*dt-a/2*dt**2#get position from previous step
    a_1=-omega**2*q_1
    for i in np.arange(0,int(t/dt)):
        BeemanPQ[:,i]=[p, q]
        # print(BeemanPQ[:,i])
        a=-omega**2*q #H-dep, force is evaluated at current timestep
        q=q+p/m*dt+(4*a-a_1)*dt**2/6 #naive position
        a1=-omega**2*q
        p=p+(2*a1+5*a-a_1)*m*dt/6 #not naive momentum
        a_1=a #saving the old acceleration(force) for propagation in the next timestep
        a=a1 
    return BeemanPQ