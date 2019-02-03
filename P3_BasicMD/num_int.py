import numpy as np
#Numerics
def Euler(p,q,t,dt,m,da):
    Euler=np.zeros((2,int(t/dt)))
    for i in np.arange(0,int(t/dt)):
        Euler[:,i]=[p, q]
        #print(Euler[:,i])
        a=da(q) #H-dep
        q=q+p/m*dt+a/2*dt**2 # position from Newton's eqtns
        p=p+m*a*dt # momentum from Newton's equations
    return Euler

def Verlet(p,q,t,dt,m,da):
    Verlet=np.zeros((2,int(t/dt)))
    a=da(q) #H-dep
    q_1=q-p/m*dt-a/2*dt**2#get position from previous step
    for i in np.arange(0,int(t/dt)):
        Verlet[:,i]=[p, q]
        # print(Verlet[:,i])
        a=da(q) #H-dep, force is evaluated at current timestep
        q1=2*q-q_1+a*dt**2 #new position at the next timestep
        p=(q1-q_1)*m/dt/2 #compute the momentum across this timestep
        q_1=q #saving the old position for propagation in the next timestep 
        q=q1 #saving the new position at the end of this timespep
    return Verlet


def vVerlet(p,q,t,dt,m,da):#velocity Verlet, same result as Verlet
    VVerlet=np.zeros((2,int(t/dt)))
    a=da(q)
    q_1=q-p/m*dt-a/2*dt**2#get position from previous step
    a_1=da(q_1)
    for i in np.arange(0,int(t/dt)):
        VVerlet[:,i]=[p, q]
        # print(vVerlet[:,i])
        a=da(q) #H-dep, force is evaluated at current timestep
        q=q+p/m*dt+a*dt**2 #naive position
        p=p+(a+a_1)*m*dt/2 #not naive momentum
        a_1=a #saving the old acceleration(force) for propagation in the next timestep 
    return VVerlet

def Beeman(p,q,t,dt,m,da):#Beeman: q's are the same as Verlet, v's are more accurate but with no t-reversal symmetry;
    Beeman=np.zeros((2,int(t/dt)))
    a=da(q)
    q_1=q-p/m*dt-a/2*dt**2#get position from previous step
    a_1=da(q_1)
    for i in np.arange(0,int(t/dt)):
        Beeman[:,i]=[p, q]
        # print(BeemanPQ[:,i])
        a=da(q) #H-dep, force is evaluated at current timestep
        q=q+p/m*dt+(4*a-a_1)*dt**2/6 #naive position
        a1=da(q)
        p=p+(2*a1+5*a-a_1)*m*dt/6 #not naive momentum
        a_1=a #saving the old acceleration(force) for propagation in the next timestep
        a=a1 
    return Beeman