#Analytics
import numpy as np
def anal_q(amp,omega,m,phi,t,dt):
    return amp*np.sin(omega*np.arange(0,t,dt)+phi)

def anal_p(amp,omega,m,phi,t,dt):
    return m*omega*amp*np.cos(omega*np.arange(0,t,dt)+phi)

#compute hamiltonian
def get_pe(q,omega,m):#H-dep
    return m/2*omega**2*q**2
def get_ke(p,m):
    return p**2/2/m
def get_ham(p,q,omega,m):
    return get_pe(q,omega,m)+get_ke(p,m)

def get_pe_Morse(q,alpha,de):
    return de*(1-np.exp(-alpha*q))

def get_ham_Morse(p,q,alpha,de,m):
    return get_pe_Morse(q,alpha,de)+get_ke(p,m)
