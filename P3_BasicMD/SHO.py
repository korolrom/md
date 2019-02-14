import numpy as np
def anal_q(amplitude,omega,mass,phi,time,time_step):
    """Returns the position of SHO as np.array at each timestep given 
    the initial condition (phi) as an angle in radians, time and SHO parameters"""
    return amplitude*np.sin(omega*np.arange(0,time,time_step)+phi)

def anal_p(amplitude,omega,mass,phi,time,time_step):
    """Returns the momentum of SHO as np.array at each timestep given 
    the initial condition (phi) as an angle in radians, time and SHO parameters"""
    return mass*omega*amplitude*np.cos(omega*np.arange(0,time,time_step)+phi)

#compute hamiltonian
def get_pe_SHO(position,omega,mass):#H-dep
    """Returns the potential energy of SHO as np.array at each timestep given 
    the initial condition (phi) as an angle in radians, time and SHO parameters"""
    return mass/2*omega**2*position**2

def get_ke(momentum,mass):
    """Returns the kinetic energy of SHO as np.array at each timestep given 
    the initial condition (phi) as an angle in radians, time and SHO parameters"""
    return momentum**2/2/mass
def get_ham(get_pe,momentum,position,omega,mass):
    return get_pe(position,omega,mass)+get_ke(momentum,mass)

def get_pe_Morse(position,alpha,De):
    return De*(1-np.exp(-alpha*position))