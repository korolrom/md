import numpy as np

def Euler(momentum,position,time,time_step,mass,new_a):
    """Straightforward Euler integrator. Gives the trajectory of a particle by propagating Classical EoM 
    with first order accuracy"""
    Euler = np.zeros((2,int(time/time_step)))
    for i in np.arange(0,int(time/time_step)):
        Euler[:,i] = [momentum, position]
        #print(Euler[:,i])
        acceleration = new_a(position) #Depends on Hamiltonian
        position   += momentum/mass*time_step + acceleration/2*time_step**2 # position from Newton's eqtns
        momentum   += mass*acceleration*time_step # momentum from Newton's equations
    return Euler

def Verlet(momentum,position,time,time_step,mass,new_a):
    """Verlet integrator [Verlet 1967] Gives the trajectory of a particle by propagating Classical EoM 
    The accuracy of position is order dt^4, velocity is dt^3. 
    Disatvantage: possible loss of numerical precision, since to get next position you add a small quantity (propto dt^2)
    Advantage: synplectic integrator (area and energy preserving)"""
    Verlet = np.zeros((2,int(time/time_step)))
    acceleration = new_a(position) #Depends on Hamiltonian
    prev_position = position - momentum/mass*time_step - acceleration/2*time_step**2#get position from previous step
    for i in np.arange(0,int(time/time_step)):
        Verlet[:,i] = [momentum, position]
        # print(Verlet[:,i])
        acceleration = new_a(position) #Depends on Hamiltonian, force is evaluated at current timestep
        next_position = 2*position - prev_position + acceleration*time_step**2 #new position at the next timestep
        momentum = (next_position - prev_position)*mass/time_step/2 #compute the momentum across this timestep
        prev_position = position #saving the old position for propagation in the next timestep 
        position = next_position #saving the new position at the end of this timespep
    return Verlet


def vVerlet(momentum,position,time,time_step,mass,new_a):
    """velocity Verlet integrator [Swope 1982]. Gives the trajectory of a particle by propagating Classical EoM 
    using velocity-Verlet version. It combines advantages of Verlet and Leapfrog algorithms."""
    VVerlet = np.zeros((2,int(time/time_step)))
    acceleration = new_a(position)
    prev_position = position - momentum/mass*time_step - acceleration/2*time_step**2 # get position from previous step
    prev_acceleration = new_a(prev_position)
    for i in np.arange(0,int(time/time_step)):
        VVerlet[:,i] = [momentum, position]
        # print(vVerlet[:,i])
        acceleration = new_a(position) #Depends on Hamiltonian, force is evaluated at current timestep
        position += momentum/mass*time_step + acceleration/2*time_step**2 #naive position
        momentum += (acceleration + prev_acceleration)*mass*time_step/2 #not naive momentum
        prev_acceleration = acceleration #saving the old acceleration(force) for propagation in the next timestep 
    return VVerlet

def leapFrog(momentum,position,time,time_step,mass,new_a): #not tested yet
    """Lepfrog integrator [Hockney, 1970] is a variation of Verlet algorithm. (formally equivalet to Verlet)
    Gives the trajectory of a particle by propagating Classical EoM.
    Advantage: numerically better (no difference of large and small numbers)
    Disatvantage: it calculates velocity at midpoints of the timesteps making it hard to evaluate total energy"""
    LeapFrog = np.zeros((2,int(time/time_step)))
    acceleration = new_a(position)
    momentum_prev = momentum - acceleration/2*time_step #get the momentum half a time step ago
    for i in np.arange(0,int(time/time_step)):
        # print(LeapFrog[:,i])
        momentum = momentum_prev+acceleration*time_step*mass # momentum for the half timestep later
        position += momentum/mass*time_step
        acceleration = new_a(position) #Depends on Hamiltonian, force is evaluated at current timestep
        LeapFrog[:,i] = [(momentum+momentum_prev)/2, position]# note that the velocity is the average between half-step back and half step forward
        momentum_prev=momentum
    return LeapFrog

def Beeman(momentum,position,time,time_step,mass,new_a):
    """Gives the trajectory of a particle by propagating Classical EoM using Beeman integrator [Beeman 1976]
    It is related to Verlet scheme, but uses a more precise calculation for velocity. Therefore, its energy conservation is usually better (KE is more accurate)"""
    Beeman = np.zeros((2,int(time/time_step)))
    acceleration = new_a(position)
    prev_position = position - momentum/mass*time_step - acceleration/2*time_step**2#get position from previous step
    prev_acceleration = new_a(prev_position)
    for i in np.arange(0,int(time/time_step)):
        Beeman[:,i] = [momentum, position]
        # print(BeemanPQ[:,i])
        acceleration = new_a(position) #Depends on Hamiltonian, force is evaluated at current timestep
        position += momentum/mass*time_step + (4*acceleration - prev_acceleration)*time_step**2/6 #naive position
        next_acceleration = new_a(position)
        momentum += (2*next_acceleration + 5*acceleration - prev_acceleration)*mass*time_step/6 #not naive momentum
        prev_acceleration = acceleration #saving the old acceleration(force) for propagation in the next timestep
        acceleration = next_acceleration 
    return Beeman

def PredCorr(momentum,position,time,time_step,mass,new_a,new_j):# not tested yet, not sure about the corrections
    """Predictor-corrector algorithm [Gear, 1971]. Gives the trajectory of a particle by propagating Classical EoM. 
    The idea is to compare the accelerations at t+dt to the values predicted by Taylor expansion around t and to correct the EoM for position and momentum accordingly."""
    PredCorr = np.zeros((2,int(time/time_step)))
    jerk = new_j(position)
    acceleration = new_a(position)
    for i in np.arange(0,int(time/time_step)):
        PredCorr[:,i] = [momentum, position]
        position += momentum*time_step/mass + acceleration*time_step**2/2 + jerk*time_step**3/6 #Taylor expansions: 7.2-7.4
        momentum += mass*(acceleration*time_step+jerk*time_step**2/2)
        accelerationDifference = (acceleration + jerk*time_step) - new_a(position) # a^C(t+dt) - a(t+dt) i.e. the difference btw Taylor expansion and evaluation of forces
        position -=accelerationDifference/6 #eq. 7.23 corrections
        momentum -= mass*5*accelerationDifference/6
        acceleration -= 2*accelerationDifference
        jerk -= 2*accelerationDifference
    return PredCorr

def PredCorrRahman(momentum,position,time,time_step,mass,new_a):# not tested yet, not sure about the corrections
    """Predictor-corrector algorithm [Raahman, 1964]. Gives the trajectory of a particle by propagating Classical EoM. 
    The idea is to compare the accelerations at t+dt to the values predicted by Taylor expansion around t and to correct the EoM for position and momentum accordingly."""
    PredCorrRahman = np.zeros((2,int(time/time_step)))
    acceleration = new_a(position)
    for i in np.arange(0,int(time/time_step)):
        PredCorrRahman[:,i] = [momentum, position]
        position += momentum*time_step/mass
        new_acceleration=new_a(position)
        new_momentum = momentum + mass*time_step*(new_acceleration+acceleration)/2 #eq. 7.28
        position += (momentum+new_momentum)*time_step/mass # eq. 7.29 (corrected position)
    return PredCorrRahman

