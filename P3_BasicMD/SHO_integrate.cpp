/* First attempt of the MD code from the reference in P3 folder. date:2018_12_02*/

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <numeric> 

using namespace std;

// void Harmonic();
// void Euler(double dt);
// void Verlet(double dt);
// void VVerlet(double dt);
// void Beeman(double dt);

//initialization
const int N=1, m=1, omega=8.944, nsteps=1e3, npoints=1e3; //number of particles, mass, frequency
const double x0=0.2, p0=0, t0=0, f0=0, a0=0, dt=1e-2, dx=1e-2;//initialization, one of p0, x0 is zero!

// int method=0; //0 for Euler, 1 for Verlet, 2 for vVerlet;
double etot, amp, phi; // analytic parameters
double Pos[nsteps], Vel[nsteps], Ekin[nsteps], Epot[nsteps];
double t, f, v, x, a, a_1, x_1, tmp, f1, a1, ek, ep, et; //simulation parameters
double Pot_pos[npoints], Pot_SHO[npoints];//evaluate potential at each point
double Time[nsteps], Position[nsteps], Velocity[nsteps], Energy[nsteps], KinEnergy[nsteps], PotEnergy[nsteps];
double Ecum[nsteps], Esquare[nsteps], Eerrorsum[nsteps], Eerror[nsteps], Eavg[nsteps];
int main(){
    
    //Analytics
    etot = m*(omega*x0)*(omega*x0)/2 + p0*p0/2/m;
    amp = sqrt(2*etot/m)/omega;
    phi = acos(x0/amp);
    //cout<<phi<<" "<<etot<<" "<<amp<<" "<<omega;   
    for (int i = 0; i < nsteps; i++) {
        Pos[i] = amp*cos(omega*i*dt+phi);//H-dependent
        Vel[i] = -omega*amp*sin(omega*i*dt+phi);//H-dependent
        Epot[i] = m*omega*omega*Pos[i]*Pos[i]/2;//H-dependent
        Ekin[i] = m*Vel[i]*Vel[i]/2;
    }
    for (int i = 0; i < npoints; i++) {
        Pot_pos[i] = (i-nsteps/2)*dx;
        Pot_SHO [i]= (m*omega*omega*((i-nsteps/2)*dx)*((i-nsteps/2)*dx))/2; 
        //compute potential on the symmetric interval around 0
    }
    
    //Numerics
    for (int method=0; method<4; method++) {
    v = p0/m; x = x0; t = t0; f = f0; a = a0;
    x_1 = x0-p0/m*dt;//for Verlet method it is necessary to propagate one step back
    for(int i = 0; i < nsteps; i++) {
      Position[i]=x; Velocity[i]=v; Time[i]=t; //save as the arrays
      KinEnergy[i] = m*v*v/2; PotEnergy[i] = m*omega*omega*x*x/2; Energy[i] = KinEnergy[i]+PotEnergy[i];//energy
      a_1 = a; // not updated acceleration (from previous Time step)
      f = -m*omega*x;//force
      a = f/m;//acceleration
      t += dt;//increment time
    switch(method) {
        case 0://Euler
            x += v*dt+a/2*dt*dt;//naive position
            v += a*dt;//naive velocity
            //   printf("%f %f %f\n", t, x, v);        
        break;
        case 1://Verlet
            tmp = x;//get position from previous step
            x = 2*x - x_1 + a*dt*dt;//new position (Verlet)
            v = (x-x_1)/2/dt;//one does not need to compute v to do the propagation, just for the final output
            x_1 = tmp;//save position from this step for the next
            break;
        case 2://velocity Verlet, same result as Verlet
            x += v*dt+a*dt*dt/2;//naive Position
            v += dt/2*(a+a_1);//not naive velocity
            break;
        case 3://Beeman: x's are the same as Verlet, v's are more accurate but with no t-reversal symmetry;
            x += v*dt+(4*a-a_1)*dt*dt/6;//naive position
            f1 = -m*omega*x;// new force
            a1 = f1/m;//new acceleration
            v += dt/6*(2*a1+5*a-a_1);//not naive velocity
        break;
  }
  }
  Ecum[0]=Energy[0]; Eavg[0]=Energy[0]; Esquare[0]=0; Eerrorsum[0]=0; Eerror[0]=0;
for (int i = 1; i<nsteps; i++){
    Ecum[i]=Energy[i]+Ecum[i-1];
    Eavg[i]=Ecum[i]/(i+1);
    Esquare[i]=(Energy[i]-Eavg[i])*(Energy[i]-Eavg[i]);
    Eerrorsum[i]=Esquare[i]+Eerrorsum[i-1];
    Eerror[i]=Eerrorsum[i]/(i+1);
}

  //Output
    ofstream outputfile; //declare an output file (a variable)
	outputfile.open("out_Exact.txt", ios::out);	// creates the file if it wasn't there and replaces if it is there
    for(int i = 0; i < nsteps; i++) {
 	outputfile<<Time[i]<<" "<<Pos[i]<<" "<<Vel[i]<<" "<<" "<<Epot[i]<<" "<<Ekin[i]<<" "<<etot<<endl;
	}
    outputfile.close();//it's a good idea to always close the file once you are done

	outputfile.open("out_Potential.txt", ios::out);	// creates the file if it wasn't there and replaces if it is there
    for(int i = 0; i < nsteps; i++) {
 	outputfile<<Pot_pos[i]<<" "<<Pot_SHO[i]<<" "<<etot<<endl;
	}
    outputfile.close();//it's a good idea to always close the file once you are done
    
    switch(method) {
    case 0:
    outputfile.open("out_Euler.txt", ios::out);	// creates the file if it wasn't there and replaces if it is there
	//you can also append using ios::app
    break;
    case 1: outputfile.open("out_Verlet.txt", ios::out); break;
    case 2: outputfile.open("out_vVerlet.txt", ios::out); break;
    case 3: outputfile.open("out_Beeman.txt", ios::out); break;
    }
    for(int i = 0; i < nsteps; i++) {
 	outputfile<<Time[i]<<" "<<Position[i]<<" "<<Velocity[i]<<" "<<PotEnergy[i]<<" "<<KinEnergy[i]<<" "<<Eavg[i]<<" "<<Eerror[i]<<endl;
    }
	outputfile.close();//it's a good idea to always close the file once you are done
    }
	
  }