/* First attempt of the MD code from the reference in P3 folder. date:2018_12_02*/
//Simple model of interacting Argon atoms
/*
Consider N atoms of argon which interact via a Lennard-Jones potential:
V(r)=4*\epsilon*[(\frac{\sigma}{r})^{12}-(\frac{\sigma}{r})^6]
where \epsilon = 1.65×10−21 J, and \sigma = 3.4 × 10−10 m is the value of r at which the energy is zero. 
The 1/r^12 term represents a repulsive hard core interaction between the argon atoms. The 1/r^6 term represents an
attractive dipole-dipole (van der Waals) interaction between the non-polar atoms. The potential has its
minimum at V (2^{1/6})=-\epsilon
One can find the force, which is just the negative derivative of the potential.
We will choose units of mass, length and energy so that m = 1, \sigma = 1, and \epsilon = 1. 
The unit of time in this system is given by \sqrt{\frac{m\sigma^2}{\epsilon}} = 2.2*10^{-12} s, 
so natural time scale for the dynamics of this system is a few picoseconds
*/
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;
const int N = 64;  // number of particles //parameter
double r[N][3];    // positions
double v[N][3];    // velocities
double a[N][3];    // accelerations
/*
We next need to set the initial positions and velocities of the particles. This is actually a complicated
problem! Because the system can be simulated only for a few nanoseconds, the starting configuration must
be very close to equilibrium to get good results. For a dense system, the atoms are usually placed at the
vertices of a face-centered cubic lattice, which is tends to minimize the potential energy. The atoms are
also given random velocities to approximate the desired temerature.
In this preliminary program we will put the system in a cubical volume of side L and place the particles at
the vertices of a simple cubic lattice:
*/
double L = 10;     // linear size of cubical volume
double vMax = 0.1; // maximum initial velocity component

void initialize()
{
    int n = int(ceil(pow(N, 1.0 / 3))); // number of atoms in each direction
    double a = L / n;                   // lattice spacing
    int pp = 0;                          // particles placed so far
    // initialize positions
    for (int x = 0; x < n; x++)
        for (int y = 0; y < n; y++)
            for (int z = 0; z < n; z++)
            {
                if (pp < N)
                {
                    r[pp][0] = (x + 0.5) * a;
                    r[pp][1] = (y + 0.5) * a;
                    r[pp][2] = (z + 0.5) * a;
                }
                ++pp;
            }
    // initialize velocities
    for (int pp = 0; pp < N; pp++)
        for (int i = 0; i < 3; i++)
            v[pp][i] = vMax * (2 * rand() / double(RAND_MAX) - 1); //gives number btw -vMax and vMax
}
/*
The following function computes the accelerations of the particles from their current positions 
using Newton's equations of motion:
*/
void computeAccelerations()
{
    for (int i = 0; i < N; i++)
        for (int kk = 0; kk < 3; kk++)
            a[i][kk] = 0;            // set all accelerations to zero
    for (int i = 0; i < N - 1; i++) // loop over all distinct pairs i,j
        for (int j = i + 1; j < N; j++)
        {
            double rij[3]; // position of i relative to j
            double rSqd = 0;
            for (int kk = 0; kk < 3; kk++)
            {
                rij[kk] = r[i][kk] - r[j][kk];
                rSqd += rij[kk] * rij[kk];
            }
            double f = 24 * (2 * pow(rSqd, -7) - pow(rSqd, -4));//from LJ potential
            for (int kk = 0; kk < 3; kk++)
            {
                a[i][kk] += rij[kk] * f;
                a[j][kk] -= rij[kk] * f;
            }
        }
}
/*
It can be shown that the errors in this algorithm are of O(dt^4), and that it is very stable in MD applications
and in particular conserves energy very well.
The following function advances the positions and velocities of the particles by one time step:
*/
void velocityVerlet(double dt) {
computeAccelerations();
for (int i = 0; i < N; i++)
for (int kk = 0; kk < 3; kk++) {
r[i][kk] += v[i][kk] * dt + 0.5 * a[i][kk] * dt * dt;
v[i][kk] += 0.5 * a[i][kk] * dt;
}
computeAccelerations();
for (int i = 0; i < N; i++)
for (int kk = 0; kk < 3; kk++)
v[i][kk] += 0.5 * a[i][kk] * dt;
}
//The instantaneous temperature
/*
This is a simulation in which the number of particles N and the volume L^3 of the system are fixed. Because
the Lennard-Jones force is conservative, the total energy of the system is also constant.
If the system is in thermal equilibrium, then Boltzmann’s Equipartition Theorem relates the absolute temperature
T to the kinetic energy:
3(N − 1)*(0.5k_BT)=<KE>
Here the angle brackets h...i represent a thermal ensemble average. The factor 3(N − 1) is the number of
internal translational degrees of freedom which contribute to thermal motion: the motion of the center of
mass of the system does not represent thermal energy!
*/
double instantaneousTemperature() {
double sum = 0;
for (int i = 0; i < N; i++)
for (int k = 0; k < 3; k++)
sum += v[i][k] * v[i][k];
return sum / (3 * (N - 1));
}
//Finally, here is the main function which steers the simulation
int main() {
initialize();
double dt = 0.01;
ofstream file;
file.open("T.txt", ios::out);
for (int i = 0; i < 1000; i++) {
velocityVerlet(dt);
file <<i*dt <<" "<<instantaneousTemperature() << "\n";
}
file.close();
}
/*
The instantaneous temperature is approximately constant for around one or two time units, 
and then it starts increasing with fluctuations.
There are several ways in which this simple program needs to be improved:
1)The volume is not really constant because the particles can move out of it! We need to impose suitable
boundary conditions, for example periodic boundary conditions.
2)The initial positions and velocities need to be chosen more carefully. We will place the particles on a
face-centered cubic lattice, and use a Maxwell-Boltzmann distribution for the velocities.
3)The system needs to be allowed to come to thermal equilibrium at the desired temperature.
4)Thermal averages of various quantities need to be measured.
*/
//Thus we turn to a more sophisticated, adv_md.cpp next