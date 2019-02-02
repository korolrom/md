/* Second attempt of the MD code from the reference in P3 folder. date:2018_12_21*/

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

int N = 64;       // number of particles
double rho = 1.0; // density (number per unit volume)
double T = 1.0;   // temperature
// function declarations
void initialize();        // allocates memory, calls following 2 functions
void initPositions();     // places particles on an fcc lattice
void initVelocities();    // initial Maxwell-Boltzmann velocity distribution
void rescaleVelocities(); // adjust the instanteous temperature to T
double gasdev();          // Gaussian distributed random numbers

//We will allocate particle arrays dynamically rather than statically
double **r; // positions
//** is a double pointer, a pointer to a pointer. It is frequently used for arrays.
double **v; // velocities
double **a; // accelerations
void initialize()
{
    r = new double *[N]; //new is an operator that allocates storage space
    v = new double *[N];
    a = new double *[N];
    for (int i = 0; i < N; i++)
    {
        r[i] = new double[3];
        v[i] = new double[3];
        a[i] = new double[3];
    }
    initPositions();
    initVelocities();
}
/*
The minimum energy configuration of this Lennard-Jones system is an fcc lattice. This has 4 lattice sites in
each conventional cubic unit cell. If the number of atoms N = 4M3, where M = 1, 2, 3, ..., then the atoms
can fill a cubical volume. So MD simulations are usually done with 32 = 4×2^3, 108 = 4×3^3, 256, 500, 864, ...
atoms.
*/
double L; // linear size of cubical volume

void initPositions()
{
    // compute side of cube from number of particles and number density
    L = pow(N / rho, 1.0 / 3);
    // find M large enough to fit N atoms on an fcc lattice
    int M = 1;
    while (4 * M * M * M < N)
        ++M;
    double a = L / M; // lattice constant of conventional cell
    /*
The 4 atoms with positions in units of a: (0, 0, 0) (0.5, 0.5, 0) (0.5, 0, 0.5) (0, 0.5, 0.5) form a basis.
In the following code we shift the basis by (0.5, 0.5, 0.5) 
so the all atoms are inside the volume and none are on its boundaries.
*/
    // 4 atomic positions in fcc unit cell
    double xCell[4] = {0.25, 0.75, 0.75, 0.25};
    double yCell[4] = {0.25, 0.75, 0.25, 0.75};
    double zCell[4] = {0.25, 0.25, 0.75, 0.75};
    //Next, the atoms are placed on the fcc lattice. If N 6= 4M3 then some of the lattice sites are left unoccupied.
    int n = 0; // atoms placed so far
    for (int x = 0; x < M; x++)
        for (int y = 0; y < M; y++)
            for (int z = 0; z < M; z++)
                for (int k = 0; k < 4; k++)
                    if (n < N)
                    {
                        r[n][0] = (x + xCell[k]) * a;
                        r[n][1] = (y + yCell[k]) * a;
                        r[n][2] = (z + zCell[k]) * a;
                        ++n;
                    }
}

//Draw initial velocities from a Maxwell-Boltzmann distribution
//The function gasdev from Numerical Recipes returns random numbers with a Gaussian probability distribution
//with center at 0 and unit variance. This function uses the Box-M¨uller algorithm
double gasdev()
{
    static bool available = false;
    static double gset;
    double fac, rsq, v1, v2;
    if (!available)
    {
        do
        {
            v1 = 2.0 * rand() / double(RAND_MAX) - 1.0;
            v2 = 2.0 * rand() / double(RAND_MAX) - 1.0;
            rsq = v1 * v1 + v2 * v2;
        } while (rsq >= 1.0 || rsq == 0.0);
        fac = sqrt(-2.0 * log(rsq) / rsq);
        gset = v1 * fac;
        available = true;
        return v2 * fac;
    }
    else
    {
        available = false;
        return gset;
    }
}

void initVelocities()
{
    // Gaussian with unit variance
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            v[n][i] = gasdev();
    /*
Since these velocities are randomly distributed around zero, the total momentum of the system will be close
to zero but not exactly zero. To prevent the system from drifting in space, the center-of-mass velocity
is computed and used to transform the atom velocities to the center-of-mass frame of reference.
*/

    // Adjust velocities so center-of-mass velocity is zero
    double vCM[3] = {0, 0, 0};
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            vCM[i] += v[n][i];
    for (int i = 0; i < 3; i++)
        vCM[i] /= N;
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            v[n][i] -= vCM[i];
    // Rescale velocities to get the desired instantaneous temperature
    rescaleVelocities();
}
//After setting the CM velocity to zero, the velocities are scaled v->\lambda v
//so that the instantaneous temperature has the desired value T
void rescaleVelocities()
{
    double vSqdSum = 0;
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            vSqdSum += v[n][i] * v[n][i];
    double lambda = sqrt(3 * (N - 1) * T / vSqdSum);
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            v[n][i] *= lambda;
}
/*
Solving Newton’s equations of motion using the same algorithms as in basic_md.cpp with two improvements:
1)periodic boundary conditions will be used to ensure that the number of particles in the simulation
volume remains constant
2)the minimum image convention is used to compute the accelerations:
Since we are using periodic boundary conditions, the system actually has an infinite number of copies of the
N particles contained in the volume L3. Thus there are an infinite number of pairs of particles, all of which
interact with one another! The forces between a particular particle and its periodic copies actually cancel,
but this is not true of pairs which are not images of one another. Since the Lennard Jones interaction is
short ranged, we can safely neglect forces between particles in volumes that are not adjacent to one another.
For adjacent volumes, we have to be more careful. It can happen that the separation is smaller for the next or
previous copy of the system; in that case we use that (smaller) separation.
*/
void computeAccelerations()
{
    for (int i = 0; i < N; i++) // set all accelerations to zero
        for (int k = 0; k < 3; k++)
            a[i][k] = 0;
    for (int i = 0; i < N - 1; i++) // loop over all distinct pairs i,j
        for (int j = i + 1; j < N; j++)
        {
            double rij[3]; // position of i relative to j
            double rSqd = 0;
            for (int k = 0; k < 3; k++)
            {
                rij[k] = r[i][k] - r[j][k];
                // closest image convention
                if (abs(rij[k]) > 0.5 * L)
                {
                    if (rij[k] > 0)
                        rij[k] -= L;
                    else
                        rij[k] += L;
                }
                rSqd += rij[k] * rij[k];
            }
            double f = 24 * (2 * pow(rSqd, -7) - pow(rSqd, -4));
            for (int k = 0; k < 3; k++)
            {
                a[i][k] += rij[k] * f;
                a[j][k] -= rij[k] * f;
            }
        }
}
void velocityVerlet(double dt)
{
    computeAccelerations();
    for (int i = 0; i < N; i++)
        for (int k = 0; k < 3; k++)
        {
            r[i][k] += v[i][k] * dt + 0.5 * a[i][k] * dt * dt;
            //Once the atom is moved, periodic boundary conditions are imposed to move it back into the system volume
            //if it has exited. This done for each component of the position as soon as it has been updated:
            // use periodic boundary conditions
            if (r[i][k] < 0)
                r[i][k] += L;
            if (r[i][k] >= L)
                r[i][k] -= L;
            v[i][k] += 0.5 * a[i][k] * dt;
        }
    computeAccelerations();
    for (int i = 0; i < N; i++)
        for (int k = 0; k < 3; k++)
            v[i][k] += 0.5 * a[i][k] * dt;
}
//The instantaneous temperature is computed as in simple_md.cpp from the equipartition formula
double instantaneousTemperature()
{
    double sum = 0;
    for (int i = 0; i < N; i++)
        for (int k = 0; k < 3; k++)
            sum += v[i][k] * v[i][k];
    return sum / (3 * (N - 1));
}
int main()
{
    initialize();
    double dt = 0.01;
    ofstream file("T2.data");
    for (int i = 0; i < 1000; i++)
    {
        velocityVerlet(dt);
        file << instantaneousTemperature() << "\n";
        if (i % 200 == 0)
            rescaleVelocities();
    }
    file.close();
}
/*The simulation is run for 1000 time steps. After every 200 steps, the velocities of the atoms are rescaled to
drive the average temperature towards the desired value. The output shows that
1)The temperature rises rapidly from the desired value T = 1.0 when the simulation is started. Why?
2)It takes a few rescaling to push the temperature back to the desired value, and then the system appears
to be in equilibrium.
*/