#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include "mersenne.h"

using namespace std;

#define ARRAYELEMCOUNT(x) sizeof(x) / sizeof(x[0])

#define TAB "\t"

//Rand01IE: Get random number from [0,1)
#define Rand01IE() mersenne()
#define SeedRand(x) seed_mersenne(x)

inline double f(const double x) {return 1.5*(1-x*x);}
//---------------------------------------------------


int main()
//--------
{
	const int seed = time(0);
	SeedRand(seed);

	ofstream fout("hw2.txt", ios::trunc);
	if(fout.fail()) {cerr << "Failed to open file." << endl;}
	ostream& oStrm = fout;
	//If you want output on screen, change fout -> cout in above line.
	oStrm.precision(12);
	oStrm << "#seed: " << seed << endl;

	const int intervalcounts[] = {10,100,1000,10000};

	for(int n = 0; n<ARRAYELEMCOUNT(intervalcounts); n++)
	{
		oStrm << "#N = " << intervalcounts[n] << endl;
		double sum = 0;
		const double Iinterval = 1.0;
		const double exactvalue = 1.0;
		const double h = Iinterval/intervalcounts[n];

		//Midpoint
		for(int i = 0; i<intervalcounts[n]; i++)
		{
			sum += f((0.5+i)*h);
		}
		sum *= h;
		oStrm << sum << TAB << fabs(sum-1.0);

		//MC
		sum = 0;
		double sum2 = 0;
		for(int i = 0; i<intervalcounts[n]; i++)
		{
			const double temp = f(Rand01IE());
			sum += temp;
			sum2 += temp*temp;
		}
		const double Fexpect = sum / intervalcounts[n];
		const double F2expect = sum2 / intervalcounts[n];
		const double error = fabs(sum*h-exactvalue);
		const double errorOneSigma = Iinterval * sqrt((F2expect - Fexpect*Fexpect)/ (intervalcounts[n]-1.0));
		oStrm << TAB << Iinterval*Fexpect << TAB << error << TAB << errorOneSigma << endl;
	}

	return 0;
}
