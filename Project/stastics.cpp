/*
	Author: Aneirin John Baker 
	Date: 23/09/17
	Description: This script is where the functions to run the stastics are housed. It will include methods
		which will calculate the mean and other such stastical variables.
*/

#include "stastics.h"
#define Stats_Flip 0
//1 = Harmonic
//0 = Anharmonic


double avgX(vector<double> results)
{	
	double sum =0;
	int length = results.size();

	for(int i=0;i<length;i++)
	{
		sum += results[i];
	}

	return sum/length;
}

double avg_X_Sqd(vector<double> results)
{	

	double sum =0;
	unsigned int length = results.size();

		for(unsigned int j=0;j<length;j++)
		{
			sum += results[j]*results[j];
		}

	return sum/(double)length;
}

double avg_X_four(vector<double> results)
{

	double sum =0;

	unsigned int length = results.size();

		for(unsigned int j=0;j<length;j++)
		{
			sum += pow(results[j],4);
		}

	return sum/(double)length;
}


double standard_Deviation(double avg_X_Sqd, double avgX,double length )
{
	double std_dev = (avg_X_Sqd - pow(avgX,2)) / (length -1.0);

	return sqrt(std_dev);
}


double error_Bars(vector<double> results)
{

	//using bootstrap algorithm to calcuate the error on the bars
	//FIND A BETTER WAY TO DO THIS

	double length =results.size() * 0.8,avgx,avgxx;
	int len = (int)length,rand_No;

	vector<double> sample(1,0);
	//THIS COULD CAUSE SOME TROUBLE IN THE STATS

	for(int i =0; i< len ; i++)
	{
		rand_No = rand() % results.size();
		//sample.insert(results[rand_No]);
	}

	avgx = avgX(sample);
	avgxx = avg_X_Sqd(sample);

	return standard_Deviation(avgx,avgxx,(double)len);

}

double lattice_Hamiltonian(vector<vector<double> > state,unsigned int length )
{
	double H=0;
	//loop for all sites which are not effected by periodic BC's
#if Stats_Flip
	for(unsigned int i=0;i<length-1;i++)
	{
		H += Harmonic_hamiltonian(state[0][i],state[1][i],state[1][i+1]);
	}
	//Periodic BC sites
	H += Harmonic_hamiltonian(state[0][length-1],state[1][length-1],state[1][0]);
#endif

#if !Stats_Flip
	for(unsigned int i=0;i<length-1;i++)
	{
		H += Anarmonic_hamiltonian(state[0][i],state[1][i],state[1][i+1]);
	}
	//Periodic BC sites
	H += Anarmonic_hamiltonian(state[0][length-1],state[1][length-1],state[1][0]);
#endif

	return H;

}
double lattice_Action(vector<double> q,unsigned int length)
{
	double S = 0;
#if Stats_Flip
	for(unsigned int i =0; i<length-1;i++)
	{
		S += Harmonic_action(q[i],q[i+1]);
	}

	S += Harmonic_action(q[length-1],q[0]);
#endif

#if !Stats_Flip
	for(unsigned int i =0; i<length-1;i++)
	{
		S += Anarmonic_action(q[i],q[i+1]);
	}

	S += Anarmonic_action(q[length-1],q[0]);
#endif

	return S/(double)length;
}

double lattice_KineticEnergy(vector<double> p,unsigned int length)
{
	double KE =0;

	for(unsigned int i =0; i<length;i++)
	{
		KE += kinetic_Energy(p[i]);
	}

	return KE/(double)length;
}


double Harmonic_hamiltonian(double p,double q,double q_plus )
{
	return (p*p*0.5) + (pow((q_plus - q),2)*0.5) + (0.5*q*q);
}

double Harmonic_action(double q, double q_plus)
{
	return (0.5*pow((q_plus - q),2) + (0.5 * pow(q,2)));
}
double Anarmonic_hamiltonian(double p,double q,double q_plus )
{
	//return (p*p*0.5) + (pow((q_plus - q),2)*0.5) + (-20 * pow(q,2)) + (0.1*pow(q,4));
	return (p*p*0.5) + (pow((q_plus - q),2)*0.5) + (pow((q*q) - 2,2));
}

double Anarmonic_action(double q, double q_plus)
{
	//return (0.5*pow((q_plus - q),2) + (-20 * pow(q,2)) + (0.1*pow(q,4)));
	return (0.5*pow((q_plus - q),2) + (pow((q*q) - 2,2)));
}
double kinetic_Energy(double p)
{
	return (p * p * 0.5);
}