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

void moving_avg_X_Sqd(vector<vector<double> > results,vector<vector<double> > &stats_data,unsigned int iterations,unsigned int length)
{	
#if 0
	FILE * output1;
	output1 = fopen("HMC_Results_x_2_10","w");

	double sum =0;
	printf("%d %d\n",results.size(),results[0].size());
	printf("%d\n",iterations/5);
	for(unsigned int i=1;i<1000;i++)
	{
		for(unsigned int j=1;j < i;j++)
		{	
			sum += results[j][500]*results[j][500];
			//printf("%d %d\n",i,j);

		}
		fprintf(output1,"%d %f \n ",i,sum/(double) i);
		//printf("%f \n",sum/(double) i);
		stats_data[i][1]=sum/(double) i;
		sum=0;

	}

	//return sum/length;
}
#endif

#if 1
	
	FILE * output1;
	output1 = fopen("HMC_Results_x_2","w");

	double sum =0,sum1=0;
	printf("%d %d\n",results.size(),results[0].size());
	printf("%d\n",iterations/5);
	for(unsigned int i=1;i<1000;i++)
	{
		for(unsigned int j=1;j < i;j++)
		{	
			for(unsigned int k=0;k<length;k++)
			{
				sum1 +=results[j][k]*results[j][k];
			}
			sum1 = sum1/length;
			sum += sum1;
			sum1=0;
			//printf("%d %d\n",i,j);
		}
		fprintf(output1,"%d %f \n ",i,sum/(double) i);
		//printf("%f \n",sum/(double) i);
		stats_data[i][1]=sum/(double) i;
		sum=0;

	}

	//return sum/length;
}
#endif
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
	//avgxx = avg_X_Sqd(sample);

	return standard_Deviation(avgx,avgxx,(double)len);

}


void autocorrelation_Time(vector<vector<double> > data,unsigned int iterations,unsigned int length)
{

	FILE * output2;
	output2 = fopen("HMC_Results_ACT","w");

	vector<double> stats(iterations/5,0);

	vector<double> ACT(iterations/5,0);

	double sum1=0,sum2=0;


	for(unsigned int i=1;i<iterations/5;i++)
	{
		stats[i-1] = data[i][1];
		//printf("%f\n",stats[i]);
	}

	double avgx2 = avgX(stats),x2_0 = stats[0];

	for(unsigned int i=0;i<iterations/5;i++)
	{
		sum2 += (stats[i]-avgx2)*(stats[i]-avgx2);
	}
	//printf("%f %f %f\n",sum2,x2_0,avgx2);
	for(unsigned int i=0;i<iterations/5;i++)
	{
		for(unsigned int j=i+1;j<iterations/5-i;j++)
		{
			sum1 += (stats[j]-avgx2)*(stats[j+i]-avgx2);
		}
		ACT[i]=sum1/sum2;
		sum1=0;
	}

	for(unsigned int i=0;i<iterations/5;i++)
	{
		//fprintf(output2,"%f\n",ACT[i]);
		fprintf(output2,"%d %f\n",i,ACT[i]);
	}

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
	return (p*p*0.5) + (pow((q_plus - q),2)*0.5) + (-4 * pow(q,2)) + (0.1*pow(q,4));
	//return (p*p*0.5) + (pow((q_plus - q),2)*0.5) + (0.25*pow((q*q) - 1,2));
}

double Anarmonic_action(double q, double q_plus)
{
	return (0.5*pow((q_plus - q),2) + (-4 * pow(q,2)) + (0.1*pow(q,4)));
	//return (0.5*pow((q_plus - q),2) + (0.25*pow((q*q) - 1,2)));
}
double kinetic_Energy(double p)
{
	return (p * p * 0.5);
}

double Energy_gap_calc(vector<vector<double> > energyArray,int length)
{
	double sum1 =0,sum2=0;

	for(unsigned int i=0;i<length;i++)
	{
		sum1 += energyArray[0][i] * energyArray[1][i];
		sum2 += energyArray[0][i] * energyArray[2][i];
	}
	sum1 = sum1/length;
	sum2 = sum2/length;

	printf("%f %f\n",sum1,sum2);

	return log(sum1/sum2);

}



