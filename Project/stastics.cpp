/*
	Author: Aneirin John Baker 
	Date: 23/09/17
	Description: This script is where the functions to run the stastics are housed. It will include methods
		which will calculate the mean and other such stastical variables.
*/

#include "stastics.h"

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
	
	FILE * output1;
	output1 = fopen("HMC_Results_x_2","w");

	double sum =0;

	for(unsigned int i=0;i<1000;i++)
	{
		for(unsigned int j=0;j<i;j++)
		{
			sum += results[j][2]*results[j][2];
		}
		fprintf(output1,"%d %f \n ",i,sum/(double) i);
		stats_data[i][1]=sum/(double) i;
		sum=0;

	}

	//return sum/length;
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

	for(unsigned int i=0;i<500;i++)
	{
		//fprintf(output2,"%f\n",ACT[i]);
		fprintf(output2,"%d %f\n",i,ACT[i]);
	}

}


