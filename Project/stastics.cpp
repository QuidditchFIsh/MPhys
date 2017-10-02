/*
	Author: Aneirin John Baker 
	Date: 23/09/17
	Description: This script is where the functions to run the stastics are housed. It will include methods
		which will calculate the mean and other such stastical variables.
*/

#include "stastics.h"

double avgX(vector<vector<vector<double> > > &results,unsigned int iterations,unsigned int length)
{	
	double sum =0;

	for(unsigned int i=0;i<length;i++)
	{
		for(unsigned int j=0;j<iterations;j++)
		{
		sum += (results[j][i][0]);
		}
		results[0][i][1] = (sum/(double)iterations);
		sum =0;

	}
}

void avg_X_Sqd(vector<vector<vector<double> > > &results,unsigned int iterations,unsigned int length)
{	
	double sum =0;

	for(unsigned int i=0;i<length;i++)
	{
		for(unsigned int j=0;j<iterations;j++)
		{
		sum += (results[j][i][0]*results[j][i][0]);
		results[j][i][1] = (sum/(double)i);
		}
		sum =0;

	}
#if 1
	for(unsigned int i=0;i<length;i++)
	{
		for(unsigned int j=0;j<iterations;j++)
		{
		sum += (results[j][i][0]);
		}
		results[0][i][1] = (sum/(double)iterations);
		sum =0;

	}
#endif
}
double standard_Deviation(double avg_X_Sqd, double avgX,double length )
{
	double std_dev = (avg_X_Sqd - pow(avgX,2)) / (length -1.0)

	return sqrt(std_dev);
}

void error_Bars(vector<vector<vector<double> > > &results,int iterations)
{

	//using bootstrap algorithm to calcuate the error on the bars
	double length = floor(0.8 * (double)iterations);

	vector<double> sample(length,0);
	




}

void autocorrelation_Time()
{

}


