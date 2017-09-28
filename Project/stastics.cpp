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
