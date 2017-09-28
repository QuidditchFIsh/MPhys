/*
	Author: Aneirin John Baker 
	Date: 23/09/17
	Description: This script is where the functions to run the stastics are housed. It will include methods
		which will calculate the mean and other such stastical variables.
*/

#include "stastics.h"

double avgX(vector<vector<vector<double> > > &results)
{
	double sum =0;

	for(unsigned int i=0;i<100;i++)
	{
		sum += results[0][i][1];
	}

	return sum/100;
}

void avg_X_Sqd(vector<vector<vector<double> > > &results)
{	
	double sum =0;

	for(unsigned int i=0;i<100;i++)
	{
		sum += (results[0][i][1]*results[0][i][1]);
		results[0][i][2] = (sum/(double)i);

	}
}