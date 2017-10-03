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

double avg_X_Sqd(vector<double> results)
{	
	double sum =0;
	int length = results.size();

	for(int i=0;i<length;i++)
	{
		sum += results[i]*results[i];
	}

	return sum/length;
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

	vector<double> sample(len,0);
	//THIS COULD CAUSE SOME TROUBLE IN THE STATS

	for(int i =0; i< (int)length ; i++)
	{
		rand_No = rand() % (int) results.size();
		sample[i] = results[rand_No];
	}

	avgx = avgX(sample);
	avgxx = avg_X_Sqd(sample);

	return standard_Deviation(avgx,avgxx,(double)len);

}

double autocorrelation_Time(vector<double> data,double avgx ,int k)
{
	double sum1=0,sum2=0;
	unsigned int length = data.size();
	for(unsigned int i = 0;i < length-k;i++)
	{
		sum1 += (data[i]-avgx) * (data[i+k] - avgx);

	}
	for(unsigned int i=0;i<length;i++)
	{
		sum2 += pow((data[i] - avgx),2);
	}

	return sum1/sum2;
}


