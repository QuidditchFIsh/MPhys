//testing the leapfrog method
#include <stdlib.h>
#include <stdio.h>
#include <random>
#include <iostream>
#include <string>
#include <ctime>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

int main()
{
	double p=1,q=0,q_old=0,t_step =0.1,sum =0,min=0,H_old=0,H_new=0,r=0;

	FILE * out;
	out = fopen("output","w");

	default_random_engine generator(random_device{}());
 	normal_distribution<double> distribution(0.0,1.0);

	for(unsigned int i = 0;i<1000;i++)
	{
		//pick a random p

		p = distribution(generator);
		H_old =(p*p*0.5) + (0.5*q_old*q_old);
		p = p - (0.5*t_step*q_old);

		for(int j=0;j<15;j++)
		{
			q = q + (t_step * p);
			p = p - (t_step * q);
		}

		p = p - (0.5*t_step*q);

		H_new =(p*p*0.5) + (0.5*q*q);

		min = (1 < exp(H_old - H_new)) ? 1 : exp(H_old - H_new);

		r = ((double) rand() / (RAND_MAX));

		if(r<min)
		{
			q_old =q;

		}
		sum+=q;
		fprintf(out,"%f %f %f %d\n",p,q,H_old-H_new,i);
	}
	printf("%f\n",sum/1000);



	

}