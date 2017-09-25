#include <stdint.h>
#include <stdio.h>
#include <string>
/*
Testing the different algorithms which can be used for simulating HMC dyunamics using a basic forward
difference differentiation method 
*/

int main()
{
	printf("Initalising System to q = 0 and p = 1 and the time step is 0.3\n");

	double p =1,q=0,t=0,t_step = 0.3,p_new=0,q_new=0,p_Half_New =0;

	FILE * fp;
	FILE * fp1;
	fp=fopen("results","w");
	fp1= fopen("results1","w");
	fprintf(fp, "%f %f\n",p,q);

	for(int i=0;i<400;i++)
	{
		p_new = p + t_step * q;
		q_new = q - t_step * p_new;
		//printf("%f %f\n ",p_new,q_new);
		fprintf(fp, "%f %f\n",p_new,q_new );

		p = p_new;
		q = q_new;
	}	
	
	p =1;
	q=0;
	t=0;
	t_step = 1.3;
	p_new=0;
	q_new=0;
	p_Half_New =0;

	for(int j=0;j<20;j++)
	{
		p_Half_New = p - (0.5 * t_step * q);

		q_new = q + (t_step * p_Half_New);

		p_new = p_Half_New - (0.5*t_step*q_new);

		fprintf(fp1, "%f %f\n",p_new,q_new );

		p = p_new;
		q = q_new;
		p_Half_New=0;


	}
	
return 0;
}