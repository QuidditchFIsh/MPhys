	vector<double> v1(3,0);
	vector<vector<double> >States(length,v1);

	double H_Old=0,H_New=0;
	int acceptance=0;
	//Random Number generator.
	//pick out new random P
	default_random_engine generator;
 	normal_distribution<double> distribution(0.0,1.0);

 	for(unsigned int i = 1; i < iterations; i++)
 	{
 		for(unsigned int j = 0; j < length;j++)
 		{
 			//initalise the p's for the state. The first q's will be 0 and evolved randomly by the HMC algorithm.
 			saved_State[1][j] = distribution(generator);
 			saved_State[0][j] = lattice[i][j];
 		}	

 		H_Old = hamiltonian(saved_State[1],saved_State[0],length,t_step);
 		for(unsigned int i = 0; i < length ; i++)
 		{
 			lattice[i][j]= hmcAlgorithm(t_step,distribution(generator),lattice[i][j-1]);	
 		}

 		H_New = hamiltonian();

 		//Now accept or reject the new state with Metroplis update 
 		if(exp(H_Old - H_New) < 1)
 		{
 			//accept the state 
 			acceptance++;
 		}
 		else
 		{
 			//reject the state and reset the state back to the old one
 			for(unsigned int i = 0; i < length; i++)
 			{
 				//set the lattice to the old state
 				lattice[j][i] = saved_State[i][0];
 			}
 		}
	}	