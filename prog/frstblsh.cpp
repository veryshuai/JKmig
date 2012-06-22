//first pass at simulation of migration model dynamics, now in C++
#include <iostream>
#include <string>
#include <math.h>
#include <cstdlib>
#include <algorithm>
#include <fstream>

using namespace std;

#define alp     0.02 //meeting hazard
#define B       0.0 //penalty for living abroad
#define theta  -0.5 //production function parameter
#define popsize 10000 //population of countries
#define rho     0.05 //time discount factor
#define cutoff  0.80 //cut-off for migration
#define fparam  15.00 //the higher this is, the better foreign is relative to home

//global parameters
int pt_num = 3000; //number of points in kernel density estimation

//calculates production given cost and location
float prod (float cost, bool home)
{
	if (home==1)
		return (pow(cost,theta));
	if (home == 0)	
		return (pow(cost,theta)/(1.0+B));
	else
		cout << "Error in GDP calculation" << endl;
}

//float init_val (float cost) //copied from Moll-Lucas
//{
//	return(pow(cost,theta)/(rho+theta*alp));
//}

void kernel_density (float hd [], float bw, int max_val, float dens [], float grid [])
{
	float gap    = (float) max_val / (float) pt_num;
	float cur_pt;
	float summation;
	float arg;
	float runsum = 0.0;;

	for (int k=0;k<pt_num;k++)
	{
		cur_pt = gap*(float) k + .001;
		grid [k] = cur_pt;
		summation = 0.0;
		for (int m=0;m<popsize;m++)
		{
			arg = (cur_pt - hd[m])/bw;
			if (arg>0)
				summation += exp(-arg);
		}
		dens[k] = summation/((float) pt_num*bw);
		runsum += dens[k];
	}

	for (int k=0; k<pt_num ; k++)
		dens[k] = dens[k]/runsum;
}

float gdp_calc (float hd [], bool ah [])
{
	float gdp = 0;
	float addme;
	for (int k=0;k<popsize;k++)
	{
		addme = prod(hd[k],ah[k]);
		if (addme<1e32)
			gdp += addme;
		else
			cout << "Error: Infs in GDP calculation." << endl;
	}
	return(gdp);
}

void cum_sum (float dens [], float cdens [])
{
	float summation = 0.0;

	for (int m=0;m<pt_num;m++)
	{
		summation += dens[m];
		cdens[m] = summation;
	}
	
	for (int m=0;m<pt_num;m++)
	{
		cdens[m] = cdens[m]/cdens[pt_num-1];
	}
}

int main () 
{
	//output stuff
	ofstream write2me_h;
	ofstream write2me_f;
	ofstream gdps;
	ofstream growth;
  	write2me_h.open("out_h.csv");
  	write2me_f.open("out_f.csv");
	gdps.open("gdps.csv");
	growth.open("growth.csv");
	
	//Initialize distributions
	float hd [(int) popsize];
	float fd [(int) popsize];
	float hcdf,fcdf;
	for (int k=0; k<popsize+1;k++)
	{
		hcdf  = (float) rand()/(RAND_MAX + 1.0);
		fcdf  = (float) rand()/(RAND_MAX + 1.0);
 		hd[k] = - (float) log(1.0-hcdf);
 		fd[k] = - (float) log(1.0-fcdf)/fparam;
	}
	
	//sort
	sort(hd,hd+popsize);
	sort(fd,fd+popsize);
        	
	//at home?
	bool ah [popsize];
	bool af [popsize];
	int placemarker = 0;
	for (int k=0; k<popsize; k++)
	{
		ah[k] = 0;
		af[k] = 1;
	}
	while (placemarker<popsize*.9)
	{
		ah[placemarker] = 1;
		placemarker++;
	}

	//kernel-density
	int max_val = 5;
	float bw = 2.0;
	float hdens [pt_num];
	float fdens [pt_num];
	float hcum [pt_num];
	float fcum [pt_num];
	float grid [pt_num];
	float co;
	kernel_density(hd,bw,max_val,hdens,grid);
	kernel_density(fd,bw,max_val,fdens,grid);
		
	//write densities to file
	for (int k = 0; k<pt_num ; k++)
	{
		write2me_h << hdens[k] << " ";
	}
	for (int k = 0; k<pt_num ; k++)
	{
		write2me_f << fdens[k] << " ";
	}
	write2me_h << endl;
	write2me_f << endl;

	//write gdp's to file
	float old_gdp_h;
	float old_gdp_f;
	float new_gdp_f;
	float new_gdp_h;
	gdps << "Home and Foreign" << endl;
	old_gdp_h = gdp_calc(hd,ah);
	old_gdp_f = gdp_calc(fd,af);
	gdps << old_gdp_h << " " << old_gdp_f << endl;

	//loop to get density progression
	for (int k = 0; k<10 ; k++)
	{
		cout << "Loop No: " << k+1 << endl;

		cum_sum(fdens,fcum);	
		cum_sum(hdens,hcum);
		int coi = popsize * cutoff; //cut off index
		float new_hd [popsize];
		float new_fd [popsize];
		
		//home at home
		for (int n = 0; n<coi; n++)
		{
			new_hd[n] = hd[n];
			if (-log(1.0-(rand()/(RAND_MAX + 1.0)))/(cutoff*alp)<1.0)
			{
				float meeted = hd[rand() % (coi-1)];
				if (meeted < hd[n])
				{
					new_hd[n] = meeted; 	
				}
			}
		}
			
		//home abroad
		for (int n = coi; n<popsize; n++)
		{	
			new_hd[n] = hd[n];
			
			//meet a foreign guy
			if (-log(1.0-(rand()/(RAND_MAX + 1.0)))/alp<1.0)
			{
				float meeted = fd[rand() % (popsize-1)];
				if (meeted < hd[n])
					new_hd[n] = meeted; 	
			}
			
			//meet a home guy
			if (-log(1.0-(rand()/(RAND_MAX + 1.0)))/((1.0-cutoff)*alp)<1.0)
			{
				float meeted = hd[(rand() % (popsize - coi - 1)) + coi];
				if (meeted < new_hd[n])
					new_hd[n] = meeted; 	
			}
		}

		//foreign
		for (int n = 0; n<popsize; n++)
		{	
			new_fd[n] = fd[n];
			//meet a foreign guy
			if (-log(1.0-(rand()/(RAND_MAX + 1.0)))/alp<1.0)
			{
				float meeted = fd[rand() % (popsize-1)];
				if (meeted < fd[n])
					new_fd[n] = meeted; 	
			}
			
			//meet a home guy
			if (-log(1.0-(rand()/(RAND_MAX + 1.0)))/((1.0-cutoff)*alp)<1.0)
			{
				float meeted = hd[rand() % (popsize - coi - 1) + coi];
				if (meeted < new_fd[n])
					new_fd[n] = meeted; 	
			}
		}

		//read in the new values
		for (int n = 0; n<popsize; n++)
		{
			hd[n] = new_hd[n];
			fd[n] = new_fd[n];
		}
		
		//sort 
		sort(hd,hd+popsize); 
		sort(fd,fd+popsize);
		 
		kernel_density(hd,bw,max_val,hdens,grid);
		kernel_density(fd,bw,max_val,fdens,grid);
		//write densities to file
		for (int m = 0; m<pt_num ; m++)
		{
			write2me_h << hdens[m] << " ";
		}
		for (int m = 0; m<pt_num ; m++)
		{
			write2me_f << fdens[m] << " ";
		}
		write2me_h << endl;
		write2me_f << endl;

		new_gdp_h = gdp_calc(hd,ah);
		new_gdp_f = gdp_calc(fd,af);
		growth << (new_gdp_h-old_gdp_h)/old_gdp_h << " " << (new_gdp_f-old_gdp_f)/old_gdp_f << endl;
		gdps << new_gdp_h << " " << new_gdp_f << endl;
		old_gdp_h = new_gdp_h;
		old_gdp_f = new_gdp_f;
	}

	return(0);
}
