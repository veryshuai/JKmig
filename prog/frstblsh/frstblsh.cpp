//I ACCIDENTALLY OVERWROTE SOME OF WHAT WAS HERE...WILL NEED TO BE MODIFIED
//BACK TO RUN.
//Value function solver for migration model
//Take foreign and domestic cost distributions as given,
//Calculate value function and cut-off.

#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <fstream>
#include <Eigen/Dense>

using namespace std::Eigen;

const double 	ALP     = 0.02;   //home meeting hazard
const double 	FALP    = 0.01;   //foreign meeting hazard
const double 	B       = 0.25;   //penalty for living abroad
const double    THETA   =-0.5;    //production function parameter
const int 	GRIDSIZE= 10000;  //number of grid points
const int       GRIDMAX = 250;    //maximum grid value
const double 	RHO     = 0.05;   //time discount factor
const double    FPARAM  = 0.02;   //foreign distribution exponential parameter
const double    HPARAM  = 0.04;   //home distribution exponential parameter

//calculates production given cost and location
double prod (double cost, bool home)
{
	if (home == 1)
		return (pow(cost,THETA));
	if (home == 0)	
		return (pow(cost,THETA)/(1.0+B));
	else
		cout << "Error in GDP calculation" << endl;
}

void kernel_density (double hd [], double bw, int max_val, double dens [], double grid [])
{
	double gap    = (double) max_val / (double) pt_num;
	double cur_pt;
	double summation;
	double arg;
	double runsum = 0.0;;

	for (int k=0;k<pt_num;k++)
	{
		cur_pt = gap*(double) k + .001;
		grid [k] = cur_pt;
		summation = 0.0;
		for (int m=0;m<POPSIZE;m++)
		{
			arg = (cur_pt - hd[m])/bw;
			if (arg>0)
				summation += exp(-arg);
		}
		dens[k] = summation/((double) pt_num*bw);
		runsum += dens[k];
	}

	for (int k=0;k<pt_num;k++)
		dens[k] = dens[k]/runsum;
}

double gdp_calc (double hd [], bool ah [])
{
	double gdp = 0;
	double addme;
	for (int k=0;k<POPSIZE;k++)
	{
		addme = prod(hd[k],ah[k]);
		if (addme<1e32)
			gdp += addme;
		else
			cout << "Error: Infs in GDP calculation." << endl;
	}
	return(gdp);
}

void cum_sum (double dens [], double cdens [])
{
	double summation = 0.0;

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
	
	//Initialize grids
	VectorXd h(GRIDSIZE); //grid values
	VectorXd f(GRIDSIZE);
	VectorXd hd(GRIDSIZE); //grid densities
	VectorXd fd(GRIDSIZE);
	double increment = (double) GRIDMAX/ (double) GRIDSIZE;
	for (int k=0;k<GRIDSIZE;k++)
	{
 		h(k)  =  increment * (k++);
		f(k)  =	 increment * (k++); 
		hd(k) =  HPARAM * (double) exp(-HPARAM * h(k));
		fd(k) =  FPARAM * (double) exp(-FPARAM * f(k));
	}
	
	// //sort
	// sort(hd,hd+POPSIZE);
	// sort(fd,fd+POPSIZE);
        // 	
	// //at home?
	// bool ah [POPSIZE];
	// bool af [POPSIZE];
	// int placemarker = 0;
	// for (int k=0;k<POPSIZE;k++)
	// {
	// 	ah[k] = 0;
	// 	af[k] = 1;
	// }
	// while (placemarker<POPSIZE*.9)
	// {
	// 	ah[placemarker] = 1;
	// 	placemarker++;
	// }

	// //kernel-density
	// int max_val = 5;
	// double bw = 2.0;
	// double hdens [pt_num];
	// double fdens [pt_num];
	// double hcum [pt_num];
	// double fcum [pt_num];
	// double grid [pt_num];
	// double co;
	// kernel_density(hd,bw,max_val,hdens,grid);
	// kernel_density(fd,bw,max_val,fdens,grid);
	// 	
	// //write densities to file
	// for (int k=0;k<pt_num;k++)
	// {
	// 	write2me_h << hdens[k] << " ";
	// }
	// for (int k=0;k<pt_num;k++)
	// {
	// 	write2me_f << fdens[k] << " ";
	// }
	// write2me_h << endl;
	// write2me_f << endl;

	// //write gdp's to file
	// double old_gdp_h;
	// double old_gdp_f;
	// double new_gdp_f;
	// double new_gdp_h;
	// gdps << "Home and Foreign" << endl;
	// old_gdp_h = gdp_calc(hd,ah);
	// old_gdp_f = gdp_calc(fd,af);
	// gdps << old_gdp_h << " " << old_gdp_f << endl;

	// //loop to get density evolution 
	// for (int k=0;k<ITERS;k++)
	// {
	// 	cout << "Loop No: " << k+1 << endl;

	// 	cum_sum(fdens,fcum);	
	// 	cum_sum(hdens,hcum);
	// 	int coi = POPSIZE * CUTOFF; //cut off index
	// 	double new_hd [POPSIZE];
	// 	double new_fd [POPSIZE];
	// 	
	// 	//home at home
	// 	for (int n=0;n<coi;n++)
	// 	{
	// 		new_hd[n] = hd[n];
	// 		if (-log(1.0-(rand()/(RAND_MAX + 1.0)))/(CUTOFF*ALP)<1.0)
	// 		{
	// 			double meeted = hd[rand() % (coi-1)];
	// 			if (meeted<hd[n])
	// 			{
	// 				new_hd[n] = meeted; 	
	// 			}
	// 		}
	// 	}
	// 		
	// 	//home abroad
	// 	for (int n=coi;n<POPSIZE;n++)
	// 	{	
	// 		new_hd[n] = hd[n];
	// 		
	// 		//meet a foreign guy
	// 		if (-log(1.0-(rand()/(RAND_MAX + 1.0)))/FALP<1.0)
	// 		{
	// 			double meeted = fd[rand() % (POPSIZE-1)];
	// 			if (meeted<hd[n])
	// 				new_hd[n] = meeted; 	
	// 		}
	// 		
	// 		// //meet a home guy
	// 		// if (-log(1.0-(rand()/(RAND_MAX + 1.0)))/((1.0-CUTOFF)*ALP)<1.0)
	// 		// {
	// 		// 	double meeted = hd[(rand() % (POPSIZE - coi - 1)) + coi];
	// 		// 	if (meeted<new_hd[n])
	// 		// 		new_hd[n] = meeted; 	
	// 		// }
	// 	}

	// 	//foreign
	// 	for (int n=0;n<POPSIZE;n++)
	// 	{	
	// 		new_fd[n] = fd[n];
	// 		//meet a foreign guy
	// 		if (-log(1.0-(rand()/(RAND_MAX + 1.0)))/ALP<1.0)
	// 		{
	// 			double meeted = fd[rand() % (POPSIZE-1)];
	// 			if (meeted<fd[n])
	// 				new_fd[n] = meeted; 	
	// 		}
	// 		
	// 		// //meet a home guy
	// 		// if (-log(1.0-(rand()/(RAND_MAX + 1.0)))/((1.0-CUTOFF)*ALP)<1.0)
	// 		// {
	// 		// 	double meeted = hd[rand() % (POPSIZE - coi - 1) + coi];
	// 		// 	if (meeted<new_fd[n])
	// 		// 		new_fd[n] = meeted; 	
	// 		// }
	// 	}

	// 	//read in the new values
	// 	for (int n=0;n<POPSIZE;n++)
	// 	{
	// 		hd[n] = new_hd[n];
	// 		fd[n] = new_fd[n];
	// 	}
	// 	
	// 	//sort 
	// 	sort(hd,hd+POPSIZE); 
	// 	sort(fd,fd+POPSIZE);
	// 	 
	// 	//write densities to file
	// 	kernel_density(hd,bw,max_val,hdens,grid);
	// 	kernel_density(fd,bw,max_val,fdens,grid);
	// 	for (int m=0;m<pt_num;m++)
	// 	{
	// 		write2me_h << hdens[m] << " ";
	// 	}
	// 	for (int m=0;m<pt_num;m++)
	// 	{
	// 		write2me_f << fdens[m] << " ";
	// 	}
	// 	write2me_h << endl;
	// 	write2me_f << endl;

	// 	//write gdps and growth rates
	// 	new_gdp_h = gdp_calc(hd,ah);
	// 	new_gdp_f = gdp_calc(fd,af);
	// 	growth << (new_gdp_h-old_gdp_h)/old_gdp_h << " " << (new_gdp_f-old_gdp_f)/old_gdp_f << endl;
	// 	gdps << new_gdp_h << " " << new_gdp_f << endl;
	// 	old_gdp_h = new_gdp_h;
	// 	old_gdp_f = new_gdp_f;
	// }

	return(0);
}
