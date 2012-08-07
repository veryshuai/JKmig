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
#include <Eigen/LU>

using namespace std;

const double 	ALP     = 0.02;   //home meeting hazard
const double 	FALP    = 0.015;   //foreign meeting hazard
const double 	B       = 0.25;   //penalty for living abroad
const double    THETA   =-0.5;    //production function parameter
const int 	GRIDSIZE= 1000;  //number of grid points
const int       GRIDMAX = 100;    //maximum grid value
const double 	RHO     = 0.05;   //time discount factor
const double    FPARAM  = 0.02;   //foreign distribution exponential parameter
const double    HPARAM  = 0.01;   //home distribution exponential parameter

const double    INCREMENT = (double) GRIDMAX/ (double) GRIDSIZE; //step size for the grid

//calculates production given cost and location
double prod (double cost, bool home)
{
	if (home == 1)
		return (pow(cost,THETA));
	if (home == 0)	
		return (pow(cost,THETA)/(1.0+B));
	else
		cout << "Error in production calculation" << endl;
}
// 
// void kernel_density (double hd [], double bw, int max_val, double dens [], double grid [])
// {
// 	double gap    = (double) max_val / (double) pt_num;
// 	double cur_pt;
// 	double summation;
// 	double arg;
// 	double runsum = 0.0;;
// 
// 	for (int k=0;k<pt_num;k++)
// 	{
// 		cur_pt = gap*(double) k + .001;
// 		grid [k] = cur_pt;
// 		summation = 0.0;
// 		for (int m=0;m<POPSIZE;m++)
// 		{
// 			arg = (cur_pt - hd[m])/bw;
// 			if (arg>0)
// 				summation += exp(-arg);
// 		}
// 		dens[k] = summation/((double) pt_num*bw);
// 		runsum += dens[k];
// 	}
// 
// 	for (int k=0;k<pt_num;k++)
// 		dens[k] = dens[k]/runsum;
// }
// 
// double gdp_calc (double hd [], bool ah [])
// {
// 	double gdp = 0;
// 	double addme;
// 	for (int k=0;k<POPSIZE;k++)
// 	{
// 		addme = prod(hd[k],ah[k]);
// 		if (addme<1e32)
// 			gdp += addme;
// 		else
// 			cout << "Error: Infs in GDP calculation." << endl;
// 	}
// 	return(gdp);
// }

void cumsum (Eigen::VectorXd& input, Eigen::VectorXd& output)
{
	double summation = 0.0;

	for (int m=0;m<GRIDSIZE;m++)
	{
		summation += input[m] * INCREMENT;
		output[m] = summation;
	}
}

int main () 
{
	// //output stuff
	// ofstream write2me_h;
	// ofstream write2me_f;
	// ofstream gdps;
	// ofstream growth;
  	// write2me_h.open("out_h.csv");
  	// write2me_f.open("out_f.csv");
	// gdps.open("gdps.csv");
	// growth.open("growth.csv");
	
	//Initialize grids
	Eigen::VectorXd h(GRIDSIZE); //grid values
	Eigen::VectorXd f(GRIDSIZE);
	Eigen::VectorXd hd(GRIDSIZE); //grid densities
	Eigen::VectorXd fd(GRIDSIZE);
	for (int k=0;k<GRIDSIZE;k++)
	{
 		h[k]  =  INCREMENT * (k+1);
		f[k]  =  INCREMENT * (k+1); 
		hd[k] =  HPARAM * (double) exp(-HPARAM * h[k]);
		fd[k] =  FPARAM * (double) exp(-FPARAM * f[k]);
	}

	//Cumulative distributions
	Eigen::VectorXd chd(GRIDSIZE); 
	Eigen::VectorXd cfd(GRIDSIZE); 
	cumsum(hd,chd);	
	cumsum(fd,cfd);	

	//Initialize value function	
	Eigen::VectorXd val0(GRIDSIZE); //initial value function
	Eigen::VectorXd val1(GRIDSIZE); //initial value function
	const double DENOM = RHO - THETA * ALP;
        for (int k=0;k<GRIDSIZE;k++)
	{
		val0(k)	= pow(h[k],THETA)/DENOM;
	}	
	val1 = val0.array() + 1;	

	//START ITERATION  
	
	//first declare all loop variables	
	Eigen::VectorXd dval0(GRIDSIZE); 
	Eigen::VectorXd fth(GRIDSIZE); //initialize first term placeholders
	Eigen::VectorXd ftf(GRIDSIZE); 
	Eigen::VectorXd cfth(GRIDSIZE); 
	Eigen::VectorXd cftf(GRIDSIZE); 
	Eigen::VectorXd Sh(GRIDSIZE); //home integral approx	
	Eigen::VectorXd Sf(GRIDSIZE); //foreign integral approx
	Eigen::VectorXd ph(GRIDSIZE);
	Eigen::VectorXd pf(GRIDSIZE);
	Eigen::MatrixXd A(GRIDSIZE,GRIDSIZE);
       	Eigen::MatrixXd B(GRIDSIZE,GRIDSIZE);	
       	Eigen::MatrixXd C(GRIDSIZE,GRIDSIZE);	
	Eigen::VectorXd b(GRIDSIZE);
	Eigen::VectorXd stayhome(GRIDSIZE);
	Eigen::VectorXd goabroad(GRIDSIZE);
	Eigen::VectorXd diff(GRIDSIZE);
	const bool home = 1;	
	const bool foreign = 0;
	int mig_cut = GRIDSIZE;
	bool found = 0;
	
	diff = val1-val0;	

	while (diff.dot(diff)/val0.dot(val0)>1e-5)
	{
		//numerical derivative
		dval0[0] = 0;
		for (int k=1;k<GRIDSIZE;k++)
			dval0[k] = (val0[k] - val0[k-1])/INCREMENT;
		
		//cumulative first terms
		fth[0] = val0[0] * hd[0]; 
		ftf[0] = val0[0] * fd[0]; 
		for (int k=1;k<GRIDSIZE;k++)	
		{
			fth[k] = fth[k-1] + val0[k] * hd[k];
			ftf[k] = ftf[k-1] + val0[k] * fd[k];
		}			
		cumsum(fth,cfth);
		cumsum(ftf,cftf);
		
		//get integral approxs
		Sh[0] = Sf[0] = 0.0; 
		for (int k=1;k<GRIDSIZE;k++) 
		{
			Sh[k] = cfth[k] - val0[k] * chd[k];
			Sf[k] = cftf[k] - val0[k] * cfd[k];
		}			
		
		//get productions
		for (int k=0;k<GRIDSIZE;k++) 
		{
			ph[k] = prod(h[k],home);
			pf[k] = prod(f[k],foreign);
		}

		//find cutoff for migration
		found = 0;
		for (int k=0;k<GRIDSIZE;k++)
		{
			stayhome[k] = ph[k] + ALP * Sh[k];
			goabroad[k] = pf[k] + FALP * Sf[k];
			if (found == 0 & stayhome[k] < goabroad[k])
			{
				found = 1;
				mig_cut = k;
			}
		}
		cout << mig_cut << endl << endl;	
		// with cut-off in hand, solve for value function using the (really cool!) Lucas Moll matrix method.
		
		for (int k=0;k<GRIDSIZE;k++)
		{
			if (k<mig_cut)
			{
				B(k,k) = (RHO-THETA*ALP+ALP*chd[k]+ALP*h[k]/INCREMENT);
				for (int m=0;m<k+1;m++)
					C(k,m) = hd[m] * ALP * INCREMENT;
				b[k] = ph[k];
				if (k<GRIDSIZE-1)
					B(k,k+1) = -ALP * h[k] / INCREMENT;
			}
			else
			{	
				B(k,k) = (RHO-THETA*ALP+FALP*cfd[k]+ALP*h[k]/INCREMENT);
				for (int m=0;m<k+1;m++)
					C(k,m) = fd[m] * FALP * INCREMENT;
				b[k] = pf[k];
				if (k<GRIDSIZE-1)
					B(k,k+1) = -FALP * f[k] / INCREMENT;
			}
		}
		A = B-C;
		//A.lu().solve(b, &val1); //This is the solve step	
		val0 = val1; //Read in old value
		val1 = A.inverse()*b; //Alternative (direct inversion)
		diff = val1-val0;	
			
		// // Display contents of matrices
		// for (int k=0;k<5;k++)
		// 	cout << B(k,0) << ' ' << B(k,1) << ' ' << B(k,2) << ' ' << B(k,3) << ' ' << B(k,4) << endl;
		// 	
		// for (int k=0;k<5;k++)
		// 	cout << C(k,0) << ' ' << C(k,1) << ' ' << C(k,2) << ' ' << C(k,3) << ' ' << C(k,4) << endl;

		// for (int k=0;k<5;k++)
		// 	cout << A(k,0) << ' ' << A(k,1) << ' ' << A(k,2) << ' ' << A(k,3) << ' ' << A(k,4) << endl;
		// for (int k=0;k<5;k++)
		// 	cout << b(k) << endl;
		for (int k=0;k<5;k++)
		cout << val1(k) << ' ' << val0(k) << endl << endl;
		
	}


	
	
	// //at home?
	// int cut = -1;
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


