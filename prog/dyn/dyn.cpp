//this program winds the clock forward on the distribution solved for in comb.cpp.

#include <iostream>
#include <fstream>
#include <Eigen/Dense>

using namespace std;

//global constants
const double 	HALP    	= 0.02;   	//home meeting hazard
const double 	FALP    	= 0.015;   	//foreign meeting hazard
const int 	GRIDSIZE 	= 3000; 	//number of grid points
const int 	PERIODS  	= 100;  	//number of periods
const double    PERIOD_LENGTH 	= 0.2;		//length of a period (say a unit is a year)


int main()
{

	//declare all variables
	Eigen::VectorXd h(GRIDSIZE); 
	Eigen::VectorXd f(GRIDSIZE); 
	Eigen::VectorXd hd(GRIDSIZE); 
	Eigen::VectorXd fd(GRIDSIZE); 
	Eigen::VectorXd chd(GRIDSIZE); 
	Eigen::VectorXd cfd(GRIDSIZE); 
	Eigen::MatrixXd fd_dyn(GRIDSIZE,PERIODS);
	Eigen::MatrixXd hd_dyn(GRIDSIZE,PERIODS);
	Eigen::VectorXd interppoint(GRIDSIZE);  //interpolation point vector
	Eigen::MatrixXd interpweights(GRIDSIZE,2);
	Eigen::VectorXd cut_off(PERIODS);
	double newspot; 	
	double propercut;
	
	//read in distribution data
	ifstream h_read;
	ifstream f_read;
	ifstream hd_read;
	ifstream fd_read;
	ifstream chd_read;
	ifstream cfd_read;
	ifstream cut_read;
	
	h_read.open("../comb/h.csv");
	f_read.open("../comb/f.csv");
	hd_read.open("../comb/hd.csv");
	fd_read.open("../comb/fd.csv");
	chd_read.open("../comb/chd.csv");
	cfd_read.open("../comb/cfd.csv");
	cut_read.open("../comb/cut.csv");
	double x;
	int k = 0;
	while(h_read >> x)
	{
	 	h[k] = x;
	 	k++;
	}
	k = 0;
	while(f_read >> x)
	{
	 	f[k] = x;
	 	k++;
	}
	k = 0;
	while(hd_read >> x)
	{
	 	hd[k] = x;
	 	k++;
	}
	k = 0;
	while(fd_read >> x)
	{
	 	fd[k] = x;
		fd_dyn(k,0) = fd[k];
	 	k++;
	}
	k = 0;
	while(chd_read >> x)
	{
	 	chd[k] = x;
	 	k++;
	}
	k = 0;
	while(cfd_read >> x)
	{
	 	cfd[k] = x;
	 	k++;
	}
	cut_off[0]  << cut_read;

	//get dynamics of the foreign distribution
	newspot = exp(HALP * PERIOD_LENGTH);
	for (int k=0;k<GRIDSIZE;k++)
	{	
		interppoint[k] = 0;
		for (int m=0;m<GRIDSIZE;m++)
		{
			if (f[m]>=newspot*f[k])
			{
				interppoint[k] = m;		
				interpweights(k,0) =  (f[m] - newspot*f[k]) / (f[m] - f[m-1]); 
				interpweights(k,1) = 1 - interpweights(k,0);
				break;
			}	
		}
	}
	// do the interpolation
	for (int t=1;t<PERIODS;t++)
	{
		for (int k=0;k<GRIDSIZE;k++)
		{
			if (interppoint[k]!=0)
			{	
				fd_dyn(k,t) = interpweights(k,0) * fd_dyn(interppoint[k]-1,t-1) 
					+ interpweights(k,1) * fd_dyn(interppoint[k],t-1); //scale over
				fd_dyn(k,t) = newspot * fd_dyn(k,t); //scale up
			}
			else
				fd_dyn(k,t) = fd[GRIDSIZE-1]; //things don't change after a while.
		}
	}

	// APPROXIMATE CUT-OFF 
	for (int t=0;t<PERIODS;t++)
	{
		//Initialize value function	
		const double DENOM = RHO - THETA * HALP;
		for (int k=0;k<GRIDSIZE;k++)
		{
			val0(k)	= pow(h[k],THETA)/DENOM;
		}	
		val1 = val0.array() + 1;	

		//START ITERATION  
		diff = val1-val0;	

		while (diff.dot(diff)/val0.dot(val0)>1e-5)
		{
			
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
			mig_cut = GRIDSIZE - 1;
			found = 0;
			for (int k=0;k<GRIDSIZE;k++)
			{
				stayhome[k] = ph[k] + HALP * Sh[k];
				goabroad[k] = pf[k] + FALP * Sf[k];
				if (found == 0 & stayhome[k] < goabroad[k])
				{
					found = 1;
					mig_cut = k;
				}
			}
			cout << "The cut_off is currently at " << h[mig_cut] << endl << endl;	
			// with cut-off in hand, solve for value function using the (really cool!) Lucas Moll matrix method.
			for (int k=0;k<GRIDSIZE-1;k++)
			{
				if (k<mig_cut)
				{
					B(k,k) = (RHO-THETA*HALP+HALP*chd[k]+HALP*h[k]/INCREMENT);
					for (int m=0;m<k+1;m++)
						C(k,m) = hd[m] * HALP * INCREMENT;
					b[k] = ph[k];
					if (k<GRIDSIZE-1)
						B(k,k+1) = -HALP * h[k] / INCREMENT;
				}
				else
				{	
					B(k,k) = (RHO-THETA*HALP+FALP*cfd[k]+HALP*h[k]/INCREMENT);
					for (int m=0;m<k+1;m++)
						C(k,m) = fd[m] * FALP * INCREMENT;
					b[k] = pf[k];
					if (k<GRIDSIZE-1)
						B(k,k+1) = -FALP * f[k] / INCREMENT;
				}
			}

			val0 = val1; //Read in old value
			A = B-C;
			A.llt().solveInPlace(b);
			val1 = b;
			//val1 = A.inverse()*b; //Alternative (direct inversion)
			diff = val1-val0;	
		}
	}
}

