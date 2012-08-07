// This program solves for a balanced growth path in the Lucas style migration model
//

#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/Cholesky>

using namespace std;

const double 	HALP    = 0.02;   //home meeting hazard
const double 	FALP    = 0.015;   //foreign meeting hazard
const double 	B       = 1.0;   //penalty for living abroad
const double    THETA   =-0.5;    //production function parameter
const int 	GRIDSIZE= 3000;  //number of grid points
const int       GRIDMAX = 300;    //maximum grid value
const double 	RHO     = 0.05;   //time discount factor
const double    FLAM    = 1;      //foreign distribution lambda
const double	HLAM 	= 0.30;      //home distribution lambda

//initial migration cut-off guess
const int  	MIG_INIT   = 100;
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
	//output stuff
	ofstream h_stream; 
	ofstream f_stream; 
	ofstream hd_stream;
	ofstream fd_stream; 
	ofstream chd_stream;
	ofstream cfd_stream; 
	ofstream mig_cut_stream;
  	h_stream.open("h.csv");
  	f_stream.open("f.csv");
  	hd_stream.open("hd.csv");
  	fd_stream.open("fd.csv");
  	chd_stream.open("chd.csv");
  	cfd_stream.open("cfd.csv");
	mig_cut_stream.open("cut.csv");
	
	//Initialize ALL required variables
	Eigen::VectorXd h(GRIDSIZE); //grid values
	Eigen::VectorXd f(GRIDSIZE);
	Eigen::VectorXd fd(GRIDSIZE); //grid densities
	Eigen::VectorXd hd(GRIDSIZE); 
	Eigen::VectorXd hd_old(GRIDSIZE); 
	Eigen::VectorXd chd(GRIDSIZE); //cumulative distributions
	Eigen::VectorXd cfd(GRIDSIZE); 
	Eigen::VectorXd val0(GRIDSIZE); //old value function
	Eigen::VectorXd val1(GRIDSIZE); //new value function
	Eigen::VectorXd fth(GRIDSIZE); //initialize first term placeholders for value function approximation
	Eigen::VectorXd ftf(GRIDSIZE); 
	Eigen::VectorXd cfth(GRIDSIZE); //cumulative first terms for value function approximation
	Eigen::VectorXd cftf(GRIDSIZE); 
	Eigen::VectorXd Sh(GRIDSIZE); //value function integral approx	
	Eigen::VectorXd Sf(GRIDSIZE); 
	Eigen::VectorXd ph(GRIDSIZE);
	Eigen::VectorXd pf(GRIDSIZE);
	Eigen::MatrixXd A(GRIDSIZE,GRIDSIZE);
       	Eigen::MatrixXd B(GRIDSIZE,GRIDSIZE);	
       	Eigen::MatrixXd C(GRIDSIZE,GRIDSIZE);	
	Eigen::VectorXd b(GRIDSIZE);
	Eigen::VectorXd stayhome(GRIDSIZE);
	Eigen::VectorXd goabroad(GRIDSIZE);
	Eigen::VectorXd diff(GRIDSIZE);
	Eigen::VectorXd dist_diff(GRIDSIZE);
	const bool home = 1;	
	const bool foreign = 0;
	int mig_cut = MIG_INIT;
	int mig_cut_old = MIG_INIT+1;
	bool found = 0;
		
	double barea_diff 	= 1; 	
	double barea 		= FALP / HALP * FLAM * h[MIG_INIT] + 1 / (1 + FLAM * h[MIG_INIT]);
	double barea_old 	= barea;
	double scale 		= chd[GRIDSIZE-1];

	for (int k=0;k<GRIDSIZE;k++)
	{
 		h[k]  		=  INCREMENT * (k+1);
		f[k]  		=  INCREMENT * (k+1); 
		fd[k] 		=  FLAM / pow(1 + FLAM * f[k],2); 
		hd[k] 		=  fd[k];
	}
	cumsum(fd,cfd);	

	hd_old = hd.array() + 1;
	dist_diff = hd - hd_old;
	while (dist_diff.dot(dist_diff)/hd_old.dot(hd_old)>1e-5 || mig_cut_old != mig_cut)
	{	
		hd_old = hd;	
		mig_cut_old = mig_cut;
		barea_diff = 1;
		
		//calculate new distribution. Stopping rule is area under mig_cut must be consistent with eq (TBA).
		while  (barea_diff>1e-3)
		{
			hd[0] = HLAM;
			chd[0] = 0;
			
			//loop to calculate the distribution
			for (int k=0;k<mig_cut;k++)
			{
				hd[k+1] = hd[k] + INCREMENT / (FLAM * h[k]) * (-hd[k] - hd[k] * chd[k] + hd[k] * (barea - chd[k]) + fd[k] * (1-barea) * FALP / HALP);
				hd[k+1] = max(hd[k+1],0.0);
				chd[k+1] = chd[k] + hd[k+1] * INCREMENT;
			}

			for (int k=mig_cut;k<GRIDSIZE-1;k++)
			{
				hd[k] = hd[k-1] + INCREMENT / (FLAM * h[k-1]) * (-hd[k-1] - hd[k-1] * FALP / HALP / (1 - FLAM * h[k-1]) + fd[k-1] * (1 - chd[k-1]) * FALP / HALP); 
				hd[k] = max(hd[k],0.0);
				chd[k] = chd[k-1] + hd[k] * INCREMENT;
			}
			
			barea_old = barea;
			barea = chd[mig_cut];
			barea_diff = abs(barea - barea_old);
			cout << "Inconsistency in area under mig_cut is " << barea_diff << endl;
		}	
		cout << endl;

		//rescale the migrant section of the distribution to get correct mass (going to be a discontinuity at hd[mig_cut].
		scale = chd[GRIDSIZE-1];
		if (scale>0)
		{
			for (int k=mig_cut;k<GRIDSIZE-1;k++)
			{
				hd[k] = hd[k] / scale;
				chd[k] = chd[k] / scale;
			}
		}

		//FIND THE NEW CUT-OFF 

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
		dist_diff = hd - hd_old;
	}

      	//write to files
 	for (int m=0;m<GRIDSIZE;m++)
 	{
 		h_stream << h[m] << endl;
 	}
 	for (int m=0;m<GRIDSIZE;m++)
 	{
 		f_stream << f[m] << endl;
 	}
 	for (int m=0;m<GRIDSIZE;m++)
 	{
 		hd_stream << hd[m] << endl;
 	}
 	for (int m=0;m<GRIDSIZE;m++)
 	{
 		fd_stream << fd[m] << endl;
 	}
 	for (int m=0;m<GRIDSIZE;m++)
 	{
 		chd_stream << chd[m] << endl;
 	}
 	for (int m=0;m<GRIDSIZE;m++)
 	{
 		cfd_stream << cfd[m] << endl;
 	}
	mig_cut_stream << mig_cut << endl;

	return(0);
}
	
