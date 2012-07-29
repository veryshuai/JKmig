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
const double 	B       = 20.0;   //penalty for living abroad
const double    THETA   =-0.5;    //production function parameter
const int 	GRIDSIZE= 2000;  //number of grid points
const int       GRIDMAX = 300;    //maximum grid value
const double 	RHO     = 0.05;   //time discount factor
const double    FLAM    = 1;      //foreign distribution lambda
const double	HLAM 	= 0.3;      //home distribution lambda

// const int 	MIG_INIT   = 0.01 * (float) GRIDSIZE; 		 //initial migration cut-off guess
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
	//Initialize ALL required variables
	Eigen::VectorXd h(GRIDSIZE); //grid values
	Eigen::VectorXd f(GRIDSIZE);
	Eigen::VectorXd fd(GRIDSIZE);
	Eigen::VectorXd hd(GRIDSIZE); //grid densities
	Eigen::VectorXd hd_old(GRIDSIZE); 
	Eigen::VectorXd chd(GRIDSIZE); //cumulative distribution
	Eigen::VectorXd cfd(GRIDSIZE); //cumulative distribution
	Eigen::VectorXd val0(GRIDSIZE); //initial value function
	Eigen::VectorXd val1(GRIDSIZE); //initial value function
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
		//
		//loop to make the area under the migration cut-off consistent
		while  (barea_diff>1e-3)
		{
			hd[0] = HLAM;
			chd[0] = 0;
			
			//main loop to calculate the distribution
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
			
			// //output
			// cout << endl;
			// for (int k=0;k<15;k++)
			// {
			// 	cout << h[k] << ' ' << hd[k] << ' ' << chd[k] << endl;
			// }
			// cout << chd[GRIDSIZE-1] << endl;
			// cout << hd[GRIDSIZE-1] << endl; 

			barea_old = barea;
			barea = chd[mig_cut];
			barea_diff = abs(barea - barea_old);
			cout << "Difference is " << barea_diff << endl;
		}	

		//rescale the migrant section of the distribution to get correct mass.
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
			cout << mig_cut << endl;	
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
			//A.lu().solve(b, &val1); //This is the solve step	
			A.llt().solveInPlace(b);
			val1 = b;
			//val1 = A.inverse()*b; //Alternative (direct inversion)
			diff = val1-val0;	
		}
		dist_diff = hd - hd_old;
	}
	for (int k=0;k<15;k++)
		cout << f[k] << "  " << fd[k] << "  " << hd[k] << "  " << chd[k] << endl;

	return(0);
}
	
