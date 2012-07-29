//Distribution solver for migration model
//Take cut-off as given,
//Calculate distribution.

#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/LU>

using namespace std;

const double 	HALP    = 0.02;   //home meeting hazard
const double 	FALP    = 0.018;  //foreign meeting hazard
const double    FLAM    = 1;      //foreign distribution lambda
const double	HLAM 	= 1;      //home distribution lambda
const int 	GRIDSIZE= 100000; //number of grid points
const int       GRIDMAX = 100;    //maximum grid value

const double    INCREMENT = (double) GRIDMAX/ (double) GRIDSIZE; //step size for the grid
const int 	MIG_CUT   = 90000; //cut-off for migratation, this will be exogenous

int main () 
{
	//Initialize grids
	Eigen::VectorXd h(GRIDSIZE); //grid values
	Eigen::VectorXd f(GRIDSIZE);
	Eigen::VectorXd fd(GRIDSIZE);
	for (int k=0;k<GRIDSIZE;k++)
	{
 		h[k]  =  INCREMENT * (k+1);
		f[k]  =  INCREMENT * (k+1); 
		fd[k] =   FLAM / pow(1 + FLAM * f[k],2);
	}

	//BEGIN DISTRIBUTION CALCULATION  

	//first declare all necessary variables	
	Eigen::VectorXd hd(GRIDSIZE); //grid densities
	Eigen::VectorXd chd(GRIDSIZE); //cumulative distribution
	Eigen::VectorXd cdiff(GRIDSIZE); //cumulative distribution finite difference
	
	//some constants
	const double BAREA  = FALP / HALP * FLAM * h[MIG_CUT] / (1 + FLAM * h[MIG_CUT]); 
	const double LSLOPE = 1 * FALP / HALP * FLAM / pow(1 + FLAM * h[MIG_CUT],2);	
	const double RSLOPE = 1 * FALP / HALP * FLAM / pow(1 + FLAM * h[MIG_CUT],2);	
	const double AAREA  = 1 - BAREA;
	
	//loop to make the area under the migration cut-off consistent
	double barea_diff = 1;	
	double barea = BAREA;
	double barea_old = BAREA;
	while  (barea_diff>1e-3)
	{
		//since the first point is zero, we need to make an assumption about the first derivative, which will be based on the density HLAM/(1+HLAM x)^2
		hd[0] = HLAM;
		//hd[MIG_CUT] = LSLOPE;
		//chd[MIG_CUT] = BAREA;
		
		//now the main loop to calculate the distribution
		for (int k=0;k<MIG_CUT;k++)
		{
			hd[k+1] = hd[k] + INCREMENT / (FLAM * h[k]) * (-hd[k] - hd[k] * chd[k] + hd[k] * (barea - chd[k]) + fd[k] * (1-barea) * FALP / HALP);
			hd[k+1] = max(hd[k+1],0.0);
			chd[k+1] = chd[k] + hd[k+1] * INCREMENT;
		}
		//going backwards
		// for (int k=MIG_CUT;k>0;k--)
		// {
		// 	hd[k-1] = hd[k] - INCREMENT / (FLAM * h[k]) * (-hd[k] - hd[k] * chd[k] + hd[k] * (BAREA - chd[k]) + fd[k] * AAREA * FALP / HALP);
		// 	hd[k-1] = max(hd[k-1],0.0);
		// 	chd[k-1] = chd[k] - hd[k-1] * INCREMENT;
		// }

		for (int k=MIG_CUT;k<GRIDSIZE;k++)
		{
			hd[k] = hd[k-1] + INCREMENT / (FLAM * h[k-1]) * (-hd[k-1] - hd[k-1] * FALP / HALP / (1 - FLAM * h[k-1]) + fd[k-1] * (1 - chd[k-1]) * FALP / HALP); 
			hd[k] = max(hd[k],0.0);
			chd[k] = chd[k-1] + hd[k] * INCREMENT;
		}
		
		//output
		cout << endl;
		for (int k=0;k<15;k++)
		{
			cout << h[k] << ' ' << hd[k] << ' ' << chd[k] << endl;
		}
		cout << chd[GRIDSIZE-1] << endl;
		cout << hd[GRIDSIZE-1] << endl; 

		barea_old = barea;
		barea = chd[MIG_CUT];
		barea_diff = abs(barea - barea_old);
		cout << "Difference is " << barea_diff << endl;
	}	

	//rescale the migrant section of the distribution to get correct mass.
	double scale = chd[GRIDSIZE-1];
	for (int k=MIG_CUT;k<GRIDSIZE;k++)
	{
		hd[k] = hd[k] / scale;
		chd[k] = chd[k] / scale;
	}
	cout << "Scaled cumulative dist ends at " << chd[GRIDSIZE-1] << endl;

	return(0);
}


