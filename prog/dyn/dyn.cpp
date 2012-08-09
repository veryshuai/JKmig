//this program winds the clock forward on the distribution solved for in comb.cpp.

#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include "mig.h"

using namespace std;

//global constants
const double 	HALP    	      = 0.02;   //home meeting hazard
const double 	FALP    	      = 0.015;  //foreign meeting hazard
const int 	  GRIDSIZE 	      = 1000;   //number of grid points
const int     GRIDMAX 	      = 300;    //maximum grid value
const int 	  PERIODS  	      = 4;    //number of periods
const double  PERIOD_LENGTH 	= 0.2;	  //length of a period (say a unit is a year)
const double 	B      	        = 1.0;    //penalty for living abroad
const double  THETA   	      =-0.5;    //production function parameter
const double 	RHO     	      = 0.05;   //time discount factor
const double  FLAM   	        = 1;      //foreign distribution lambda
const double	HLAM 		        = 0.15;    //home distribution lambda

const double  INCREMENT = (double) GRIDMAX/ (double) GRIDSIZE; //step size for the grid

struct Distributions {
	Eigen::VectorXd h; 
	Eigen::VectorXd f; 
	Eigen::VectorXd hd; 
	Eigen::VectorXd fd; 
	Eigen::VectorXd chd; 
	Eigen::VectorXd cfd; 
	Eigen::VectorXd cut_off;
	Eigen::MatrixXd fd_dyn;
	Eigen::MatrixXd hd_dyn;
	Eigen::MatrixXd chd_dyn; 
	Eigen::MatrixXd cfd_dyn; 
	Eigen::MatrixXd val_dyn;

	Distributions(int gridSize, int periods) {
    h       =  Eigen::VectorXd(gridSize); 
    f       =  Eigen::VectorXd(gridSize); 
    hd      =  Eigen::VectorXd(gridSize); 
    fd      =  Eigen::VectorXd(gridSize); 
    chd     =  Eigen::VectorXd(gridSize); 
    cfd     =  Eigen::VectorXd(gridSize); 
	  cut_off =  Eigen::VectorXd(periods); 
    fd_dyn  =  Eigen::MatrixXd(gridSize,periods);
    hd_dyn  =  Eigen::MatrixXd(gridSize,periods);
    chd_dyn =  Eigen::MatrixXd(gridSize,periods); 
    cfd_dyn =  Eigen::MatrixXd(gridSize,periods); 
    val_dyn =  Eigen::MatrixXd(gridSize,periods); 
  }
};

void read (Eigen::VectorXd& target, const char* filename, Eigen::MatrixXd& target2) {
  ifstream readme;
  readme.open(filename);
	double x;
	int k = 0;
	while(readme >> x)
	{
	  	target[k] = x;
	  	target2(k,0) = target[k];
	  	k++;
	}
}

void dynfor (Eigen::MatrixXd& fd_dyn, Eigen::MatrixXd& cfd_dyn, Eigen::VectorXd& f){
  Eigen::VectorXd fd(GRIDSIZE);
  Eigen::VectorXd cfd(GRIDSIZE);
  Eigen::VectorXd interppoint(GRIDSIZE);  
  Eigen::MatrixXd interpweights(GRIDSIZE,2);
  double newspot = exp(HALP * PERIOD_LENGTH);
  for (int k=0;k<GRIDSIZE;k++)
  {	
    interppoint[k] = 0;
    for (int m=1;m<GRIDSIZE;m++)
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
        fd[k] = fd_dyn(k,t); //use this to write cfd_dyn
      }
      else
        fd_dyn(k,t) = fd[GRIDSIZE-1]; //things don't change after a while.
    }
    cumsum(fd,cfd);
    for (int k=0;k<GRIDSIZE;k++)
      cfd_dyn(k,t) = cfd[k];
  }
}

int main() {
	Distributions distr(GRIDSIZE,PERIODS);
	
	//read in distribution data
  Eigen::MatrixXd junk_mat(GRIDSIZE,1);
  Eigen::VectorXd junk_vec(GRIDSIZE);
  read(distr.h,"../comb/h.csv",junk_mat);
  read(distr.f,"../comb/f.csv",junk_mat);
  read(distr.hd,"../comb/hd.csv",distr.hd_dyn);
  read(distr.fd,"../comb/fd.csv",distr.fd_dyn);
  read(distr.chd,"../comb/chd.csv",distr.chd_dyn);
  read(distr.cfd,"../comb/cfd.csv",distr.cfd_dyn);
  read(distr.cut_off,"../comb/cut.csv",junk_mat);
  read(junk_vec,"../comb/val.csv",distr.val_dyn);

  // get foreign distribution dynamics
  dynfor(distr.fd_dyn,distr.cfd_dyn,distr.f);
  
  // APPROXIMATE CUT-OFF 

  //first get productions
	Eigen::VectorXd ph(GRIDSIZE);
	Eigen::VectorXd pf(GRIDSIZE);
	const bool HOME = 1;	
	const bool FOREIGN = 0;
  for (int k=0;k<GRIDSIZE;k++) 
  {
    ph[k] = prod(distr.h[k],HOME);
    pf[k] = prod(distr.f[k],FOREIGN);
  }

	//define time iteration variables
	Eigen::VectorXd val0(GRIDSIZE);
	Eigen::VectorXd val1(GRIDSIZE);
	Eigen::VectorXd diff(GRIDSIZE);
	Eigen::VectorXd fth(GRIDSIZE); 
	Eigen::VectorXd ftf(GRIDSIZE); 
	Eigen::VectorXd cfth(GRIDSIZE); 
	Eigen::VectorXd cftf(GRIDSIZE); 
	Eigen::VectorXd Sh(GRIDSIZE); 
	Eigen::VectorXd Sf(GRIDSIZE); 
	Eigen::VectorXd stayhome(GRIDSIZE);
	Eigen::VectorXd goabroad(GRIDSIZE);
	Eigen::MatrixXd A(GRIDSIZE,GRIDSIZE);
  Eigen::MatrixXd B(GRIDSIZE,GRIDSIZE);	
  Eigen::MatrixXd C(GRIDSIZE,GRIDSIZE);	
	Eigen::VectorXd b(GRIDSIZE);
	bool found;
	double mig_cut;
	double scale;

  //begin time iteration
  for (int t=1;t<PERIODS;t++)
  {
    cout << "NOW ON ITERATION " << t << endl;
    //read in new val0 and val1
    for (int k=0;k<GRIDMAX;k++)
    {
      val0[k] = distr.val_dyn(k,t-1);
      val1[k] = 0;
      diff[k] = val0[k];
    }

    //value iteration  
    while (diff.dot(diff)/val0.dot(val0)>1e-5)
    {

      //cumulative first terms
      fth[0] = val0[0] * distr.hd_dyn(0,t-1); 
      ftf[0] = val0[0] * distr.fd_dyn(0,t-1); 
      for (int k=1;k<GRIDSIZE;k++)	
      {
        fth[k] = fth[k-1] + val0[k] * distr.hd_dyn(k,t-1);
        ftf[k] = ftf[k-1] + val0[k] * distr.fd_dyn(k,t-1);
      }			
      cumsum(fth,cfth);
      cumsum(ftf,cftf);

      //get integral approxs
      Sh[0] = Sf[0] = 0.0; 
      for (int k=1;k<GRIDSIZE;k++) 
      {
        Sh[k] = cfth[k] - val0[k] * distr.chd_dyn(k,t-1);
        Sf[k] = cftf[k] - val0[k] * distr.cfd_dyn(k,t-1);
      }			

      //find cutoff for migration
      mig_cut = GRIDSIZE - 1;
      found = 0;
      for (int k=0;k<GRIDSIZE;k++)
      {
        stayhome[k] = ph[k] + HALP * Sh[k];
        goabroad[k] = pf[k] + FALP * Sf[k];
        if ((found == 0) & (stayhome[k] < goabroad[k]))
        {
          found = 1;
          mig_cut = k;
        }
      }
      cout << "The cut_off is currently at " << distr.h[mig_cut] << endl << endl;	
      // with cut-off in hand, solve for value function 
      // using the (really cool!) Lucas Moll matrix method.
      for (int k=0;k<GRIDSIZE;k++)
      {
        if (k<mig_cut)
        {
					B(k,k) = (RHO-THETA*HALP+HALP*distr.chd_dyn(k,t-1)+HALP*distr.h[k]/INCREMENT);
					for (int m=0;m!=k+1;m++)
						C(k,m) = distr.hd_dyn(m,t-1) * HALP * INCREMENT;
					b[k] = ph[k];
					if (k<GRIDSIZE-1)
						B(k,k+1) = -HALP * distr.h[k] / INCREMENT;
				}
				else
				{	
					B(k,k) = (RHO-THETA*HALP+FALP*distr.cfd_dyn(k,t-1)+HALP*distr.h[k]/INCREMENT);
					for (int m=0;m!=k+1;m++)
						C(k,m) = distr.fd_dyn(m,t-1) * FALP * INCREMENT;
					b[k] = pf[k];
					if (k<GRIDSIZE-1)
						B(k,k+1) = -FALP * distr.f[k] / INCREMENT;
				}
			}

			A = B-C;
			A.llt().solveInPlace(b);
			val1 = b;
			//val1 = A.inverse()*b; //Alternative (direct inversion)
			diff = val1-val0;	
			val0 = val1; //Read in old value
		}
		//read in new val0 and val1
		for (int k=0;k<GRIDMAX;k++)
		{
			distr.val_dyn(k,t) = val0[k];
		}
    distr.cut_off[t] = mig_cut;
		
		//update distributions, given cut-off THINK ABOUT hd[0]
			
		//calculate the distribution 
		for (int k=0;k<distr.cut_off[t];k++)
		{
			distr.hd_dyn(k,t) = distr.hd_dyn(k,t-1) + PERIOD_LENGTH * 
			  (-HALP * distr.hd_dyn(k,t-1) * distr.chd_dyn(k,t-1) + 
			   HALP * distr.hd_dyn(k,t-1) * 
			   (distr.chd_dyn(distr.cut_off[t],t-1)-distr.chd_dyn(k,t-1)) 
			  + FALP * (1-distr.chd_dyn(distr.cut_off[t],t-1)) * distr.fd_dyn(k,t-1));
			if (k!=0)
			  distr.chd_dyn(k,t) = distr.chd_dyn(k-1,t) + distr.hd_dyn(k,t) * INCREMENT;
      else 
        distr.chd_dyn(k,t) = 0;
		}

		for (int k=distr.cut_off[t];k<GRIDSIZE;k++)
		{
			distr.hd_dyn(k,t) = distr.hd_dyn(k,t-1) + PERIOD_LENGTH *
			  (-FALP * distr.hd_dyn(k,t-1) * distr.cfd_dyn(k,t-1)  +
			  FALP * distr.fd_dyn(k,t-1) * (1 - distr.chd_dyn(k,t-1)));
  			distr.chd_dyn(k,t) = distr.chd_dyn(k-1,t) + distr.hd_dyn(k,t) * INCREMENT;
	  }

		//rescale the distribution to get correct mass 
		scale = distr.chd_dyn(GRIDSIZE-1,t);
		if (scale>0)
		{
		 	// for (int k=mig_cut;k<GRIDSIZE;k++)
		 	for (int k=0;k<GRIDSIZE;k++)
		 	{
				distr.hd_dyn(k,t) = distr.hd_dyn(k,t) / scale;
		 		distr.chd_dyn(k,t) = distr.chd_dyn(k,t) / scale;
		 	}
		}
	}  //end time loop

	//WRITE TO FILE
	ofstream hd_dyn_stream;
	ofstream fd_dyn_stream;
	ofstream chd_dyn_stream;
	ofstream cfd_dyn_stream;
	hd_dyn_stream.open("hd_dyn.csv");
	fd_dyn_stream.open("fd_dyn.csv");
	chd_dyn_stream.open("chd_dyn.csv");
	cfd_dyn_stream.open("cfd_dyn.csv");
  hd_dyn_stream << distr.hd_dyn;  
  fd_dyn_stream << distr.fd_dyn;
  chd_dyn_stream << distr.chd_dyn;  
  cfd_dyn_stream << distr.cfd_dyn;
	
}

