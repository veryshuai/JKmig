// comb.cpp
// This program solves for a balanced growth path in the Lucas style migration model

#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/Cholesky>
#include "mig.h"

using namespace std;

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
  ofstream val_stream;
  ofstream test_stream;
  h_stream.open("h.csv");
  f_stream.open("f.csv");
  hd_stream.open("hd.csv");
  fd_stream.open("fd.csv");
  chd_stream.open("chd.csv");
  cfd_stream.open("cfd.csv");
  val_stream.open("val.csv");
  mig_cut_stream.open("cut.csv");
  test_stream.open("test.csv");

  //Initialize ALL required variables
  Eigen::VectorXd h(constants::GRIDSIZE); //grid values
  Eigen::VectorXd f(constants::GRIDSIZE);
  Eigen::VectorXd fd(constants::GRIDSIZE); //grid densities
  Eigen::VectorXd hd(constants::GRIDSIZE);
  Eigen::VectorXd hd_old(constants::GRIDSIZE);
  Eigen::VectorXd chd(constants::GRIDSIZE); //cumulative distributions
  Eigen::VectorXd cfd(constants::GRIDSIZE);
  Eigen::VectorXd val0(constants::GRIDSIZE); //old value function
  Eigen::VectorXd val1(constants::GRIDSIZE); //new value function
  Eigen::VectorXd fth(constants::GRIDSIZE); //initialize first term placeholders for value function approximation
  Eigen::VectorXd ftf(constants::GRIDSIZE);
  Eigen::VectorXd cfth(constants::GRIDSIZE); //cumulative first terms for value function approximation
  Eigen::VectorXd cftf(constants::GRIDSIZE);
  Eigen::VectorXd Sh(constants::GRIDSIZE); //value function integral approx	
  Eigen::VectorXd Sf(constants::GRIDSIZE);
  Eigen::VectorXd ph(constants::GRIDSIZE);
  Eigen::VectorXd pf(constants::GRIDSIZE);
  Eigen::MatrixXd A(constants::GRIDSIZE,constants::GRIDSIZE);
  Eigen::MatrixXd B(constants::GRIDSIZE,constants::GRIDSIZE);
  Eigen::MatrixXd C(constants::GRIDSIZE,constants::GRIDSIZE);
  Eigen::VectorXd b(constants::GRIDSIZE);
  Eigen::VectorXd stayhome(constants::GRIDSIZE);
  Eigen::VectorXd goabroad(constants::GRIDSIZE);
  Eigen::VectorXd diff(constants::GRIDSIZE);
  Eigen::VectorXd dist_diff(constants::GRIDSIZE);
  const bool home = 1;
  const bool foreign = 0;
  int mig_cut = constants::MIG_INIT;
  int mig_cut_old = constants::MIG_INIT+1;
  bool found = 0;
  int counter;
  double dynamicslodown;

  double barea_diff = 1;
  double barea = 0.5;
  double barea_old = barea;
  double scale;

  for (int k=0;k<constants::GRIDSIZE;k++)
  {
    h[k]  =  constants::INCREMENT * (k+1);
    f[k]  =  constants::INCREMENT * (k+1);
    fd[k] =  constants::FLAM / pow(1 + constants::FLAM * f[k],2);
    hd[k] =  fd[k];
    if (k>0)
      cfd[k] = cfd[k-1] + fd[k] * constants::INCREMENT;
    else
      cfd[k] = fd[k] * constants::INCREMENT;
  }

  //Initialize value function
  const double DENOM = constants::RHO + constants::THETA * constants::HALP;
  for (int k=1;k<constants::GRIDSIZE;k++)
  {
    val0(k) = pow(h[k],constants::THETA)/DENOM;
  }
  val0[0] = val0[1];
  val1 = val0.array() + 1;

  hd_old = hd.array() + 1;
  dist_diff = hd - hd_old;
  while (dist_diff.dot(dist_diff)/hd_old.dot(hd_old)>1e-6 || mig_cut_old != mig_cut)
  {
    hd_old = hd;
    mig_cut_old = mig_cut;
    barea_diff = 1;

    //calculate new distribution. Stopping rule is area under mig_cut must be consistent with eq.
    while  (barea_diff>1e-5)
    {
      //put in proper hd[0]
      hd[0] = constants::HLAM + constants::INCREMENT / constants::HLAM
        * (-constants::HLAM + constants::HLAM * barea + constants::FLAM
        * (1-barea) * constants::FALP / constants::HALP);
      hd[0] = min(hd[0],fd[0]);
      chd[0] = hd[0] * constants::INCREMENT;

      //loop to calculate the distribution
      for (int k=0;k<mig_cut;k++)
      {
        hd[k+1] = hd[k] + constants::INCREMENT / h[k]
          * (-hd[k] - hd[k] * chd[k] + hd[k] * (barea - chd[k]) + fd[k]
              * (1-barea) * constants::FALP / constants::HALP);
        if (hd[k+1] < 0){
          hd[k+1] = 0;
          cout << "Returning negative probability densities" << endl;
        }
        chd[k+1] = chd[k] + hd[k+1] * constants::INCREMENT;
        if (chd[k+1] > cfd[k+1]){
          hd[k+1] = fd[k+1];
          chd[k+1] = chd[k] + hd[k+1] * constants::INCREMENT;
          cout << "Foreign distribution stoch. dominance correction" << endl;
        }
      }

      for (int k=mig_cut+1;k<constants::GRIDSIZE;k++)
      {
        hd[k] = hd[k-1] + constants::INCREMENT / (h[k-1])
          * (-hd[k-1] - hd[k-1] * constants::FALP / constants::HALP
          * cfd[k-1] + fd[k-1] * (1 - chd[k-1]) * constants::FALP / constants::HALP);
        if (hd[k] < 0){
          hd[k] = 0;
          cout << "Returning negative probability densities" << endl;
        }
        chd[k] = chd[k-1] + hd[k] * constants::INCREMENT;
        if (chd[k] > cfd[k]){
          hd[k] = fd[k];
          chd[k] = chd[k-1] + hd[k] * constants::INCREMENT;
          cout << "Foreign distribution stoch. dominance correction" << endl;
        }
      }

      barea_old = barea;
      barea = min(chd[mig_cut],.99);
      barea_diff = abs(barea - barea_old);
      cout << "Inconsistency in area under mig_cut is " << barea_diff << endl;
    }
    cout << endl;

    //if mass > 1, rescale the distribution to get an allowable mass
    scale = chd[constants::GRIDSIZE-1];
    cout << scale << endl;

    //FIND THE NEW CUT-OFF

    diff = val0;
    counter = 1;
    dynamicslodown = 1;

    while (diff.dot(diff) / val0.dot(val0)>1e-5)
    {

      //cumulative first terms
      for (int k=0;k<constants::GRIDSIZE;k++)
      {
        fth[k] = val0[k] * hd[k];
        ftf[k] = val0[k] * fd[k];
      }
      cumsum(fth,cfth);
      cumsum(ftf,cftf);

      //get integral approxs
      for (int k=0;k<constants::GRIDSIZE;k++)
      {
        Sh[k] = cfth[k] - val0[k] * chd[k];
        Sf[k] = cftf[k] - val0[k] * cfd[k];
      }
      cout << "Expected Moving up" << endl;
      for (int k=0;k<10;k++)
        cout << constants::HALP * Sh[k] << "  " << constants::FALP * Sf[k] << endl;
      cout << endl;
      cout << "Cumdists" << endl;
      for (int k=0;k<10;k++)
        cout << chd[k] << "  " << cfd[k] << endl;
      cout << endl;

      //get productions
      for (int k=0;k<constants::GRIDSIZE;k++) 
      {
        ph[k] = prod(h[k],home);
        pf[k] = prod(f[k],foreign);
      }
      cout << "products" << endl;
      for (int k=0;k<10;k++)
        cout << ph[k] << "  " << pf[k] << endl;
      cout << endl;

      //find cutoff for migration
      mig_cut = constants::GRIDSIZE-1;
      found = 0;
      for (int k=0;k<constants::GRIDSIZE;k++)
      {
        stayhome[k] = ph[k] + constants::HALP * Sh[k];
        goabroad[k] = pf[k] + constants::FALP * Sf[k];
        // cout << stayhome[k] << "  " << goabroad[k] << endl;
        if ((found == 0) & (stayhome[k] < goabroad[k]))
        {
          found = 1;
          mig_cut = k;
        }
      }
      // for (int k=0;k<10;k++)
      //  cout << stayhome[k] << "  " << goabroad[k] << endl;
      cout << "The cut_off is currently at " << h[mig_cut] << endl << endl;
      // with cut-off in hand, solve for value function using the (really cool!) Lucas Moll matrix method.
      for (int k=0;k<constants::GRIDSIZE;k++)
      {
        if (k<mig_cut)
        {
          B(k,k) = (constants::RHO+constants::THETA*constants::HALP+constants::HALP*chd[k]+constants::HALP*h[k]/constants::INCREMENT);
          for (int m=0;m!=k+1;m++)
            C(k,m) = hd[m] * constants::HALP * constants::INCREMENT;
          b[k] = ph[k];
          if (k<constants::GRIDSIZE-1)
            B(k,k+1) = -constants::HALP * h[k] / constants::INCREMENT;
        }
        else
        {
          B(k,k) = constants::RHO+constants::THETA*constants::HALP+constants::FALP*cfd[k]+constants::HALP*h[k]/constants::INCREMENT;
          for (int m=0;m!=k+1;m++)
            C(k,m) = fd[m] * constants::FALP * constants::INCREMENT;
          b[k] = pf[k];
          if (k<constants::GRIDSIZE-1)
            B(k,k+1) = -constants::HALP * h[k] / constants::INCREMENT;
        }
      }

      A = B-C;
      val1 = A.inverse()*b; //slow direct inversion
      //A.llt().solveInPlace(b);//Alternative much faster, but unreliable inversion
      //val1 = b;

      if (counter % 3 == 0)
        dynamicslodown = dynamicslodown + 1;
      cout << "Cycle break at " << dynamicslodown << endl;
      counter = counter + 1;
      val1 = (((double) dynamicslodown - 1) * val0 + val1) / dynamicslodown;
      diff = val1-val0;
      val0 = val1; //Read in old value
    }
    dist_diff = hd - hd_old;
  }

  //write to files
  for (int m=0;m<constants::GRIDSIZE;m++)
  {
    for (int k=0;k<constants::GRIDSIZE;k++)
    {
      test_stream << C(m,k) << ", ";
    }
    test_stream << endl;
  }
  for (int m=0;m<constants::GRIDSIZE;m++)
  {
    h_stream << h[m] << endl;
  }
  for (int m=0;m<constants::GRIDSIZE;m++)
  {
    f_stream << f[m] << endl;
  }
  for (int m=0;m<constants::GRIDSIZE;m++)
  {
    hd_stream << hd[m] << endl;
  }
  for (int m=0;m<constants::GRIDSIZE;m++)
  {
    fd_stream << fd[m] << endl;
  }
  for (int m=0;m<constants::GRIDSIZE;m++)
  {
    chd_stream << chd[m] << endl;
  }
  for (int m=0;m<constants::GRIDSIZE;m++)
  {
    cfd_stream << cfd[m] << endl;
  }
  for (int m=0;m<constants::GRIDSIZE;m++)
  {
    val_stream << val1[m] << endl;
  }
  mig_cut_stream << mig_cut << endl;

  return(0);
}

