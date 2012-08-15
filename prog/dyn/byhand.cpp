// byhand.cpp
// This program writes by hand distributions for use by dyn.cpp 

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

int main() {
  //output stuff
  ofstream h_stream;
  ofstream f_stream;
  ofstream hd_stream;
  ofstream fd_stream;
  ofstream chd_stream;
  ofstream cfd_stream;
  ofstream mig_cut_stream;
  ofstream val_stream;
  h_stream.open("h.csv");
  f_stream.open("f.csv");
  hd_stream.open("hd.csv");
  fd_stream.open("fd.csv");
  chd_stream.open("chd.csv");
  cfd_stream.open("cfd.csv");
  mig_cut_stream.open("cut.csv");
  val_stream.open("val.csv");

  //Initialize ALL required variables
  Eigen::VectorXd h(constants::GRIDSIZE); //grid values
  Eigen::VectorXd f(constants::GRIDSIZE);
  Eigen::VectorXd fd(constants::GRIDSIZE); //grid densities
  Eigen::VectorXd hd(constants::GRIDSIZE);
  Eigen::VectorXd hd_old(constants::GRIDSIZE);
  Eigen::VectorXd chd(constants::GRIDSIZE); //cumulative distributions
  Eigen::VectorXd cfd(constants::GRIDSIZE); //cumulative distributions
  Eigen::VectorXd val(constants::GRIDSIZE); //cumulative distributions

  for (int k=0;k<constants::GRIDSIZE;k++)
  {
    h[k]  =  constants::INCREMENT * (k+1);
    f[k]  =  constants::INCREMENT * (k+1);
    hd[k] =  constants::HLAM / pow(1 + constants::HLAM * h[k],2);
    fd[k] =  constants::FLAM / pow(1 + constants::FLAM * f[k],2);
    if (k>0) {
      chd[k] = chd[k-1] + hd[k] * constants::INCREMENT;
      cfd[k] = cfd[k-1] + fd[k] * constants::INCREMENT;
    }
    else {
      chd[k] = hd[k] * constants::INCREMENT;
      cfd[k] = fd[k] * constants::INCREMENT;
    }
  }

  //scale
  double hscale = chd[constants::GRIDSIZE-1];
  double fscale = cfd[constants::GRIDSIZE-1];
  for (int k=0;k<constants::GRIDSIZE;k++) {
    hd[k]  = hd[k] / hscale;
    fd[k]  = fd[k] / fscale;
    chd[k] = chd[k] / hscale;
    cfd[k] = cfd[k] / fscale;
  }

  const double DENOM = constants::RHO + constants::THETA * constants::HALP;
  for (int k=1;k<constants::GRIDSIZE;k++)
  {
    val[k] = pow(h[k],constants::THETA)/DENOM;
  }

  //write to files
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
    val_stream << val[m] << endl;
  }
  mig_cut_stream << constants::GRIDSIZE-1 << endl;

  return(0);
}

