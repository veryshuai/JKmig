//dyn.cpp--this program winds the clock forward on the distribution solved for in comb.cpp.

#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include "mig.h"
#include <sys/stat.h> //for file existance

using namespace std;

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
  Eigen::VectorXd hGDP;
  Eigen::VectorXd fGDP;
  Eigen::VectorXd hGrowth;
  Eigen::VectorXd fGrowth;

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
    hGDP    =  Eigen::VectorXd(periods);
    fGDP    =  Eigen::VectorXd(periods);
    hGrowth =  Eigen::VectorXd(periods);
    fGrowth =  Eigen::VectorXd(periods);
  }
};

struct ValVars {
  Eigen::VectorXd vdiff;
  Eigen::VectorXd fth;
  Eigen::VectorXd ftf;
  Eigen::VectorXd cfth;
  Eigen::VectorXd cftf;
  Eigen::VectorXd Sh;
  Eigen::VectorXd Sf;
  Eigen::VectorXd stayhome;
  Eigen::VectorXd goabroad;
  Eigen::MatrixXd A;
  Eigen::MatrixXd B;
  Eigen::MatrixXd C;
  Eigen::VectorXd b;
  bool found;
  double mig_cut;

  ValVars(int gridSize, int periods) {
    vdiff    =  Eigen::VectorXd(gridSize);
    fth      =  Eigen::VectorXd(gridSize);
    ftf      =  Eigen::VectorXd(gridSize);
    cfth     =  Eigen::VectorXd(gridSize);
    cftf     =  Eigen::VectorXd(gridSize);
    Sh       =  Eigen::VectorXd(gridSize);
    Sf       =  Eigen::VectorXd(gridSize);
    stayhome =  Eigen::VectorXd(gridSize);
    goabroad =  Eigen::VectorXd(gridSize);
    A        =  Eigen::MatrixXd(gridSize,gridSize);
    B        =  Eigen::MatrixXd(gridSize,gridSize);
    C        =  Eigen::MatrixXd(gridSize,gridSize);
    b        =  Eigen::VectorXd(gridSize);
  }
};

bool fileExists(const std::string& filename)
{
    struct stat buf;
    if (stat(filename.c_str(), &buf) != -1)
    {
        return true;
    }
    return false;
}

bool read (Eigen::VectorXd& target, const char* filename, Eigen::MatrixXd& target2) {
  bool exist;
  exist = fileExists(filename);
  if (exist == 1){
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
  return exist;
}

void write (Eigen::VectorXd& input, const char* filename) {
  ofstream writeme;
  writeme.open(filename);
  writeme << input;
}
void write (Eigen::MatrixXd& input, const char* filename) {
  ofstream writeme;
  writeme.open(filename);
  writeme << input;
}

void dynFor (Eigen::MatrixXd& fd_dyn, Eigen::MatrixXd& cfd_dyn, Eigen::VectorXd& f){
  Eigen::VectorXd fd(constants::GRIDSIZE);
  Eigen::VectorXd cfd(constants::GRIDSIZE);
  Eigen::VectorXd interppoint(constants::GRIDSIZE);
  Eigen::MatrixXd interpweights(constants::GRIDSIZE,2);
  for (int t=1;t<constants::PERIODS;t++)
  {
    double newspot = exp(constants::HALP * constants::PERIOD_LENGTH * t);
    for (int k=0;k<constants::GRIDSIZE;k++)
    {
      interppoint[k] = 0;
      for (int m=1;m<constants::GRIDSIZE;m++)
      {
        if (f[m]>=newspot*f[k])
        {
          interppoint[k] = m;
          interpweights(k,0) =  (f[m] - newspot*f[k]) / constants::INCREMENT;
          interpweights(k,1) = 1 - interpweights(k,0);
          break;
        }
      }
    }
    // do the interpolation
    for (int k=0;k<constants::GRIDSIZE;k++)
    {
      if (interppoint[k]!=0)
      {
        fd_dyn(k,t) = interpweights(k,0) * fd_dyn(interppoint[k]-1,0)
          + interpweights(k,1) * fd_dyn(interppoint[k],0); //scale over
        fd_dyn(k,t) = newspot * fd_dyn(k,0); //scale up
      }
      else
        fd_dyn(k,t) = fd_dyn(constants::GRIDSIZE-1,0); //things don't change after a while.
      if (k!=0){
        cfd_dyn(k,t) = cfd_dyn(k-1,t) + fd_dyn(k,t) * constants::INCREMENT;
      }
      else
        cfd_dyn(k,t) = fd_dyn(k,t) * constants::INCREMENT;
    }
  }
}

void distCalc (int t, Eigen::MatrixXd& hd_dyn, Eigen::MatrixXd& fd_dyn,
    Eigen::MatrixXd& chd_dyn, Eigen::MatrixXd& cfd_dyn, Eigen::VectorXd& cut_off){
  for (int k=0;k<cut_off[t];k++)
  {
    hd_dyn(k,t) = hd_dyn(k,t-1) + constants::PERIOD_LENGTH *
      (-constants::HALP * hd_dyn(k,t-1) * chd_dyn(k,t-1) +
       constants::HALP * hd_dyn(k,t-1) * (chd_dyn(cut_off[t],t-1)-chd_dyn(k,t-1))
       + constants::FALP * (1-chd_dyn(cut_off[t],t-1)) * fd_dyn(k,t-1));
    if (k!=0)
      chd_dyn(k,t) = chd_dyn(k-1,t) + hd_dyn(k,t) * constants::INCREMENT;
    else
      chd_dyn(k,t) = hd_dyn(k,t) * constants::INCREMENT;
  }

  for (int k=cut_off[t];k<constants::GRIDSIZE;k++)
  {
    hd_dyn(k,t) = hd_dyn(k,t-1) + constants::PERIOD_LENGTH *
      (-constants::FALP * hd_dyn(k,t-1) * cfd_dyn(k,t-1)  +
       constants::FALP * fd_dyn(k,t-1) * (1 - chd_dyn(k,t-1)));
    chd_dyn(k,t) = chd_dyn(k-1,t) + hd_dyn(k,t) * constants::INCREMENT;
  }
}

int valIter (int t, Eigen::VectorXd& val0, Eigen::VectorXd& val1, Eigen::VectorXd ph,
    Eigen::VectorXd pf, Eigen::MatrixXd& val_dyn, Eigen::MatrixXd& hd_dyn,
    Eigen::MatrixXd& fd_dyn, Eigen::MatrixXd& chd_dyn, Eigen::MatrixXd& cfd_dyn,
    Eigen::VectorXd& cut_off, Eigen::VectorXd& h, Eigen::VectorXd& f){
  ValVars ti(constants::GRIDSIZE,constants::PERIODS);

  ti.vdiff = val0;
  while (ti.vdiff.dot(ti.vdiff)/val0.dot(val0)>1e-5)
  {
    //cumulative first terms
    ti.fth[0] = val0[0] * hd_dyn(0,t-1);
    ti.ftf[0] = val0[0] * fd_dyn(0,t-1);
    for (int k=1;k<constants::GRIDSIZE;k++)
    {
      ti.fth[k] = val0[k] * hd_dyn(k,t-1);
      ti.ftf[k] = val0[k] * fd_dyn(k,t-1);
    }
    cumsum(ti.fth,ti.cfth);
    cumsum(ti.ftf,ti.cftf);

    //get integral approxs
    ti.Sh[0] = ti.Sf[0] = 0.0; 
    for (int k=1;k<constants::GRIDSIZE;k++) 
    {
      ti.Sh[k] = ti.cfth[k] - val0[k] * chd_dyn(k,t-1);
      ti.Sf[k] = ti.cftf[k] - val0[k] * cfd_dyn(k,t-1);
    }

    //find cutoff for migration
    ti.mig_cut = 0;
    ti.found = 0;
    for (int k=constants::GRIDSIZE-1;k>-1;k--)
    {
      ti.stayhome[k] = ph[k] + constants::HALP * ti.Sh[k];
      ti.goabroad[k] = pf[k] + constants::FALP * ti.Sf[k];
      if ((ti.found == 0) & (ti.stayhome[k] > ti.goabroad[k]))
      {
        ti.found = 1;
        ti.mig_cut = k;
      }
    }
    cout << "The cut_off is currently at " << h[ti.mig_cut] << endl << endl;
    // with cut-off in hand, solve for value function
    // using the (really cool!) Lucas Moll matrix method.
    ti.A.setZero(); ti.B.setZero(); ti.C.setZero();
    for (int k=0;k<constants::GRIDSIZE;k++)
    {
      if (k<ti.mig_cut)
      {
        ti.B(k,k) = (constants::RHO+constants::THETA*constants::HALP
            +constants::HALP*chd_dyn(k,t-1)+constants::HALP*h[k]/constants::INCREMENT);
        for (int m=0;m<k+1;m++)
          ti.C(k,m) = hd_dyn(m,t-1) * constants::HALP * constants::INCREMENT;
        ti.b[k] = ph[k];
        if (k<constants::GRIDSIZE-1)
          ti.B(k,k+1) = -constants::HALP * h[k] / constants::INCREMENT;
      }
      else
      {
        ti.B(k,k) = (constants::RHO+constants::THETA*constants::HALP
            +constants::FALP*cfd_dyn(k,t-1)+constants::HALP*h[k]/constants::INCREMENT);
        for (int m=0;m<k+1;m++)
          ti.C(k,m) = fd_dyn(m,t-1) * constants::FALP * constants::INCREMENT;
        ti.b[k] = pf[k];
        if (k<constants::GRIDSIZE-1)
          ti.B(k,k+1) = -constants::HALP * f[k] / constants::INCREMENT;
      }
    }
    ti.A = ti.B-ti.C;
    // ti.A.llt().solveInPlace(ti.b); //Alternative, much faster, but unreliable
    val1 = ti.b;
    val1 = ti.A.inverse()*ti.b; //slow as !@#@ direct inversion
    ti.vdiff = val1-val0;
    val0 = val1; //Read in old value
  }
  //read in new val0 and val1
  for (int k=0;k<constants::GRIDMAX;k++)
  {
    val_dyn(k,t) = val0[k];
  }
  return ti.mig_cut;
}

void dynHom (Eigen::MatrixXd& val_dyn, Eigen::MatrixXd& hd_dyn, Eigen::MatrixXd& fd_dyn,
    Eigen::MatrixXd& chd_dyn, Eigen::MatrixXd& cfd_dyn, Eigen::VectorXd& cut_off,
    Eigen::VectorXd& h, Eigen::VectorXd& f){
  Eigen::VectorXd ph(constants::GRIDSIZE);
  Eigen::VectorXd pf(constants::GRIDSIZE);
  Eigen::VectorXd val0(constants::GRIDSIZE);
  Eigen::VectorXd val1(constants::GRIDSIZE);

  //get productions
  for (int k=0;k<constants::GRIDSIZE;k++)
  {
    ph[k] = prod(h[k],1);
    pf[k] = prod(f[k],0);
  }

  //begin time iteration
  for (int t=1;t<constants::PERIODS;t++)
  {
    cout << "NOW ON ITERATION " << t << endl;

    //read in new initial val0
    for (int k=0;k<constants::GRIDSIZE;k++)
    {
      val0[k] = val_dyn(k,t-1);
    }
    //value iteration
    cut_off[t] = valIter(t,val0,val1,ph,pf,val_dyn,hd_dyn,fd_dyn,
        chd_dyn,cfd_dyn, cut_off,h,f);

    //update distributions, given cut-off
    distCalc(t,hd_dyn,fd_dyn,chd_dyn,cfd_dyn,cut_off);
  }
}

void growth (Eigen::VectorXd& hGDP, Eigen::VectorXd& fGDP, Eigen::VectorXd& hGrowth,
    Eigen::VectorXd& fGrowth, Eigen::MatrixXd& hd_dyn, Eigen::MatrixXd& fd_dyn,
    Eigen::VectorXd& h, Eigen::VectorXd& f, Eigen::VectorXd& cut_off){
  //get productions
  Eigen::VectorXd ph(constants::GRIDSIZE);
  Eigen::VectorXd pf(constants::GRIDSIZE);
  for (int k=0;k<constants::GRIDSIZE;k++) {
    ph[k] = prod(h[k],1);
    pf[k] = prod(f[k],0);
  }

  for (int t=0;t<constants::PERIODS;t++) {
    hGDP[t] = ph[0] * hd_dyn(0,t) * constants::INCREMENT;
    fGDP[t] = ph[0] * fd_dyn(0,t) * constants::INCREMENT;
    for (int k=1;k<cut_off[t];k++) {
      hGDP[t] = hGDP[t] + ph[k] * hd_dyn(k,t) * constants::INCREMENT;
      fGDP[t] = fGDP[t] + ph[k] * fd_dyn(k,t) * constants::INCREMENT;
    }
    for (int k=cut_off[t];k<constants::GRIDSIZE;k++) {
      hGDP[t] = hGDP[t] + pf[k] * hd_dyn(k,t) * constants::INCREMENT;
      fGDP[t] = fGDP[t] + ph[k] * fd_dyn(k,t) * constants::INCREMENT;
    }
    if (t>0){
      hGrowth[t] = ((hGDP[t] - hGDP[t-1]) / hGDP[t-1]) / constants::PERIOD_LENGTH;
      fGrowth[t] = ((fGDP[t] - fGDP[t-1]) / fGDP[t-1]) / constants::PERIOD_LENGTH;
    }
  }
}

int main() {
  Distributions distr(constants::GRIDSIZE,constants::PERIODS);

  //read in distribution data
  Eigen::MatrixXd junk_mat(constants::GRIDSIZE,1);
  Eigen::VectorXd junk_vec(constants::GRIDSIZE);
  read(distr.h,"h.csv",junk_mat);
  read(distr.f,"f.csv",junk_mat);
  read(distr.hd,"hd.csv",distr.hd_dyn);
  read(distr.fd,"fd.csv",distr.fd_dyn);
  read(distr.chd,"chd.csv",distr.chd_dyn);
  read(distr.cfd,"cfd.csv",distr.cfd_dyn);
  read(distr.cut_off,"cut.csv",junk_mat);
  read(junk_vec,"val.csv",distr.val_dyn);

  // get foreign distribution dynamics
  // dynFor(distr.fd_dyn,distr.cfd_dyn,distr.f);
  Eigen::VectorXd fake_cut;
  fake_cut.setConstant(constants::GRIDSIZE,constants::GRIDSIZE-1);
  for (int t=1;t<constants::PERIODS;t++){
    distCalc(t,distr.fd_dyn,distr.fd_dyn,distr.cfd_dyn,distr.cfd_dyn,fake_cut);
  }

  // get home distribution dynamics
  dynHom(distr.val_dyn,distr.hd_dyn,distr.fd_dyn,
      distr.chd_dyn,distr.cfd_dyn,distr.cut_off,distr.h,distr.f);

  // get growth
  growth(distr.hGDP,distr.fGDP,distr.hGrowth,distr.fGrowth,distr.hd_dyn,
      distr.fd_dyn,distr.h,distr.f,distr.cut_off);

  //WRITE TO FILE
  write(distr.hd_dyn,"hd_dyn.csv");
  write(distr.fd_dyn,"fd_dyn.csv");
  write(distr.chd_dyn,"chd_dyn.csv");
  write(distr.cfd_dyn,"cfd_dyn.csv");
  write(distr.hGDP,"hGDP.csv");
  write(distr.fGDP,"fGDP.csv");
  write(distr.hGrowth,"hGrowth.csv");
  write(distr.fGrowth,"fGrowth.csv");
  write(distr.val_dyn,"val_dyn.csv");
  write(distr.cut_off,"cut_dyn.csv");

}

