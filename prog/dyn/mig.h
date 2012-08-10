// This is my first header file.  It contains several functions which I use through the sequence of c++ programs for use with my migration growth project.

#ifndef MIG_H
#define MIG_H
 
namespace constants{
  const double  HALP            = 0.02;   //home meeting hazard
  const double  FALP            = 0.015;  //foreign meeting hazard
  const int     GRIDSIZE        = 1000;   //number of grid points
  const int     GRIDMAX         = 300;    //maximum grid value
  const int     PERIODS         = 10;    //number of periods
  const double  PERIOD_LENGTH   = 0.2;  //length of a period (say a unit is a year)
  const double  B               = 1;    //penalty for living abroad
  const double  THETA           =-0.5;    //production function parameter
  const double  RHO             = 0.05;   //time discount factor
  const double  FLAM            = 1;      //foreign distribution lambda
  const double  HLAM            = 0.35;    //home distribution lambda
  const int     MIG_INIT        = 100;    //initial cut off guess
  const double  INCREMENT = (double) GRIDMAX/ (double) GRIDSIZE; //step size for the grid
}

double prod (double cost, bool home);
void cumsum (Eigen::VectorXd& input, Eigen::VectorXd& output);

#endif