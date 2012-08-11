#include <iostream>
#include <Eigen/Dense>
#include "mig.h"

double cost;
double home;
Eigen::VectorXd input;
Eigen::VectorXd output;

//calculates production given cost, whether you live at home, and global parameters THETA and B.
double prod (double cost, bool home)
{
  if (home == 1)
    return (pow(cost,constants::THETA));
  if (home == 0)
    return (pow(cost,constants::THETA)/(1.0+constants::B));
  else
  {
    std::cout << "Error in production calculation" << std::endl;
    return 0;
  }
}

//"integration" of Eigen vector
void cumsum (Eigen::VectorXd& input, Eigen::VectorXd& output)
{
  double summation = 0;

  for (int m=0;m<constants::GRIDSIZE;m++)
  {
    summation += input[m] * constants::INCREMENT;
    output[m] = summation;
  }
}
