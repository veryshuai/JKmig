#include <iostream>
#include <Eigen/Dense>
#include "mig.h"

using namespace std;

double THETA;
double B;
int GRIDSIZE;
double INCREMENT;
double cost;
double home;
Eigen::VectorXd input;
Eigen::VectorXd output;

//calculates production given cost, whether you live at home, and global parameters THETA and B.
double prod (double cost, bool home)
{
	if (home == 1)
		return (pow(cost,THETA));
	if (home == 0)	
		return (pow(cost,THETA)/(1.0+B));
	else
		cout << "Error in production calculation" << endl;
}

//running summation of Eigen vector
void cumsum (Eigen::VectorXd& input, Eigen::VectorXd& output)
{
	double summation = 0.0;

	for (int m=0;m<GRIDSIZE;m++)
	{
		summation += input[m] * INCREMENT;
		output[m] = summation;
	}
}
