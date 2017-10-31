#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_



#include "math.h"
#include "lines.h"

#include "stdio.h"
#include "stdlib.h"
#include <vector>

#include "Param.h"

#include <algorithm>
using namespace std;



bool goToNextCharInFile(FILE *fp, char c);

int binomialCoeff(int n, int k);

// clustering using tanh
double stretchingTanh(double t, double C, double E);

// clustering using atanh
double stretchingAtanh(double t, double C, double E);

void firstLengthDistr(double fL, int np, double *distr);

void firstLastLengthDistr(double fL, double lL, int np, double *distr);

//void arbitraryClustering(vector<Param> *distrPtsFile, SPLINE &distr);


//void curvatureDistr(SPLINE edge, int np, double *distr);

#endif
