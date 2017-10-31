/*
 * functions.cpp
 *
 *  Created on: Mar 20, 2012
 *      Author: enrico
 */
#include "functions.h"

bool goToNextCharInFile(FILE *fp, char c)
{
  char c_compare;

  while (EOF != fscanf(fp,"%c",&c_compare))
  {
//    printf("%c",c_compare);
    if(c_compare == c) return true;
  }

  return false;
}

int binomialCoeff(int n, int k)
{
  if (k > n - k) k = n-k;
  int c = 1;
  for (int i=0;i<k;i++)
  {
    c = c*(n-(k-(i+1)));
    c = c/(i+1);
  }
  return c;
}


// clustering used for wall, uses tanh
double stretchingTanh(double t, double C, double E)
{
  double y, y1, y2;
  E = min( max(0.0,E), 1.0) * 2.0;

  y  = 1.0 + tanh(C*(2.0*t-E));
  y1 = 1.0 + tanh(C*(2.0*0.0-E));
  y2 = 1.0 + tanh(C*(2.0*1.0-E));

  y = (y-y1)/(y2-y1);

  return(y);
}


// clustering used for wall, uses atanh
double stretchingAtanh(double t, double C, double E)
{
  double y, y1, y2;
  E = min( max(0.0,E), 1.0) * 2.0;

  y  = atanh(C*(2.0*t-E));
  y1 = atanh(C*(2.0*(0.0)-E));
  y2 = atanh(C*(2.0*1.0-E));

  y = (y-y1)/(y2-y1);

  return(y);
}



void firstLengthDistr(double fL, int np, double *distr)
{
	double toll = 1e-6, err = 1.;
	int nmax = 50, k = 0;
	double R = 1.5, dR;
	double f, df;

	int nInt = np-1;

	while ((err>toll) && (k<nmax))
	{
      f = -1./fL;
      df = 0.0;

      for (int i=1; i<nInt+1; i++)
          f += pow(R, i-1);
      for (int i=2; i<nInt+1; i++)
          df = df + (i-1)*pow(R, i-2);

      err = fabs(f/df) / R;
      k++;

      R = R - f/df;
	}

	double dL=fL;
	distr[0] = 0.0;
	distr[1] = fL;

	for (int i=2; i<np; i++)
	{
		distr[i] = distr[i-1] + dL*R;
		dL = distr[i] - distr[i-1];
	}
}


void firstLastLengthDistr(double fL, double lL, int np, double *distr)
{
	if (np%2!=0)
	{

		int k = 0, nmax = 100;
		double toll = 1e-10, diff = 1.;
		double a = 0.5;
		double s1, s2=1.;

		int npHalf = (np+1)/2;
		double distrHalfL[npHalf], distrHalfR[npHalf];

		while ((fabs(diff)>toll) && (k<nmax) && (a<1.0))
		{
			firstLengthDistr(fL/a, npHalf, distrHalfL);
			s2 = (1-a);
			firstLengthDistr(lL/s2, npHalf, distrHalfR);

			diff = a*(distrHalfL[npHalf-1]-distrHalfL[npHalf-2]) - s2*(distrHalfR[npHalf-1]-distrHalfR[npHalf-2]);
			a += -diff;
			k++;
		}
		for (int i=0; i<npHalf; i++)
			distr[i] = a*distrHalfL[i];
		for (int i=npHalf; i>0; i--)
			distr[npHalf+(npHalf-i)] = (1-s2*distrHalfR[i-2]);
	}
	else
	{
		int k = 0, nmax = 100;
		double toll = 1e-10, diff = 1.;
		double a = 0.5;
		double s1, s2=1.;

		int npHalf = np/2+1;
		double distrHalfL[npHalf], distrHalfR[npHalf];

		while ((fabs(diff)>toll) && (k<nmax) && (a<1.0))
		{
			firstLengthDistr(fL/a, npHalf, distrHalfL);
			s2 = (1-a)+a*(distrHalfL[npHalf-1]-distrHalfL[npHalf-2]);
			firstLengthDistr(lL/s2, npHalf, distrHalfR);

			diff = a*(distrHalfL[npHalf-1]-distrHalfL[npHalf-2]) - s2*(distrHalfR[npHalf-1]-distrHalfR[npHalf-2]);
			a += -diff;
			k++;
		}
		for (int i=0; i<npHalf; i++)
			distr[i] = a*distrHalfL[i];
		for (int i=npHalf; i>1; i--)
			distr[npHalf-1+(npHalf-i)] = (1-s2*distrHalfR[i-2]);
	}
}



