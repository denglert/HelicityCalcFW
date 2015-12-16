#include <iostream>
#include <cmath>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TGenPhaseSpace.h>
#include <TCanvas.h>
#include "HelicityFW.h"
#include "FortranInterface.h"
#include "PhysConst.h"
#include "UtilFunctions.h"
#include "cuba.h"

double fxn(double *x)
{
	double fx = 6;
	return fx;
};

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#if REALSIZE == 16
#include "cubaq.h"
#elif REALSIZE == 10
#include "cubal.h"
#else
#include "cuba.h"
#endif


static inline cubareal Sq(cubareal x) {
  return x*x;
}


static int Integrand(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define x xx[0]
#define y xx[1]
#define z xx[2]
#define f ff[0]

#ifndef FUN
#define FUN 1
#endif

#define rsq (Sq(x) + Sq(y) + Sq(z))

#if FUN == 1
  // f = sin(x)*cos(y)*exp(z);
  f = 6;
#elif FUN == 2
  f = 1/(Sq(x + y) + .003)*cos(y)*exp(z);
#elif FUN == 3
  f = 1/(3.75 - cos(M_PI*x) - cos(M_PI*y) - cos(M_PI*z));
#elif FUN == 4
  f = fabs(rsq - .125);
#elif FUN == 5
  f = exp(-rsq);
#elif FUN == 6
  f = 1/(1 - x*y*z + 1e-10);
#elif FUN == 7
  f = sqrt(fabs(x - y - z));
#elif FUN == 8
  f = exp(-x*y*z);
#elif FUN == 9
  f = Sq(x)/(cos(x + y + z + 1) + 5);
#elif FUN == 10
  f = (x > .5) ? 1/sqrt(x*y*z + 1e-5) : sqrt(x*y*z);
#else
  f = (rsq < 1) ? 1 : 0;
#endif

  return 0;
}

/*********************************************************************/

#define NDIM 3
#define NCOMP 1
#define USERDATA NULL
#define NVEC 1
#define EPSREL 1e-3
#define EPSABS 1e-12
#define VERBOSE 2
#define LAST 4
#define SEED 0
#define MINEVAL 0
#define MAXEVAL 50000

#define NSTART 1000
#define NINCREASE 500
#define NBATCH 1000
#define GRIDNO 0
#define STATEFILE NULL
#define SPIN NULL

#define NNEW 1000
#define NMIN 2
#define FLATNESS 25.

#define KEY1 47
#define KEY2 1
#define KEY3 1
#define MAXPASS 5
#define BORDER 0.
#define MAXCHISQ 10.
#define MINDEVIATION .25
#define NGIVEN 0
#define LDXGIVEN NDIM
#define NEXTRA 0

#define KEY 0

int main() {
  int comp, nregions, neval, fail;
  cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];

#if 1
  printf("-------------------- Vegas test --------------------\n");

  Vegas(NDIM, NCOMP, Integrand, USERDATA, NVEC,
    EPSREL, EPSABS, VERBOSE, SEED,
    MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
    GRIDNO, STATEFILE, SPIN,
    &neval, &fail, integral, error, prob);

  printf("VEGAS RESULT:\tneval %d\tfail %d\n",
    neval, fail);
  for( comp = 0; comp < NCOMP; ++comp )
    printf("VEGAS RESULT:\t%.8f +- %.8f\tp = %.3f\n",
      (double)integral[comp], (double)error[comp], (double)prob[comp]);
#endif

#if 1
  printf("\n-------------------- Suave test --------------------\n");

  Suave(NDIM, NCOMP, Integrand, USERDATA, NVEC,
    EPSREL, EPSABS, VERBOSE | LAST, SEED,
    MINEVAL, MAXEVAL, NNEW, NMIN, FLATNESS,
    STATEFILE, SPIN,
    &nregions, &neval, &fail, integral, error, prob);

  printf("SUAVE RESULT:\tnregions %d\tneval %d\tfail %d\n",
    nregions, neval, fail);
  for( comp = 0; comp < NCOMP; ++comp )
    printf("SUAVE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
      (double)integral[comp], (double)error[comp], (double)prob[comp]);
#endif

#if 1
  printf("\n------------------- Divonne test -------------------\n");

  Divonne(NDIM, NCOMP, Integrand, USERDATA, NVEC,
    EPSREL, EPSABS, VERBOSE, SEED,
    MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
    BORDER, MAXCHISQ, MINDEVIATION,
    NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
    STATEFILE, SPIN,
    &nregions, &neval, &fail, integral, error, prob);

  printf("DIVONNE RESULT:\tnregions %d\tneval %d\tfail %d\n",
    nregions, neval, fail);
  for( comp = 0; comp < NCOMP; ++comp )
    printf("DIVONNE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
      (double)integral[comp], (double)error[comp], (double)prob[comp]);
#endif

#if 1
  printf("\n-------------------- Cuhre test --------------------\n");

  Cuhre(NDIM, NCOMP, Integrand, USERDATA, NVEC,
    EPSREL, EPSABS, VERBOSE | LAST,
    MINEVAL, MAXEVAL, KEY,
    STATEFILE, SPIN,
    &nregions, &neval, &fail, integral, error, prob);

  printf("CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n",
    nregions, neval, fail);
  for( comp = 0; comp < NCOMP; ++comp )
    printf("CUHRE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
      (double)integral[comp], (double)error[comp], (double)prob[comp]);
#endif

	///////////////////
	// --- VEGAS --- //
	///////////////////
	
	// Input parameters
	int ndim = 3;
	const int ndmx = 50;
	const int mxdim = 10;

	double *region = new double[2*ndim];
	region[0] =  0.0;
	region[1] = +1.0;
	region[2] =  0.0;
	region[3] = +1.0;
	region[4] =  0.0;
	region[5] = +1.0;

	int    init  = -1;
	int    ncall = 5000;
	int    itmx = 10;
	int    nprn  = 1;
	double tgral = 0;
	double sd = 0;
	double acc = 0.01;

	// output
	double chi2a ;

	double **xi = new double*[mxdim];
	for (int i= 0; i<mxdim; i++)
	{
		xi[i] = new double[ndmx];
	}

	int it;
	int ndo ;
	double si;
	double swgt;
	double schi;


 vegas_(region, &ndim, &fxn, &init, &ncall,
			 &itmx, &nprn, &tgral, &sd, &chi2a, &acc,
			 xi, &it, &ndo, &si, &swgt, &schi );

  return 0;

}

