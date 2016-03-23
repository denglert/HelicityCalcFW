#include <iostream>
#include <cmath>
#include <TH3D.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include "HelicityFW.h"
#include "FortranInterface.h"
#include "PhysConst.h"
#include "UtilFunctions.h"
#include "PhaseSpaceTools.h"
#include "cuba.h"

//////////////////////////////
// -- Cuba configuration -- //
//////////////////////////////

#define USERDATA NULL
#define NCOMP 1
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

#define NDIM 5
#define xA1 xx[0]
#define xA2 xx[1]
#define xA3 xx[2]
#define xA4 xx[3]
#define xA5 xx[4]
#define f ff[0]


// -- Masses -- //
const double M   = 100.00;
const double mA  = 5.00;
const double mB  = 5.00;
const double mA1 = 0.00;
const double mA2 = 0.00;
const double mA3 = 0.00;
const double mB1 = 0.00;
const double mB2 = 0.00;
const double mB3 = 0.00;

// -- Mother particle momentum -- //
const double Px = 0.00;
const double Py = 0.00;
const double Pz = 0.00;
const double E  = sqrt( Px*Px + Py*Py + Pz*Pz + M*M );

const double xAB1 = 0.50;
const double xAB2 = 0.50;
const double  xB1 = 0.50;
const double  xB2 = 0.50;
const double  xB3 = 0.50;
const double  xB4 = 0.50;
const double  xB5 = 0.50;

// -- BoostBack flag -- //
const bool BitBoostBack = true;

DecayChain126 decay(M, mA,  mB,
					 	     mA1, mA2, mA3,
							  mB1, mB2, mB3);


// -- Defining pointers to momenta
TLorentzVector *P; 
TLorentzVector *pA;
TLorentzVector *pB;
TLorentzVector *pA1;
TLorentzVector *pA2;
TLorentzVector *pA3;
TLorentzVector *pB1;
TLorentzVector *pB2;
TLorentzVector *pB3;


static int Integrand(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata)
{
	

  decay.SetPhaseSpace( xAB1, xAB2,
						 	   xA1, xA2, xA3, xA4, xA5,
							   xB1, xB2, xB3, xB4, xB5);

   double PSWeight_A123 = decay.DecayA123->GetPhaseSpaceWeight(xA1, xA2, xA3, xA4, xA5);

	decay.DisplayAllInfo();

  	f = PSWeight_A123;
  	return 0;

}


////////////////////////////

int main() 
{

  printf("\n");
  printf("##########################################################\n");
  printf("### --- Phi -> 2f -> 6f, single branch integration --- ###\n");
  printf("##########################################################\n");
  printf("\n");


  P = decay.P;
  P->SetPxPyPzE(Px,Py,Pz,E);
  decay.SetBitBoostBack( BitBoostBack );
  printf("\n");
  decay.SetTag("DecayChain126");

  
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
  comp = 0;
    printf("VEGAS RESULT:\t%.8f +- %.8f\tp = %.3f\n",
      (double)integral[comp], (double)error[comp], (double)prob[comp]);

	 double num_result = (double)integral[comp];
	 printf("VEGAS numerical result: %10.5f\n", num_result);


  // Get phase space constant from the DecayChain126 class
  double PSConst = decay.DecayA123->GetPSConst();
  // Final result 
  double result = num_result * PSConst;

  // Expected value in the massless scenario:
  double threebodyA123 = mA*mA/(256.0*M_PI*M_PI*M_PI);
  double expected_result = threebodyA123 ;

  printf("\n");
  printf("                                           final result: %10.10e\n", result);
  printf("Expexted result in massless (m1 = m2 = m3 = 0) scenario: %10.10e\n", expected_result);
  printf("                                            their ratio: %10.10f\n", result/expected_result);
  printf("\n");

  double PSWeight_B123 = decay.DecayB123->GetPhaseSpaceWeight(xB1, xB2, xB3, xB4, xB5);
  printf("PSWeight_B123 = %10.5e \n", PSWeight_B123);
	
#endif
#if 0
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

#if 0
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

#if 0
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

  return 0;

}

