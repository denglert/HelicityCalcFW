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
#define x1 xx[0]
#define x2 xx[1]
#define x3 xx[2]
#define x4 xx[3]
#define x5 xx[4]
#define f ff[0]

//const double  M  = 1.77682;
//const double  M  = 0.1056583715;
const double  M  = 1.0;
//const double m1  = 0.000510998928;
const double m1  = 0.00000;
const double m2  = 0.00000;
const double m3  = 0.00000;

const double BR  = 1.00000;

const double gamma_decay_PDG = hbar_c / c_tau_muon;
const double gamma_decay_analytic = pow(M,5.0) * G_Fermi * G_Fermi / 192.0 / pow(M_PI,3.0);
const double gamma_decay_analytic_reduced = pow(M,5.0) / 192.;

static int Integrand(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata)
{
	

   ThreeBodyDecay muon(M, m1, m2, m3);

	TLorentzVector *P = muon.P;
	TLorentzVector *p1 = muon.p[0];
	TLorentzVector *p2 = muon.p[1];
	TLorentzVector *p3 = muon.p[2];
	double weight = muon.GetPhaseSpaceWeight(x1,x2,x3,x4,x5);
  	double amp = M*(p3->E()) * p2->Dot( (*p1) );

  	//f = amp*weight;
  	f = weight;
  	return 0;

}


////////////////////////////

int main() 
{

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
#endif

//  double x1_ = 1.00;
//  double x2_ = 0.0;
//  double x3_ = 0.0;
//  double x4_ = 1.0;
//  double x5_ = 0.75;

  std::cout << "(double)integral[comp]: " << (double)integral[comp] << std::endl;
  std::cout << "Gamma(PDG): " << gamma_decay_PDG << std::endl;

  ThreeBodyDecay muon(M, m1, m2, m3);
  //double weight = muon.GetPhaseSpaceWeight(x1_,x2_,x3_,x4_,x5_);
//  muon.DisplayAll();
  
  double PSConst = muon.GetPSConst();
  double unit = PSConst * (M-m1-m2-m3) * pow(2*M_PI,3) * 16; 

  int i = 0;
  int j = 0;

  double result = (double)integral[comp];
  double result_reduced;
//
//  result = 32 * result * pow(2,i) * pow(M_PI,j) * BR  ;
//  result_reduced = (M-m1-m2-m3)*result/128.;
//
//  result = result * PSConst * G_Fermi * G_Fermi;
  result = result * PSConst * pow(4*M_PI,3) * 4;

  double ratio_reduced = result_reduced/gamma_decay_analytic_reduced;
  double ratio = result/gamma_decay_analytic;

// const double gamma_decay_analytic = pow(M,5.0) * G_Fermi * G_Fermi / 192.0 / pow(M_PI,3.0);
// const double gamma_decay_analytic_reduced = pow(M,5.0) / 192.0;

  std::cout << "i: " << i << std::endl;
  std::cout << "j: " << j << std::endl;
  std::cout << "Gamma(PDG): " << gamma_decay_PDG << std::endl;
  std::cout << "Gamma(analytic): " << gamma_decay_analytic << std::endl;
  std::cout << "Gamma(our result): " << result << std::endl;
  std::cout << "PSConst * (M-m1-m2-m3) * pow(2*M_PI,3) * 16 = " << unit << std::endl;
  std::cout << "ratio: (our result)/gamma_decay_analytic: " << ratio << std::endl;
  std::cout << "ratio_reduced: (our result)/gamma_decay_analytic: " << ratio_reduced << std::endl;
  std::cout << "1./ratio_reduced: " << 1./ratio_reduced << std::endl;


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

