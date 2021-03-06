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
const double  M  = 0.1056583715;
//const double  M  = 1.0;
//const double m1  = 0.000510998928;
const double m1  = 0.00000;
const double m2  = 0.00000;
const double m3  = 0.00000;

const double BR  = 1.00000;

const double gamma_PDG = hbar_c / c_tau_muon;
const double gamma_formula = pow(M,5.0) * G_Fermi * G_Fermi / 192.0 / pow(M_PI,3.0);

static int Integrand(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata)
{
	

   ThreeBodyDecay muon(M, m1, m2, m3);

	TLorentzVector *P = muon.P;
	TLorentzVector *p1 = muon.p[0];
	TLorentzVector *p2 = muon.p[1];
	TLorentzVector *p3 = muon.p[2];
	double weight = muon.GetPhaseSpaceWeight(x1,x2,x3,x4,x5);
  	double amp = 32*(p3->E()) * p1->Dot( (*p2) );

  	f = amp*weight;
  	//f = weight;
  	return 0;

}


////////////////////////////

int main() 
{

  int comp, nregions, neval, fail;
  cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];

  ThreeBodyDecay muon(M, m1, m2, m3);
  double PSConst = muon.GetPSConst();

#if 1
  printf("####################################################\n");
  printf("### --- test_PhaseSpaceIntegration               ###\n");
  printf("####################################################\n");
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


  double result = (double)integral[comp];

  result = result * PSConst * G_Fermi * G_Fermi;

  double ratio_formula = result/gamma_formula;
  double ratio_PDG = result/gamma_PDG;

  printf("\n\n");
  printf("--------------------------\n");
  printf("--- Physical constants ---\n");
  printf("--------------------------\n");

  printf("\n");
  printf("hbarc:   %12.6f [GeV fm]\n", hbar_c);
  printf("G_Fermi: %12.6f [GeV^{-2}]\n", G_Fermi);

  printf("\n");
  printf("Muon constants:\n");
  printf("m (muon):                  %12.6f [GeV]\n", M);
  printf("c_tau:                     %12.6e  [fm]\n", c_tau_muon);
  printf("Gamma(PDG) = hbarc/c_tau = %12.6e [GeV]\n", gamma_PDG);


  printf("\n");
  printf("Other masses:\n");
  printf("m (elec):       %12.6f [GeV]\n", m1);
  printf("m (nu_e):       %12.6f [GeV]\n", m2);
  printf("m (nu_m):       %12.6f [GeV]\n", m3);

  printf("\n");
  printf("-------------------------\n");
  printf("---   Decay process   ---\n");
  printf("-------------------------\n");

  printf("\n");
  printf("(muon)- ---> (electron)- (nu_mu) (nu_electronbar)\n");

  printf("\n");
  printf("Amplitude:\n");
  printf("(general form)\n");
  printf("128*G_Fermi^{2}*(p k2)*(q k1)\n");
  printf("(in the rest frame of muon)\n");
  printf("128*G_Fermi^{2}*M*E2*(q k1)\n");

  printf("\n");
  printf("Other factors:\n");
  printf("- From normalization:\n");
  printf("  1.0/(2*E) = 1.0/2*M (in the rest frame)\n");

  printf("- Spin averaging for the muon:\n");
  printf("  1.0/2.0\n");

  printf("\n");
  printf("ThreeBodyDecay class\n");
  printf("PSConstant (formula)  s23_length/pow(M_PI,3.0)/128.0\n");
  printf("PSConstant  (numval): %12.6e \n", PSConst);

  printf("\n");
  printf("Gamma(formula): G_Fermi^{2}*m_mu^{5}/(192*pi^{3})\n");
  
  printf("\n");
  printf("-------------------------\n");
  printf("--- Numerical results ---\n");
  printf("-------------------------\n");

  printf("Gamma(PDG):        %12.6e [GeV] (note: this is the total gamma!)\n", gamma_PDG);
  printf("Gamma(formula):    %12.6e [GeV]\n", gamma_formula);
  printf("Gamma(our result): %12.6e [GeV]\n", result);

  printf("ratio:\n");
  printf("(our result)/gamma_formula: %12.6f \n", ratio_formula);
  printf("(our result)/gamma_PDG:     %12.6f \n", ratio_PDG);


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

