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

const double muon_p     = 0.065;
const double muon_theta = 0.234*M_PI;
const double muon_phi   = 0.923*M_PI;

const double gamma_PDG = hbar_c / c_tau_muon;
const double gamma_formula = pow(M,5.0) * G_Fermi * G_Fermi / 192.0 / pow(M_PI,3.0);

// Flags and bits
const int k0_flag           = 1; 	// 1 - custom, 2 - physical
const int polvec_flag       = 2; 	// 1 - custom, 2 - physical
const bool bit_addpolarized = true;

static int Integrand(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata)
{
	

   ThreeBodyDecay muon(M, m1, m2, m3);

	muon.SetMotherMPThetaPhi(M,muon_p,muon_theta,muon_phi);
	muon.SetBitBoostBack(true);
	
	TLorentzVector *P  = muon.P;
	TLorentzVector *p1 = muon.p[0];
	TLorentzVector *p2 = muon.p[1];
	TLorentzVector *p3 = muon.p[2];

	// Auxiliary vector(s)
	TLorentzVector k0;

	TLorentzVector k0_custom (0.0, 0.0, 1.0, 1.0);
	TLorentzVector k0_physical;
	
	double Energy = P->E();
	double pmagnitude = P->P();

	k0_physical[3] = 1;
	for (int i = 0; i<3; i++)
	{ k0_physical[i] = (*P)[i]/pmagnitude; }

	for (int nu = 0; nu<4; nu++)
	{ k0_physical[nu] = k0_physical[nu]/(Energy+pmagnitude); }
	
	if ( k0_flag == 1)
	{ k0 = k0_custom; }

	if ( k0_flag == 2)
	{ k0 = k0_physical; }

	double Pk0 = P->Dot(k0);

	// Spin polarization vector
	TLorentzVector polvec;

	// Custom spin polarization vector with p and k0
	if (polvec_flag == 1)
	{
		for(int nu = 0; nu<4; nu++)
		{ polvec[nu] = ( (*P)[nu]/M) - (M/Pk0)*k0[nu]; }
	}

	// Helicity spin polarization vector
	// Note:
	// This should be equivalent to the custom pol. vector
	// when choosing physical k0
	if (polvec_flag == 2)
	{
		polvec[3] = pmagnitude*pmagnitude;

		for (int i = 0; i<3; i++)
		{
			polvec[i] = Energy*(*P)[i];
		}

		for(int nu = 0; nu<4; nu++)
		{ polvec[nu] = polvec[nu]/M/pmagnitude; }
	}

	///////////////////////////////////////////////
	double weight = muon.GetPhaseSpaceWeight(x1,x2,x3,x4,x5);
  	double amp_unpolarized = 128*P->Dot( (*p3) )  * p1->Dot( (*p2) );
  	double amp_polarized   = M*(polvec*(*p3)) * ( (*p1) * (*p2) );
	//double amp = amp_unpolarized + amp_polarized;
	
	double amp;

	if ( bit_addpolarized == false)
	{
		amp = amp_unpolarized;
	}
	
	if ( bit_addpolarized == true)
	{
		amp = amp_unpolarized + amp_polarized;
	}

  	f = amp*weight;
  	//f = weight;
  	return 0;

};


////////////////////////////

int main() 
{

  int comp, nregions, neval, fail;
  cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];

  ThreeBodyDecay muon(M, m1, m2, m3);
  muon.SetMotherMPThetaPhi(M,muon_p,muon_theta,muon_phi);
  TLorentzVector *P = muon.P;

  double PSConst = muon.GetPSConst();

#if 1
  printf("###############################################\n");
  printf("### --- Muon decay lifetime calculation --- ###\n");
  printf("###############################################\n");
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

  result = result * PSConst * G_Fermi * G_Fermi / P->E() / 2.0 / 2.0 ;

  double ratio_formula = result/gamma_formula;
  double ratio_PDG = result/gamma_PDG;

  double  tau = hbar/result;
  double ctau = tau*c;

  printf("\n\n");

  printf("--------------------------\n");
  printf("--- Physical constants ---\n");
  printf("--------------------------\n");

  printf("\n");
  printf("hbar:    %12.6e [GeV s]\n", hbar);
  printf("c:       %12.6f [m/s]\n", G_Fermi);
  printf("hbar*c:  %12.6f [GeV fm]\n", hbar_c);
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
  printf("  p     --->      q        k1           k2 \n");

  printf("\n");
  printf("Initial Muon Configuration:\n");
  printf("|pvec|: %12.6f [GeV/c]\n", muon_p);
  printf("theta:  %12.6f \n", muon_theta);
  printf("phi:    %12.6f \n", muon_phi);
  printf("gamma:  %12.6f \n", P->Gamma());
  printf("beta:   %12.6f \n", P->Beta());
  printf("Four-momentum (p):\n");
  displayTLorentzVector(P);
  printf("\n");

	TLorentzVector k0_custom (0.0, 0.0, 1.0, 1.0);
	TLorentzVector k0_physical;
	TLorentzVector k0;
	
	double Energy = P->E();
	double pmagnitude = P->P();

	k0_physical[3] = 1;
	for (int i = 0; i<3; i++)
	{ k0_physical[i] = (*P)[i]/pmagnitude; }

	for (int nu = 0; nu<4; nu++)
	{ k0_physical[nu] = k0_physical[nu]/(Energy+pmagnitude); }

  if ( k0_flag == 1)
  { k0 = k0_custom; }
  
  if ( k0_flag == 2)
  { k0 = k0_physical; }

  printf("Auxiliary vector k0\n");
  displayTLorentzVector(&k0);
  printf("\n");

	// Spin polarization vector
	TLorentzVector polvec;
	double Pk0 = P->Dot(k0);

	// Custom spin polarization vector with p and k0
	if (polvec_flag == 1)
	{
		for(int nu = 0; nu<4; nu++)
		{ polvec[nu] = ( (*P)[nu]/M) - (M/Pk0)*k0[nu]; }
	}

	// Helicity spin polarization vector
	// Note:
	// This should be equivalent to the custom pol. vector
	// when choosing physical k0
	if (polvec_flag == 2)
	{
		polvec[3] = pmagnitude*pmagnitude;

		for (int i = 0; i<3; i++)
		{
			polvec[i] = Energy*(*P)[i];
		}

		for(int nu = 0; nu<4; nu++)
		{ polvec[nu] = polvec[nu]/M/pmagnitude; }
	}

  printf("Spin polarization vector\n");
  if ( polvec_flag == 1 )
  {
  printf("s^{mu} = P^{mu}/M - (m/Pk0)*k0^{mu}\n");
  }
  if ( polvec_flag == 2 )
  {
  printf("s^{mu} = ( |pvec|^2 , p0 pvec) / (m |pvec|))\n");
  }
  displayTLorentzVector(&polvec);

  printf("\n");
  printf("Unpolarized amplitude:\n");
  printf("(general form)\n");
  printf("128*G_Fermi^{2}*(p k2)*(q k1)\n");
  printf("(in the rest frame of muon)\n");
  printf("128*G_Fermi^{2}*M*E2*(q k1)\n");

  printf("Polarized amplitude:\n");
  printf("+/- M*(s k2) (q k1)\n");

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
  printf("tau(our result): hbar/Gamma(our result)\n");
  
  printf("\n");
  printf("-------------------------\n");
  printf("--- Numerical results ---\n");
  printf("-------------------------\n");

  printf("Note: our result are quoted in the LAB frame, while the PDG\n");
  printf("      and the formula are calculated in the muon rest frame!\n");
  printf("\n");

  printf("Important flags:\n", bit_addpolarized);
  printf("k0_type:                    %d (1 = custom (see above), 2 = physical)\n", k0_flag);
  printf("spin_polarization:          %d (1 = computed with k0,   2 = helicity)\n", polvec_flag);
  printf("Polarized component added?: %d (0 = no, 1 = yes)\n", bit_addpolarized);
  printf("\n");

  printf("Gamma(PDG):        %12.6e [GeV] (note: this is the total gamma!)\n", gamma_PDG);
  printf("Gamma(formula):    %12.6e [GeV]\n", gamma_formula);
  printf("Gamma(our result): %12.6e [GeV]\n", result);

  printf("\n");
  printf(" tau(our result):  %12.6e s\n", tau);
  printf("ctau(our result):  %12.6e m\n", ctau);
  printf("ctau(PDG):         %12.6e m\n", c_tau_muon*1e-15);
  printf("\n");
  printf("Muon\n");
  printf("beta:                            %12.6f\n", P->Beta());
  printf("gamma:                           %12.6f\n", P->Gamma());
  printf("ctau ratio(our result)/rest_PDG: %12.6f\n", ctau/c_tau_muon/1e-15);
  printf("!!! COMPARE THE ABOVE !!! gamma vs. ctau ratio\n");

  //printf("\n");
  //printf("(our result)/Gamma(formula): %12.6f \n", ratio_formula);
  //printf("(our result)/Gamma(PDG):     %12.6f \n", ratio_PDG);


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

