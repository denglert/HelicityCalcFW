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
#define MAXEVAL 100000

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

////////////////////////
// -- Input values -- //
////////////////////////

// -- Physical constansts and computation
//const double  M  = 1.77682;
const double  M  = 0.1056583715;
//const double  M  = 10.0;
//const double  M  = 1.0;
//const double m1  = 0.000510998928;
const double m1  = 0.00000;
const double m2  = 0.00000;
const double m3  = 0.00000;

const double BR  = 1.00000;

//const double muon_p     = 0.05;
const double muon_p     = 0.05;
const double muon_theta = 0.234*M_PI;
const double muon_phi   = 0.923*M_PI;

//const double phat_theta = muon_theta;
//const double phat_phi   = muon_phi;

const double phat_theta = 0.0*M_PI;
const double phat_phi   = 0.3*M_PI;

//const double phat34_theta = muon_theta;
//const double phat34_phi   = muon_phi;

const double phat34_theta =  0.555*M_PI;
const double phat34_phi   = -0.412*M_PI;

const double gamma_PDG = hbar_c / c_tau_muon;
const double gamma_formula = pow(M,5.0) * G_Fermi * G_Fermi / 192.0 / pow(M_PI,3.0);
const double ctau_formula = hbar_c/gamma_formula;

// -- Flags and bits
const int k0_flag        = 2; 	// 1 - custom, 2 - physical
const int polvec_flag    = 2; 	// 1 - custom, 2 - physical
const int ampltiude      = 34;    // 0 - unpolarized, +1 polarized, -1 polarized, 3,4,34
const bool BitBoostBack  = true;

//////////////////////////////////////
// -- Declaring global variables -- //
//////////////////////////////////////
// -- main() and integrand() share common variables
// -- and functions that are accessible from both

// Declaring global ThreeBodyDecay class
ThreeBodyDecay muon(M, m1, m2, m3);

// Auxiliary vector
TLorentzVector k0;           // Final auxiliary vector used
TLorentzVector k0_custom;    // Arbitrary auxiliary vector
TLorentzVector k0_physical;  // Physical auxiliary vector

// Pointers to the momenta
TLorentzVector *P;   // Muon    p
TLorentzVector *p1;  // e-      q
TLorentzVector *p2;  // v_mu    k
TLorentzVector *p3;  // v_ebar  k'

// Spin ampltiude vector
TLorentzVector polvec;
TLorentzVector polvec34;
TVector3 phat;
TVector3 phat34;
 
/////////////////////
// -- Integrand -- //
/////////////////////

static int Integrand(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata)
{

//	TVector3 Boostvec = P->BoostVector();
//	p3->Boost(-Boostvec);
   
//	p1->Boost(-Boostvec);
//	p2->Boost(-Boostvec);
	
	double amp;
	double weight = muon.GetPhaseSpaceWeight(x1,x2,x3,x4,x5);
	
	// -- Default case -- //
	double amp_unpolarized = 128*P->Dot( (*p3) )  * p1->Dot( (*p2) );
	double amp_polarized1  =  64*P->Dot( (*p3) )  * p1->Dot( (*p2) ) - 64*M*(polvec*(*p3)) * ( (*p1) * (*p2) );
	double amp_polarized2  =  64*P->Dot( (*p3) )  * p1->Dot( (*p2) ) + 64*M*(polvec*(*p3)) * ( (*p1) * (*p2) );
	
	// -- Custom case -- //
	double lambda3 =  1.0;
	double amp_polarized3 =         64*P->Gamma()*(1-lambda3*P->Beta()) * M * p3->E()             * p1->Dot((*p2))-
			  						lambda3*64*P->Gamma()*(1-lambda3*P->Beta()) * M * polvec34.Dot((*p3)) * p1->Dot((*p2));
	
	double lambda4 = -1.0;
	double amp_polarized4 =         64*P->Gamma()*(1-lambda4*P->Beta()) * M * p3->E()             * p1->Dot((*p2))-
			  						lambda4*64*P->Gamma()*(1-lambda4*P->Beta()) * M * polvec34.Dot((*p3)) * p1->Dot((*p2));
	
	double amp_unpolarized34 = amp_polarized3 + amp_polarized4;
	
	// -- Choosing amplitude -- //
	switch(ampltiude)
	{ 
		// Default
		case  0: amp = amp_unpolarized; break;
		case  1: amp = amp_polarized1;  break;
		case -1: amp = amp_polarized2;  break;
	
		// Custom
		case   3: amp = amp_polarized3;     break;
		case   4: amp = amp_polarized4;     break;
		case  34: amp = amp_unpolarized34;  break;
	}

	f = amp*weight;
	//f = weight;
	
  	return 0;
};


////////////////////////
// -- Main program -- //
////////////////////////

int main() 
{

	// -- Declaring Cuba variables -- //
	int comp, nregions, neval, fail;
	cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];
	
	// -- ThreeBodyDecay initialization-- // 
	muon.SetMotherMPThetaPhi(M,muon_p,muon_theta,muon_phi);
	muon.SetBitBoostBack(BitBoostBack);

	// -- Set up pointers to the momenta
	P  = muon.P;
	p1 = muon.p[0];
	p2 = muon.p[1];
	p3 = muon.p[2];

	// -- Get muon kinematics
   double P_gamma  = P->Gamma();
   double P_beta   = P->Beta();
   double P_M      = P->M();
   double P_MomMag = P->P();

	// -- Polarization vector -- //
	// -- Custom vector construction
   k0_custom = TLorentzVector(1.0, 0.0, 0.0, 1.0);
   // k0_custom = TLorentzVector(3, 2, 0.0, 1.0);
   // k0_custom = TLorentzVector(1.0, 0.0, 0.0, 1.0);
	double alpha = 1;

	// -- Phyicsal vector construction
	k0_physical[3] = 1;
	for (int i = 0; i<3; i++)
	{ k0_physical[i] = (*P)[i]/P_MomMag; }

	for (int nu = 0; nu<4; nu++)
	{ k0_physical[nu] = alpha*k0_physical[nu]/(P_M*P_gamma*(1+P_beta)); }
	
	// -- Choosing k0 auxiliary vector
	if ( k0_flag == 1)
	{ k0 = k0_custom; }

	if ( k0_flag == 2)
	{ k0 = k0_physical; }

	// Get pk0 scalar product
	double Pk0 = P->Dot(k0);

	// Custom spin polarization vector with p and k0
	if (polvec_flag == 1)
	{
		for(int nu = 0; nu<4; nu++)
		{ polvec[nu] = ( (*P)[nu]/M) - (M/Pk0)*k0[nu]; }
	}

	// Helicity spin ampltiude vector
	// Note:
	// This should be equivalent to the custom pol. vector
	// when choosing physical k0
	// Default:
	if (polvec_flag == 2)
	{

		TVector3 phat(1.0,0.0,0.0);
		phat.SetPhi(muon_phi);
		phat.SetTheta(muon_theta);

		polvec = TLorentzVector(phat,P->Beta());

		for(int nu = 0; nu<4; nu++)
		{ polvec[nu] = polvec[nu]*P_gamma; }
	}

	phat34 = TVector3(1.0,0.0,0.0);
   phat34.SetPhi(phat34_phi);
   phat34.SetTheta(phat34_theta);
	polvec34 = TLorentzVector(phat34,0);



/////////////////////////////////////////////////////////////////////////

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

#if 0
  printf("\n-------------------- Suave test --------------------\n");

  Suave(NDIM, NCOMP, Integrand, USERDATA, NVEC,
    EPSREL, EPSABS, VERBOSE | LAST, SEED,
    MINEVAL, MAXEVAL, NNEW, NMIN, FLATNESS,
    STATEFILE, SPIN,
    &nregions, &neval, &fail, integral, error, prob);

  printf("SUAVE RESULT:\tnregions %d\tneval %d\tfail %d\n",
    nregions, neval, fail);
  comp = 0;
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

  double result = (double)integral[comp];

  // If result is unpolarized divide by 2.0
  if ( (ampltiude == 0) || (ampltiude == 34) )
  { result = result / 2.0;}

  double PSConst = muon.GetPSConst();
  double result_wogamma = result * PSConst * G_Fermi * G_Fermi / M / 2.0  ;
  result = result * PSConst * G_Fermi * G_Fermi / P->E() / 2.0 ;

  double ratio_formula = result/gamma_formula;
  double ratio_PDG = result/gamma_PDG;

  double  tau = hbar/result;
  double ctau = tau*c;

  double  tau_wogamma = hbar/result_wogamma;
  double ctau_wogamma = tau_wogamma*c;


  // -- Displaying info -- //

  printf("\n");
  printf("###############################################\n");
  printf("### --- Muon decay lifetime calculation --- ###\n");
  printf("###############################################\n");

  printf("\n");

  printf("--------------------------\n");
  printf("--- Physical constants ---\n");
  printf("--------------------------\n");

  printf("\n");
  printf("hbar:    %12.6e [GeV s]\n", hbar);
  printf("c:       %12.6e [m/s]\n", c);
  printf("hbar*c:  %12.6e [GeV fm]\n", hbar_c);
  printf("G_Fermi: %12.6e [GeV^{-2}]\n", G_Fermi);

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
  printf("  p     --->      q        k            k' \n");
  printf("  p     --->      q        k1           k2 \n");
  printf("  P     --->     p1        p2           p3 \n");
  printf("\n");

  printf("-------------------------\n");
  printf("---   Configuration   ---\n");
  printf("-------------------------\n");

  printf("\n");
  printf("############\n");
  printf("### Muon ###\n");
  printf("############\n");
  printf("\n");
  printf("|mom|:  %12.6f [GeV/c]\n", muon_p);
  printf("theta:  %12.6f [rad]\n", muon_theta);
  printf("phi:    %12.6f [rad]\n", muon_phi);
  printf("gamma:  %12.6f \n", P->Gamma());
  printf("beta:   %12.6f \n", P->Beta());
  printf("\n");
  printf("Momentum:\n");
  displayTLorentzVector(P);
  printf("\n");

  printf("Unitvector pointing in the direction of momentum:\n");
  printf("x: %12.6f\n", (*P)[0]/P->P());
  printf("y: %12.6f\n", (*P)[1]/P->P());
  printf("z: %12.6f\n", (*P)[1]/P->P());
  printf("\n");


  printf("################################\n");
  printf("### Spin Polarization vector ###\n");
  printf("################################\n");
  printf("### - Denoted by s^{mu}         \n");
  printf("\n");


  if ( (ampltiude == -1) || (ampltiude == 0) || (ampltiude == 1) )
  {

  if ( polvec_flag == 1 )
  {

   printf("polvec_flag: %d\n", polvec_flag);
  	printf("Spin polarization defined with auxiliary vector k0:\n");
  	printf("s^{mu} = P^{mu}/M - (m/Pk0)*k0^{mu}\n");
   printf("\n");

  	printf("k0_flag: %d\n", k0_flag);
   if ( k0_flag == 1 )
	{ printf("Custom k0:\n"); }

   if ( k0_flag == 2 )
	{ 
		printf("Choice of k0, yielding helicity polarization vector.\n");
		printf("Physical k0:\n");
	}

  	displayTLorentzVector(&k0);
  	printf("\n");

  }

  if ( polvec_flag == 2 )
  {
   printf("polvec_flag: %d\n", polvec_flag);
  	printf("Spin polarization vector = helicity polarization vector\n");
 	printf("s^{mu} = ( |pvec|^2 , p0 pvec) / (m |pvec|)) = \n");
 	printf("       = gamma (beta,phat vector)\n");
   printf("\n");
  }

  printf("Spin polarization vector components:\n");
  displayTLorentzVector(&polvec);

  }

  if ( (ampltiude == 3) || (ampltiude == 4) || (ampltiude == 34) )
  {
  	 printf("Spin polarization vector defined in rest frame of the muon\n");
  	 printf("s^{mu} = (0,phat34)\n");

   printf("\n");
  	 printf("where the phat34 is pointing towards:\n");
  	 printf("phat34 theta: %12.6f\n", phat34_theta);
  	 printf("phat34 phi:   %12.6f\n", phat34_phi);

   printf("\n");
  	 printf("phat34 components:\n");
  	 printf("x: %12.6f\n", phat34[0]);
  	 printf("y: %12.6f\n", phat34[1]);
  	 printf("z: %12.6f\n", phat34[2]);

  }

  printf("\n");
  printf("--------------------------\n");
  printf("--- Consistency checks ---\n");
  printf("--------------------------\n");

  printf("\n");
  printf("- Lorentz scalar product (ps) = 0 (orthogonality check)\n");
  printf("(ps): %12.6e     (should give really small value)\n", (*P)*polvec);

  printf("\n");
  printf("- k0^{2} = 0 (massless auxuliary vector)\n");
  printf("(k0)^2: %12.6e     (should give really small value)\n", k0.M2());


  printf("\n");
  printf("------------------\n");
  printf("--- Ampltidude ---\n");
  printf("------------------\n");

  printf("\n");
  printf("Amplitude formula chosen: --->  %d  <--- \n", ampltiude);
  printf("\n");

  printf("Amplitude formula list:\n");
  printf("(No):  ---------------------- Formula ---------------------------- (polarization)\n");
  printf("- Default case:                              \n");
  printf(" (0): 128 * (p k') * (q k)                                         (unpolarized)\n");
  printf("(-1):  64 * (p k') * (q k) + 64 * M * (s k') * (q k)               (polarized)\n");
  printf("(+1):  64 * (p k') * (q k) - 64 * M * (s k') * (q k)               (polarized)\n");

  printf("- Custom:                              \n");
  printf(" (+3): 64 * gamma * ( 1 - beta ) *     (M * (k')^{0}) * (q k) -    (polarized)\n");
  printf("      -64 * gamma * ( 1 - beta ) * M *    (s k')      * (q k)                 \n");
  printf(" (+4): 64 * gamma * ( 1 + beta ) *     (M * (k')^{0}) * (q k) +    (polarized)\n");
  printf("      +64 * gamma * ( 1 + beta ) * M *    (s k')      * (q k)                 \n");
  printf("(+34): sum of (+3) and (+4)                                        (unpolarized)\n");

  printf("\n");
  printf("Notations:\n");
  printf("(muon)- ---> (electron)- (nu_mu) (nu_electronbar)\n");
  printf("  p     --->      q        k            k' \n");
  printf("  p     --->      q        k1           k2 \n");
  printf("  P     --->     p1        p2           p3 \n");
  printf("\n");

  printf("Leftover factors multiplying the result:\n");
  printf("- Coupling constant\n");
  printf("  G_Fermi^{2}\n");

  printf("- Initial particle state normalization:\n");
  printf("  1.0/(2*E) = 1.0/2*M (in the rest frame)\n");

  printf("- Spin averaging for the muon (only if unpolarized):\n");
  printf("  1.0/2.0\n");

  printf("\n");
  printf("ThreeBodyDecay class\n");
  printf("PSConstant (formula)  s23_length/pow(M_PI,3.0)/128.0\n");
  printf("PSConstant  (numval): %12.6e \n", PSConst);

  printf("\n");
  printf("---------------------------\n");
  printf("--- Other formulas used ---\n");
  printf("---------------------------\n");

  printf("\n");
  printf("Textbook result of the integration:\n");
  printf("Gamma: G_Fermi^{2}*m_mu^{5}/(192*pi^{3})\n");

  printf("\n");
  printf("tau  = hbar/Gamma\n");
  printf("ctau = c*tau\n");
  
  printf("\n");
  printf("-------------------------\n");
  printf("--- Numerical results ---\n");
  printf("-------------------------\n");

  printf("Don't forget: our result are quoted in the LAB frame, while the PDG and\n");
  printf("              the textbook formula are calculated in the muon rest frame!\n");
  printf("\n");

  printf("# Important flags:\n", ampltiude);
  printf("- k0_type:                    %d (see above at 'Spin polarization', 1 = custom, 2 = physical)\n", k0_flag);
  printf("- spin_polarization:          %d (see above at 'Spin polarization', 1 = computed with k0, 2 = helicity)\n", polvec_flag);
  printf("- Amplitude formula chosen:   %d (see above at 'Amplitude' )\n", ampltiude);
  printf("\n");
  printf("Note: Choosing (k0_type = 2) and (spin_polarization = 1) should give the same result as\n");
  printf("      choosing (spin_polarization = 2) by definition\n");
  printf("\n");

  printf("Gamma(PDG):                             %12.6e [GeV] (note: this is the total gamma!)\n", gamma_PDG);
  printf("Gamma(textbook):                        %12.6e [GeV]\n", gamma_formula);
  printf("Gamma(our result w/o gamma factor):     %12.6e [GeV]\n", result_wogamma);
  printf("Gamma(our result):                      %12.6e [GeV]\n", result);

  printf("\n");
  printf("ctau(PDG):                              %12.6e [m]\n", c_tau_muon*1e-15);
  printf("ctau(textbook):                         %12.6e [m]\n", ctau_formula*1e-15);
  printf("ctau(our result w/o gamma factor):      %12.6e [m]\n", ctau_wogamma);
  printf("ctau(our result):                       %12.6e [m]\n", ctau);
  printf("\n");
  printf("# Muon configuration\n");
  printf("beta:                                                %12.6f\n", P->Beta());
  printf("gamma:                                               %12.6f <-\n", P->Gamma());

  printf("\n");
  printf("# Ratios\n");
  printf("ratio of ctau's: (our result)/(textbook rest frame)  %12.6f <-\n", ctau/ctau_formula/1e-15);
  printf("ratio of the above two values:                       %12.6f\n", P->Gamma()/(ctau/ctau_formula/1e-15));
  printf("\n");
  printf("Note: You should compare the values of 'gamma' with the 'ratio of ctau's',\n");
  printf("      as they should be equal.\n");
  printf("      The 'ratio of the above two values' gives a measure how well they agree.\n");
  printf("      It should be close to unity, but there could be some tiny deviation as\n");
  printf("      this is after all a numerical integration.\n");

  return 0;

}

