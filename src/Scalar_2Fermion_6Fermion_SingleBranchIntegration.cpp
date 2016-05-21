#include <iostream>
#include <cmath>
#include <TROOT.h>
#include <TStyle.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include "HelicityFW.h"
#include "FortranInterface.h"
#include "PhysConst.h"
#include "UtilFunctions.h"
#include "PhaseSpaceTools.h"
#include "HelicityTools.h"
#include "cuba.h"

//////////////////////////////
// -- Cuba configuration -- //
//////////////////////////////

#define USERDATA NULL
#define NCOMP 6
#define NVEC 1
#define EPSREL 1e-3
#define EPSABS 1e-12
#define VERBOSE 2
#define LAST 4
#define SEED 0
#define MINEVAL 0
#define MAXEVAL 50000
//#define MAXEVAL 50

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
#define f_combined1 ff[0]
#define f_combined2 ff[1]
#define f_h_tautau1 ff[2]
#define f_h_tautau2 ff[3]
#define f_taudecay1 ff[4]
#define f_taudecay2 ff[5]


// -- Masses -- //
const double M   = m_higgs;
//const double M   = 10.0;
const double mA  = m_tau;
const double mB  = m_tau;
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

//const double xAB1 = 0.50;
//const double xAB2 = 0.50;
//const double  xB1 = 0.50;
//const double  xB2 = 0.50;
//const double  xB3 = 0.50;
//const double  xB4 = 0.50;
//const double  xB5 = 0.50;
//
//double xAB1 = 0.39;
//double xAB2 = 0.81;
double xAB1 = 0.90;
double xAB2 = 0.50;

const double  xB1 = 0.73;
const double  xB2 = 0.12;
const double  xB3 = 0.71;
const double  xB4 = 0.23;
const double  xB5 = 0.12;

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

double pA_[4];
double pB_[4];
double pA1_[4];
double pA2_[4];
double pA3_[4];
double pB1_[4];
double pB2_[4];
double pB3_[4];

// Eucledian to Minkowski indeces mapping
int e2m[4] = {1, 2, 3, 0};


// -- Tau matrices
TauMatrix h_tautau;
TauMatrix taum;
TauMatrix taup;
TauMatrix c7_568;

// -- Tau matrices, normalized by (pik0)
// -- Matrix elements correspond to actual physical amplitudes
// -- Might be missing a factor of '8' or '16'
TauMatrix c_amp_res;
TauMatrix c_amp_htautau;
TauMatrix c_amp_dec_taum;
TauMatrix c_amp_dec_taup;
TauMatrix c_amp_7_568;

// Tau decay formula
// const double gamma_PDG = hbar_c / c_tau_muon;
const double Gamma_formula = pow(mA,5.0) * G_Fermi * G_Fermi / 192.0 / pow(M_PI,3.0);
const double ctau_formula = hbar_c/Gamma_formula;

// Integrated values
double integral_combined_pol_1;
double integral_combined_pol_2;
double integral_htautau_pol_1;
double integral_htautau_pol_2;
double integral_taudecay_pol_1;
double integral_taudecay_pol_2;


/////////////////////////
//// -- Integrand -- //// 
/////////////////////////

static int Integrand(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata)
{
	

  decay.SetPhaseSpace( xAB1, xAB2,
						 	   xA1, xA2, xA3, xA4, xA5,
							   xB1, xB2, xB3, xB4, xB5);


  // Convert to Minkowski metric
	for (int i = 0; i < 4; i++)
	{
		pA_[e2m[i]]  = (*pA)[i];
		pB_[e2m[i]]  = (*pB)[i];
		pA1_[e2m[i]] = (*pA1)[i];
		pA2_[e2m[i]] = (*pA2)[i];
		pA3_[e2m[i]] = (*pA3)[i];
		pB1_[e2m[i]] = (*pB1)[i];
		pB2_[e2m[i]] = (*pB2)[i];
		pB3_[e2m[i]] = (*pB3)[i];
	}


	rh_6f_(pA1_,pA2_,pB2_,pB1_,pA3_,pB3_);

	     c_amp_res.ReadInCMatrix_2_2( tau_amplitudes_.c_amp_res      );
	 c_amp_htautau.ReadInCMatrix_2_2( tau_amplitudes_.c_amp_htautau  );
	c_amp_dec_taum.ReadInCMatrix_2_2( tau_amplitudes_.c_amp_dec_taum );
	c_amp_dec_taup.ReadInCMatrix_2_2( tau_amplitudes_.c_amp_dec_taup );
	   c_amp_7_568.ReadInCMatrix_2_2( tau_amplitudes_.c_amp_7_568    );


   double PSWeight_A123 = decay.DecayA123->GetPhaseSpaceWeight(xA1, xA2, xA3, xA4, xA5);

//	c_amp_res.Show();
//	c_amp_htautau.Show();
//	c_amp_dec_taum.Show();

//	f_combined = PSWeight_A123;
//	f_h_tautau = PSWeight_A123;
//	f_taudecay = PSWeight_A123;

	double w = PSWeight_A123;

	f_combined1 = w * std::abs(  c_amp_7_568.m[1][0] ); // spinB is fixed, spinA is summed
	f_combined2 = w * std::abs(  c_amp_7_568.m[1][1] ); // spinB is fixed, spinA is summed
	f_h_tautau1 = std::abs(c_amp_htautau.m[0][0]);
	f_h_tautau2 = std::abs(c_amp_htautau.m[0][1]);   // spinB is fixed, spinA is summed
//	f_h_tautau = std::abs(c_amp_htautau.m[1][0]) + std::abs(c_amp_htautau.m[1][1]);   // spinB is fixed, spinA is summed
	f_taudecay1 = w * (std::abs(c_amp_dec_taum.m[1][0]));
	f_taudecay2 = w * (std::abs(c_amp_dec_taum.m[1][1])); // unpolarized amplitude for 

//	printf("c_amp_res: %.2e\n", std::abs(c_amp_res.m[1][1]));
//	printf("c_amp_res: %.2e\n", std::abs(c_amp_res.m[1][1]));
	
//	taum.Show();
//	decay.DisplayAllInfo();

  	return 0;

}

/////////////////////////////
//// -- Main Function -- //// 
/////////////////////////////
int main() 
{

  
  // Histogram for values
  //
  const double x1nBins = 0;
  const double x2nBins = 0;
  const double x1Min = 0.49;
  const double x1Max = 0.51;
  const double x2Min = 0.00;
  const double x2Max = 0.01;
//  const double x2Min = 0.49;
//  const double x2Max = 0.51;
  const double x1BinWidth = (x1Max-x1Min)/x1nBins;
  const double x2BinWidth = (x2Max-x2Min)/x2nBins;

  const double thetaMin = x1Min*M_PI;
  const double thetaMax = x1Max*M_PI;
  const double phiMin = x2Min*2.0*M_PI;
  const double phiMax = x2Max*2.0*M_PI;
  TH2D *hist_Gamma_Combined_full_pol_1 = new TH2D ("Gamma_Combined_full_pol_1",";#theta [rad];#phi [rad];value",x1nBins,thetaMin,thetaMax,x2nBins,phiMin,phiMax);
  TH2D *hist_Gamma_Combined_full_pol_2 = new TH2D ("Gamma_Combined_full_pol_2",";#theta [rad];#phi [rad];value",x1nBins,thetaMin,thetaMax,x2nBins,phiMin,phiMax);
  TH2D *hist_Gamma_Prod_BR_full = new TH2D ("Gamma_Prod_BR_full",";#theta [rad];#phi [rad];value",x1nBins,thetaMin,thetaMax,x2nBins,phiMin,phiMax);

    
//  xAB1 = 0.9;
//  xAB2 = 0.9;

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

  pA = decay.pA;
  pB = decay.pB;
  pA1 = decay.pA1;
  pA2 = decay.pA2;
  pA3 = decay.pA3;
  pB1 = decay.pB1;
  pB2 = decay.pB2;
  pB3 = decay.pB3;

  h_tautau.SetName("h_tautau");
  taum.SetName("taum");
  taup.SetName("taup");
  c7_568.SetName("c7_568");

  c_amp_htautau.SetName("c_amp_htautau");
  c_amp_dec_taum.SetName("c_amp_dec_taum");
  c_amp_dec_taup.SetName("c_amp_dec_taup");
  c_amp_7_568.SetName("c_amp_7_568");
  c_amp_res.SetName("c_amp_res");

	// Fortran COMMON blocks
	masses_.rmtau             = m_tau;
	couplings_.wcl            = 1.0;
	couplings_.gh_tautau      = 1.0;
	amplitudes_.rh_6f_tautau  = 0.0;
	amplitudes_.rh_6f_taum    = 0.0;
	amplitudes_.rh_6f_taup    = 0.0;
	amplitudes_.rh_6f_res     = 0.0;
	amplitudes_.rh_6f_res_nwa = 0.0;

  int comp, nregions, neval, fail;
  cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];

// Test suite
#if 1


  printf("\n");
  printf("###########################################################\n");
  printf("### --- Testing DecayChain126 phase space generator --- ###\n");
  printf("###########################################################\n");
  printf("\n");

  const double xAB1_ = xAB1;
  const double xAB2_ = xAB2;
  const double xA1_  = 0.50;
  const double xA2_  = 0.50;
  const double xA3_  = 0.50;
  const double xA4_  = 0.50;
  const double xA5_  = 0.50;
  const double xB1_  = 0.50;
  const double xB2_  = 0.50;
  const double xB3_  = 0.50;
  const double xB4_  = 0.50;
  const double xB5_  = 0.50;

  decay.SetPhaseSpace( xAB1_, xAB2_,
						 	   xA1_, xA2_, xA3_, xA4_, xA5_,
							   xB1_, xB2_, xB3_, xB4_, xB5_);


  printf("xAB1_ : %.2f\n", xAB1_);
  printf("xAB2_ : %.2f\n", xAB2_);
  printf("xA1_  : %.2f\n", xA1_);
  printf("xA2_  : %.2f\n", xA2_);
  printf("xA3_  : %.2f\n", xA3_);
  printf("xA4_  : %.2f\n", xA4_);
  printf("xA5_  : %.2f\n", xA5_);
  printf("xB1_  : %.2f\n", xB1_);
  printf("xB2_  : %.2f\n", xB2_);
  printf("xB3_  : %.2f\n", xB3_);
  printf("xB4_  : %.2f\n", xB4_);
  printf("xB5_  : %.2f\n", xB5_);
  printf("\n");

  // Convert to Minkowski metric
	for (int i = 0; i < 4; i++)
	{
		pA_[e2m[i]]  = (*pA)[i];
		pB_[e2m[i]]  = (*pB)[i];
		pA1_[e2m[i]] = (*pA1)[i];
		pA2_[e2m[i]] = (*pA2)[i];
		pA3_[e2m[i]] = (*pA3)[i];
		pB1_[e2m[i]] = (*pB1)[i];
		pB2_[e2m[i]] = (*pB2)[i];
		pB3_[e2m[i]] = (*pB3)[i];

	}

	decay.DisplayMomenta();

   double PSWeight_A123 = decay.DecayA123->GetPhaseSpaceWeight(xA1_, xA2_, xA3_, xA4_, xA5_);

	// H(p) -> e-(p3) vebar(p4) vmu(p5) mu+(p6) vtau(p7) vtaubar(p8)                 
	// double rh_6f_val = rh_6f_(p3_,p4_,p5_,p6_,p7_,p8_);
	//
	rh_6f_(pA1_,pA2_,pB2_,pB1_,pA3_,pB3_);

	h_tautau.ReadInCMatrix_2_2(taumatrices_.ch_tautau);
	taum.ReadInCMatrix_2_2(    taumatrices_.cdec_taum);
	taup.ReadInCMatrix_2_2(    taumatrices_.cdec_taup);
	c7_568.ReadInCMatrix_2_2(  taumatrices_.c7_568);

	c_amp_htautau.ReadInCMatrix_2_2(  tau_amplitudes_.c_amp_htautau  );
   c_amp_dec_taum.ReadInCMatrix_2_2( tau_amplitudes_.c_amp_dec_taum );
   c_amp_dec_taup.ReadInCMatrix_2_2( tau_amplitudes_.c_amp_dec_taup );
	c_amp_7_568.ReadInCMatrix_2_2(    tau_amplitudes_.c_amp_7_568    );
   c_amp_res.ReadInCMatrix_2_2(      tau_amplitudes_.c_amp_res    );

	h_tautau.Show();
   taum.Show();
   taup.Show();
   c7_568.Show();

	c_amp_htautau.Show();
	c_amp_dec_taum.Show();
	c_amp_dec_taup.Show();
	c_amp_7_568.Show();
	c_amp_res.Show();

//	decay.DisplayAllInfo();


#endif

// theta,phi sweep
// Start of loop
//  for (int x1Bin = 0; x1Bin < x1nBins; x1Bin++ )
//  for (int x2Bin = 0; x2Bin < x2nBins; x2Bin++ )
//  {
//  
//  xAB1 = x1Min + x1BinWidth*x1Bin;
//  xAB2 = x2Min + x2BinWidth*x2Bin;

  printf("xAB1: %.2f\n", xAB1);
  printf("xAB2: %.2f\n", xAB2);



#if 1

  printf("-------------------- Vegas test --------------------\n");

  // Calling Vegas
  Vegas(NDIM, NCOMP, Integrand, USERDATA, NVEC,
    EPSREL, EPSABS, VERBOSE, SEED,
    MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
    GRIDNO, STATEFILE, SPIN,
    &neval, &fail, integral, error, prob);

  printf("VEGAS RESULT:\tneval %d\tfail %d\n",
    neval, fail);

  	 printf("Component 0:\n");
    printf("VEGAS RESULT:\t%.8f +- %.8f\tp = %.3f\n",
      (double)integral[0], (double)error[0], (double)prob[0]);
  	 printf("Component 1:\n");
    printf("VEGAS RESULT:\t%.8f +- %.8f\tp = %.3f\n",
      (double)integral[1], (double)error[1], (double)prob[1]);
  	 printf("Component 2:\n");
    printf("VEGAS RESULT:\t%.8f +- %.8f\tp = %.3f\n",
      (double)integral[2], (double)error[2], (double)prob[2]);
  	 printf("Component 3:\n");
    printf("VEGAS RESULT:\t%.8f +- %.8f\tp = %.3f\n",
      (double)integral[3], (double)error[3], (double)prob[3]);
    printf("Component 4:\n");
    printf("VEGAS RESULT:\t%.8f +- %.8f\tp = %.3f\n",
      (double)integral[4], (double)error[4], (double)prob[4]);
    printf("Component 5:\n");
    printf("VEGAS RESULT:\t%.8f +- %.8f\tp = %.3f\n",
      (double)integral[5], (double)error[5], (double)prob[5]);



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

  	 printf("Component 0:\n");
    printf("SUAVE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
      (double)integral[0], (double)error[0], (double)prob[0]);
  	 printf("Component 1:\n");
    printf("SUAVE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
      (double)integral[1], (double)error[1], (double)prob[1]);
  	 printf("Component 2:\n");
    printf("SUAVE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
      (double)integral[2], (double)error[2], (double)prob[2]);
  	 printf("Component 3:\n");
    printf("SUAVE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
      (double)integral[3], (double)error[3], (double)prob[3]);
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

	// Get the results of the integrals
	integral_combined_pol_1 = integral[0];
	integral_combined_pol_2 = integral[1];
	integral_htautau_pol_1  = integral[2];
	integral_htautau_pol_2  = integral[3];
	integral_taudecay_pol_1 = integral[4];
	integral_taudecay_pol_2 = integral[5];

  ////////////////////////////////////////////////////////////////////////////
  // -- Cross-checking the Gamma(tau) value result with the textbook formula
  // (similar to the case of muon decay)
  double integral_taudecay_unpol = (integral_taudecay_pol_1 + integral_taudecay_pol_2);
  double DecayA123_PSConst = decay.DecayA123->GetPSConst();
  double Gamma_tau_unpol_woutgamma = 4.0 * integral_taudecay_unpol * DecayA123_PSConst * G_Fermi * G_Fermi / mA / 2.0  ;
  double Gamma_tau_pol_1_woutgamma = 8.0 * integral_taudecay_pol_1 * DecayA123_PSConst * G_Fermi * G_Fermi / mA / 2.0  ;
  double Gamma_tau_pol_2_woutgamma = 8.0 * integral_taudecay_pol_2 * DecayA123_PSConst * G_Fermi * G_Fermi / mA / 2.0  ;
  double ratio_Gamma_tau_unpol_woutgamma_result_per_expected = Gamma_tau_unpol_woutgamma/Gamma_formula;
  double ratio_Gamma_tau_pol_1_woutgamma_result_per_expected = Gamma_tau_pol_1_woutgamma/Gamma_formula;
  double ratio_Gamma_tau_pol_2_woutgamma_result_per_expected = Gamma_tau_pol_1_woutgamma/Gamma_formula;

  printf("\n");
  printf("############################################################\n");
  printf("### Decay process: h --> A + B --> (a1 a2 a3) (b1 b2 b3) ###\n");
  printf("############################################################\n");
  printf("\n");


  printf("\n");
  printf("--------------------------\n");
  printf("--- Physical constants ---\n");
  printf("--------------------------\n");
  printf("\n");
  printf("hbar:    %12.6e [GeV s]\n", hbar);
  printf("c:       %12.6e [m/s]\n", c);
  printf("hbar*c:  %12.6e [GeV fm]\n", hbar_c);
  printf("G_Fermi: %12.6e [GeV^{-2}]\n", G_Fermi);
  printf("m (tau): %12.6e [GeV]\n", m_tau);

  printf("\n");
  printf("\n");
  printf("##################################################\n");
  printf("### --- h -> (tau) (tau) consistency check --- ###\n");
  printf("#################################################\n");
  printf("# Checking 'M(h-> tau tau)' value\n");


  //double integral_htautau_unpol = (integral_htautau_pol_1 + integral_htautau_pol_2);
  double integral_htautau_unpol = (integral_htautau_pol_1 + integral_htautau_pol_2);
  double htautau_exp   = 0.5*((2*M*M)-(8*mA*mA));
  double htautau_ratio = (integral_htautau_unpol) /htautau_exp;

  printf("\n");
  printf("pol_1 = [0][0]; pol_2 = [0][1]");
  printf("integral_htautau_pol_1:        %12.6f [a.u.]\n", integral_htautau_pol_1);
  printf("integral_htautau_pol_2:        %12.6f [a.u.]\n", integral_htautau_pol_2);
  printf("integral_htautau_unpol:        %12.6f [a.u.]\n", integral_taudecay_unpol);
  printf("\n");

  printf("\n");
  printf("expected amplitude: 0.5*(2*(M_h)^(2) - 8*(m_tau)^2)\n");
  printf("\n");
  printf("h tautau amplitude value:  %12.6f\n", integral_htautau_unpol);
  printf("expected amplitude value:  %12.6f\n", htautau_exp);
  printf("\n");
  printf("ratio:\n");
  printf("ampl/expected:  %12.6e\n", htautau_ratio);


  printf("\n");
  printf("\n");
  printf("############################################\n");
  printf("### --- Gamma(tau) consistency check --- ###\n");
  printf("############################################\n");
  printf("# Comparing 'Gamma(tau) numerical int. result' with\n");
  printf("#       the 'Gamma(tau) textbook formula'.\n");

  printf("\n");
  printf("Unpolarized case\n");
  printf("integral_taudecay_unpol value:  %12.6f\n", integral_taudecay_unpol);
  printf("integral_taudecay_pol_1 value:  %12.6f\n", integral_taudecay_pol_1);
  printf("integral_taudecay_pol_2 value:  %12.6f\n", integral_taudecay_pol_2);
  printf("mA (daughter particle A): %12.6f [GeV]\n", mA);
  printf("\n");
  printf("Textbook result of the integration:\n");
  printf("Gamma(tau): G_Fermi^{2}*m_tau^{5}/(192*pi^{3})\n");
  printf("\n");
  printf("Our result of the integration:\n");
  printf("Gamma(tau): 4*integral_taudecay * G_Fermi * G_Fermi / mA / 2.0\n");
  printf("\n");

  printf("Gamma(tau) textbook:                                  			   %12.6e [GeV]\n", Gamma_formula);
  printf("Gamma(tau) unpolarized our result w/o gamma factor:   			   %12.6e [GeV]\n", Gamma_tau_unpol_woutgamma);

  printf("ratio [Gamma(tau)_unpol_result_woutgamma / Gamma(tau)_expected]: %12.6e\n", ratio_Gamma_tau_unpol_woutgamma_result_per_expected);
  printf("ratio [Gamma(tau)_pol_1_result_woutgamma / Gamma(tau)_expected]: %12.6e\n", ratio_Gamma_tau_pol_1_woutgamma_result_per_expected);
  printf("ratio [Gamma(tau)_pol_2_result_woutgamma / Gamma(tau)_expected]: %12.6e\n", ratio_Gamma_tau_pol_2_woutgamma_result_per_expected);



  printf("\n");
  printf("\n");

  ////////////////////////////////////////////////////////////////////////////
  // -- 
  // -- 
  
  // Note: These Gamma's are not weighted with the appropriate couplings
  //       and several other factors might be missing
  double Gamma_Combined_redu_pol_1 = 2*integral_combined_pol_1*(M_PI/(mA*integral_taudecay_unpol));
  double Gamma_Combined_redu_pol_2 = 2*integral_combined_pol_2*(M_PI/(mA*integral_taudecay_unpol));
  double Gamma_Combined_full_pol_1 = 2*integral_combined_pol_1*(4.0* DecayA123_PSConst * G_Fermi * G_Fermi )*(M_PI/(mA*Gamma_tau_unpol_woutgamma));
  double Gamma_Combined_full_pol_2 = 2*integral_combined_pol_2*(4.0* DecayA123_PSConst * G_Fermi * G_Fermi )*(M_PI/(mA*Gamma_tau_unpol_woutgamma));
  //double Gamma_Combined_NWA = integral_combined*(1.0/(integral_taudecay));
  double Gamma_Prod_BR_redu        = integral_htautau_unpol*M_PI/mA;
  double Gamma_Prod_BR_full        = integral_htautau_unpol*2.0*M_PI;

  printf("#################################################################\n");
  printf("### --- Comparison of 'Combined result' and '(Prod)x(BR)' --- ###\n");
  printf("#################################################################\n");
  printf("# Here we compare the result of the following two formulae:\n");
  printf("# o A) 'Combined result'\n");
  printf("#   (Formula A) = |M_{combined}|^2 * (pi)/(Gamma(tau)*m(tau))\n");
  printf("#   where Gamma(tau) is the numerical integration result.\n");
  printf("# o B) '(Production) x (BR)'\n");
  printf("#   (Formula B) = 2.0*pi |M_{prod}|^2 * BR(tau->3f)\n");
  printf("#   Both of them use NWA.\n");
  printf("#   These two formulae are integrated over the phase space\n");
  printf("#   of the daughter particle of A (a1, a2 and a3).\n");

  printf("\n");
  printf("integral_combined_pol_1: %12.6f [a.u.]\n", integral_combined_pol_1);
  printf("integral_combined_pol_2: %12.6f [a.u.]\n", integral_combined_pol_2);
  printf("integral_htautau_pol_1:        %12.6f [a.u.]\n", integral_htautau_pol_1);
  printf("integral_htautau_pol_2:        %12.6f [a.u.]\n", integral_htautau_pol_2);
  printf("integral_htautau_unpol:        %12.6f [a.u.]\n", integral_htautau_unpol);
  printf("integral_taudecay_pol_1:       %12.6f [a.u.]\n", integral_taudecay_pol_1);
  printf("integral_taudecay_pol_2:       %12.6f [a.u.]\n", integral_taudecay_pol_2);
  printf("mA:                      %12.6f [GeV]\n", mA);

  printf("\n");
  printf("----------------------------------------------\n");
  printf("--- With reduced calculation on Gamma(tau) ---\n");
  printf("-----------------------------------------------\n");
  printf("  no normalization to Gamma(tau) and |M| with\n");
  printf("  coupling constants, phase space constants, etc.\n");
  printf("  i.e taking directly the 'integral_taudecay' value\n");
  printf("  Note: The ratios should be equal to the 'full calculation' case!\n");
  printf("\n");
  printf("Gamma_Combined_redu_pol_1 value:     %12.6e [a.u.]\n", Gamma_Combined_redu_pol_1);
  printf("Gamma_Combined_redu_pol_2 value:     %12.6e [a.u.]\n", Gamma_Combined_redu_pol_2);
  printf("Gamma_Prod_BR_redu value:            %12.6e [a.u.]\n", Gamma_Prod_BR_redu);
  printf("\n");
  printf("ratios:\n");
  printf("(Gamma_Combined_redu_pol_1)/(Gamma_Prod_BR_redu): %12.6e [a.u.]\n", (Gamma_Combined_redu_pol_1)/Gamma_Prod_BR_redu);
  printf("(Gamma_Combined_redu_pol_2)/(Gamma_Prod_BR_redu): %12.6e [a.u.]\n", (Gamma_Combined_redu_pol_2)/Gamma_Prod_BR_redu);

  printf("\n");
  printf("\n");
  printf("-------------------------------------------\n");
  printf("--- With full calculation on Gamma(tau) ---\n");
  printf("-------------------------------------------\n");
  printf("  all normalization to Gamma(tau) and |M| with\n");
  printf("  coupling constants, phase space constants, etc.\n");
  printf("  i.e taking directly the 'Gamma(tau)' value\n");
  printf("  Note: The ratios should be equal to the 'reduced calculation' case!\n");

  printf("\n");
  printf("Gamma_Combined_full_pol_1 value:     %12.6e [a.u.]\n", Gamma_Combined_full_pol_1);
  printf("Gamma_Combined_full_pol_2 value:     %12.6e [a.u.]\n", Gamma_Combined_full_pol_2);
  printf("Gamma_Prod_BR_full value:            %12.6e [a.u.]\n", Gamma_Prod_BR_full);
  printf("\n");
  printf("ratios:\n");
  printf("(Gamma_Combined_full_pol_1)/(Gamma_Prod_BR_full): %12.6e [a.u.]\n", (Gamma_Combined_full_pol_1)/Gamma_Prod_BR_full);
  printf("(Gamma_Combined_full_pol_2)/(Gamma_Prod_BR_full): %12.6e [a.u.]\n", (Gamma_Combined_full_pol_2)/Gamma_Prod_BR_full);

// --- Expected result if the weight is only PSWeight --- ///
//// Get phase space constant from the DecayChain126 class
////double PSConst = decay.DecayA123->GetPSConst();
//
//// Expected value in the massless scenario:
//double threebodyA123 = mA*mA/(256.0*M_PI*M_PI*M_PI);
//double expected_result = threebodyA123 ;
//
//  printf("\n");
//  printf("                                           final result: %10.10e\n", result);
//  printf("Expexted result in massless (m1 = m2 = m3 = 0) scenario: %10.10e\n", expected_result);
//  printf("                                            their ratio: %10.10f\n", result/expected_result);
//  printf("\n");

//  double PSWeight_B123 = decay.DecayB123->GetPhaseSpaceWeight(xB1, xB2, xB3, xB4, xB5);
//  printf("PSWeight_B123 = %10.5e \n", PSWeight_B123);

///////////////////
//  --- Root --- //
///////////////////

//  hist_Gamma_Combined_full_pol_1 -> SetBinContent( (x1nBins-x1Bin), 1+x2Bin, Gamma_Combined_full_pol_1);
//  hist_Gamma_Combined_full_pol_2 -> SetBinContent( (x1nBins-x1Bin), 1+x2Bin, Gamma_Combined_full_pol_2);
//  hist_Gamma_Prod_BR_full        -> SetBinContent( (x1nBins-x1Bin), 1+x2Bin, Gamma_Prod_BR_full );


  // End of loops
  // }

//	hist_Gamma_Combined_full_pol_1 -> Divide( hist_Gamma_Prod_BR_full);
//	hist_Gamma_Combined_full_pol_2 -> Divide( hist_Gamma_Prod_BR_full);
//
//   gStyle->SetOptStat(0);
//   TCanvas canv ("canvas","canv", 800, 600);
//
//	canv.SetRightMargin(0.20);
//
//	hist_Gamma_Combined_full_pol_1 -> Draw("COLZ");
////	hist_Gamma_Combined_full_pol_2 			
////	hist_Gamma_Prod_BR_full        
//	
//	canv.SaveAs("histo_Gamma_Combined_full_pol1.pdf");
//
//	canv.Clear();
//	hist_Gamma_Combined_full_pol_2 -> Draw("COLZ");
//	canv.SaveAs("histo_Gamma_Combined_full_pol2.pdf");
//
//	canv.Clear();
//	hist_Gamma_Prod_BR_full -> Draw("COLZ");
//	canv.SaveAs("histo_Gamma_Prod_BR_full.pdf");

////////////////////////////////////////////

	//hist_Gamma_Combined_full_pol_1 -> SaveAs("./hist_Gamma_Combined_full_pol_1.pdf");
	//hist_Gamma_Combined_full_pol_2 -> SaveAs("./hist_Gamma_Combined_full_pol_2.pdf");  
	//hist_Gamma_Prod_BR_full        -> SaveAs("./hist_Gamma_Prod_BR_full.pdf");

  return 0;

}
