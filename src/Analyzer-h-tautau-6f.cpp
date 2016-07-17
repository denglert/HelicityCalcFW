#include <iostream>
#include <fstream>
#include <cmath>
#include <TROOT.h>
#include <TStyle.h>
#include <TH2D.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include "HelicityFW.h"
#include "FortranInterface.h"
#include "PhysConst.h"
#include "UtilFunctions.h"
#include "PhaseSpaceTools.h"
#include "HelicityTools.h"


///////////////////////////
// --- Configuration --- //
///////////////////////////

// -- Masses
const double M   = m_higgs;
const double mA  = m_tau;
const double mB  = m_tau;
const double mA1 = 0.00;
const double mA2 = 0.00;
const double mA3 = 0.00;
const double mB1 = 0.00;
const double mB2 = 0.00;
const double mB3 = 0.00;

// -- Mother particle momentum 
const double Px = 0.00;
const double Py = 0.00;
const double Pz = 0.00;
const double E  = sqrt( Px*Px + Py*Py + Pz*Pz + M*M );

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


/////////////////////////////
//// -- Main Function -- //// 
/////////////////////////////
//

int main() 
{

  // -- Banner -- //
  printf("\n");
  printf("############################################\n");
  printf("### --- h -> tau tau -> 6f  Analyzer --- ###\n");
  printf("############################################\n");
  printf("\n");


  // -- Initialization -- //
        h_tautau.SetName("h_tautau");
          c7_568.SetName("c7_568");
            taum.SetName("taum");
            taup.SetName("taup");
   c_amp_htautau.SetName("c_amp_htautau ");
  c_amp_dec_taum.SetName("c_amp_dec_taum");
  c_amp_dec_taup.SetName("c_amp_dec_taup");
     c_amp_7_568.SetName("c_amp_7_568");
       c_amp_res.SetName("c_amp_res");

  decay.SetTag("DecayChain126");

  // Settings Fortran COMMON blocks
  masses_.rmtau             = m_tau;
  couplings_.wcl            = 1.0;
  couplings_.gh_tautau      = 1.0;
//  amplitudes_.rh_6f_tautau  = 0.0;
//  amplitudes_.rh_6f_taum    = 0.0;
//  amplitudes_.rh_6f_taup    = 0.0;
//  amplitudes_.rh_6f_res     = 0.0;
//  amplitudes_.rh_6f_res_nwa = 0.0;

  P = decay.P;
  P->SetPxPyPzE(Px,Py,Pz,E);
  decay.SetBitBoostBack( BitBoostBack );

  // -- Set up pointers
  pA = decay.pA;
  pB = decay.pB;
  pA1 = decay.pA1;
  pA2 = decay.pA2;
  pA3 = decay.pA3;
  pB1 = decay.pB1;
  pB2 = decay.pB2;
  pB3 = decay.pB3;

  TRandom3 *r = new TRandom3();
 
  //  -- Phase Space Generator test -- //
  printf("\n");
  printf("###########################################################\n");
  printf("### --- Testing DecayChain126 phase space generator --- ###\n");
  printf("###########################################################\n");
  printf("\n");


  const double xAB1_ = 0.00;
  const double xAB2_ = 0.00;
  const double xA1_  = 0.001;
  const double xA2_  = 0.001;
  const double xA3_  = 0.001;
  const double xA4_  = 0.001;
  const double xA5_  = 0.001;
  const double xB1_  = 0.001;
  const double xB2_  = 0.001;
  const double xB3_  = 0.001;
  const double xB4_  = 0.001;
  const double xB5_  = 0.001;

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

	// H(p) -> e-(p3) vebar(p4) vmu(p5) mu+(p6) vtau(p7) vtaubar(p8)                 
	// double rh_6f_val = rh_6f_(p3_,p4_,p5_,p6_,p7_,p8_);
	rh_6f_(pA1_,pA2_,pB2_,pB1_,pA3_,pB3_);

	// Read In TauMatrices
	h_tautau.ReadInCMatrix_2_2(taumatrices_.ch_tautau);
	taum.ReadInCMatrix_2_2(    taumatrices_.cdec_taum);
	taup.ReadInCMatrix_2_2(    taumatrices_.cdec_taup);
	c7_568.ReadInCMatrix_2_2(  taumatrices_.c7_568);

	c_amp_htautau.ReadInCMatrix_2_2(  tau_amplitudes_.c_amp_htautau  );
   c_amp_dec_taum.ReadInCMatrix_2_2( tau_amplitudes_.c_amp_dec_taum );
   c_amp_dec_taup.ReadInCMatrix_2_2( tau_amplitudes_.c_amp_dec_taup );
	c_amp_7_568.ReadInCMatrix_2_2(    tau_amplitudes_.c_amp_7_568    );
   c_amp_res.ReadInCMatrix_2_2(      tau_amplitudes_.c_amp_res    );

	// Show TauMatrices
	h_tautau.Show();
   taum.Show();
   taup.Show();
   c7_568.Show();
	c_amp_htautau.Show();
	c_amp_dec_taum.Show();
	c_amp_dec_taup.Show();
	c_amp_7_568.Show();
	c_amp_res.Show();

	
	double amp_h_taum_pol_1;
	double amp_h_taum_pol_2;

	double amp_h_taup_pol_1;
	double amp_h_taup_pol_2;
	
	double amp_taum_dec_pol_1;
	double amp_taum_dec_pol_2;

	double amp_taup_dec_pol_1;
	double amp_taup_dec_pol_2;

	double prop_full_unpol;
	double prop_h_taum_unpol;
	double prop_h_taup_unpol;
	double prop_taum_dec_unpol;
	double prop_taup_dec_unpol;

	double amp_full_intact;

	double amp_taum_dec_unpol;
	double amp_taup_dec_unpol;

	double amp_h_taum_unpol;
	double amp_h_taup_unpol;


  /////////////////////////////////////
  // -- Random phase space points -- //
  /////////////////////////////////////

  // Output file
  std::ofstream output;	
  output.open("h_6f_random_pts.dat");


//  int nCalls = 10;
  int nCalls = 10000;
  for (int i = 0; i < nCalls; i++)
  {

	  double xAB1 = r->Uniform(0.0,1.0);
//	  double xAB1 = 0.98;
	  double xAB2 = r->Uniform(0.0,1.0);
	  double xA1  = r->Uniform(0.0,1.0);
	  double xA2  = r->Uniform(0.0,1.0);
	  double xA3  = r->Uniform(0.0,1.0);
	  double xA4  = r->Uniform(0.0,1.0);
	  double xA5  = r->Uniform(0.0,1.0);
	  double xB1  = r->Uniform(0.0,1.0);
	  double xB2  = r->Uniform(0.0,1.0);
	  double xB3  = r->Uniform(0.0,1.0);
	  double xB4  = r->Uniform(0.0,1.0);
	  double xB5  = r->Uniform(0.0,1.0);
	
	  double PSConst  = decay.PSConst;
	  double PSWeight = decay.GetPhaseSpaceWeight( xAB1, xAB2,
							 	   xA1,  xA2, xA3, xA4, xA5,
								   xB1,  xB2, xB3, xB4, xB5);

//	  printf("xAB1: %.3f\n", xAB1);
//	  printf("xAB2: %.3f\n", xAB2);
//	  printf("xA1:  %.3f\n", xA1);
//	  printf("xA2:  %.3f\n", xA2);
//	  printf("xA3:  %.3f\n", xA3);
//	  printf("xA4:  %.3f\n", xA4);
//	  printf("xA4:  %.3f\n", xA5);
//	  printf("xB1:  %.3f\n", xB1);
//	  printf("xB2:  %.3f\n", xB2);
//	  printf("xB3:  %.3f\n", xB3);
//	  printf("xB4:  %.3f\n", xB4);
//	  printf("xB5:  %.3f\n", xB5);

//	  decay.DisplayMomenta();

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


	// Phase space weights
	  double PSWeight_A123    = decay.PSWeight_DecayA123;
	  double PSWeight_B123    = decay.PSWeight_DecayB123;

	  double PSConst_AB = decay.DecayAB->GetPhaseSpaceWeight();
	  double PSConst_A123    = decay.DecayA123->GetPSConst();
	  double PSConst_B123    = decay.DecayB123->GetPSConst();

	// Call Phact code
	
	rh_6f_(pA1_,pA2_,pB2_,pB1_,pA3_,pB3_);

	// Read in tau matrices
	c_amp_htautau.ReadInCMatrix_2_2(  tau_amplitudes_.c_amp_htautau  );
   c_amp_dec_taum.ReadInCMatrix_2_2( tau_amplitudes_.c_amp_dec_taum );
   c_amp_dec_taup.ReadInCMatrix_2_2( tau_amplitudes_.c_amp_dec_taup );
	c_amp_7_568.ReadInCMatrix_2_2(    tau_amplitudes_.c_amp_7_568    );
   c_amp_res.ReadInCMatrix_2_2(      tau_amplitudes_.c_amp_res    );

//	c_amp_dec_taum.Show();
//	c_amp_dec_taup.Show();


	////////////////////////
	// --- Amplitudes --- //
	////////////////////////
	
	// -- Assigning amplitudes
	
	// Full amplitude containing spin correlation cross terms
	amp_full_intact     = std::abs( c_amp_res.m[1][1]) ;

	// Polarized amplitudes
	amp_taum_dec_pol_1 = std::abs( c_amp_dec_taum.m[1][0] );
	amp_taum_dec_pol_2 = std::abs( c_amp_dec_taum.m[1][1] );
	amp_taup_dec_pol_1 = std::abs( c_amp_dec_taup.m[0][1] );
	amp_taup_dec_pol_2 = std::abs( c_amp_dec_taup.m[1][1] );

	// Not sure which is taum which is taup
	amp_h_taum_pol_1 = std::abs( c_amp_htautau.m[0][0] );
	amp_h_taum_pol_2 = std::abs( c_amp_htautau.m[0][1] );
	amp_h_taup_pol_1 = std::abs( c_amp_htautau.m[1][0] );
	amp_h_taup_pol_2 = std::abs( c_amp_htautau.m[1][1] );

	double amp_h_tautau_11 = std::abs( c_amp_htautau.m[0][0] );
	double amp_h_tautau_12 = std::abs( c_amp_htautau.m[0][1] );
	double amp_h_tautau_21 = std::abs( c_amp_htautau.m[1][0] );
	double amp_h_tautau_22 = std::abs( c_amp_htautau.m[1][1] );

	// Still questionable
	double amp_A_1 = std::abs( c_amp_dec_taum.m[1][0] );
	double amp_A_2 = std::abs( c_amp_dec_taum.m[1][1] );
//	double amp_B_1 = std::abs( c_amp_dec_taup.m[0][1] );
//	double amp_B_2 = std::abs( c_amp_dec_taup.m[1][1] );
	double amp_B_1 = std::abs( c_amp_dec_taup.m[0][1] );
	double amp_B_2 = std::abs( c_amp_dec_taup.m[1][1] );


	// Unpolarized amplitudes
	amp_taum_dec_unpol = amp_taum_dec_pol_1 + amp_taum_dec_pol_2;
	amp_taup_dec_unpol = amp_taup_dec_pol_1 + amp_taup_dec_pol_2;
	// Is this the right way to do it ????????
	// Is this the right way to do it ????????
	// Is this the right way to do it ????????
	amp_h_taum_unpol   = amp_h_taum_pol_1 + amp_h_taum_pol_2;
	amp_h_taup_unpol   = amp_h_taup_pol_1 + amp_h_taup_pol_2;


	// Scaling with coupling constants & etc.
	prop_full_unpol     = G_Fermi*G_Fermi*G_Fermi*G_Fermi*PSConst      * PSWeight      * amp_full_intact;
	prop_h_taum_unpol   = PSConst_AB                   * amp_h_taum_unpol;
	prop_h_taup_unpol   = PSConst_AB                   * amp_h_taup_unpol;
	prop_taum_dec_unpol = G_Fermi*G_Fermi*PSConst_A123 * PSWeight_A123 * amp_taum_dec_unpol;
	prop_taup_dec_unpol = G_Fermi*G_Fermi*PSConst_B123 * PSWeight_B123 * amp_taup_dec_unpol;

	double full_calculation = prop_full_unpol * (1.0/(mA*Gamma_formula)) * (1.0/(mB*Gamma_formula)) / 4.0;
	double prod_x_BR        = (prop_h_taum_unpol + prop_h_taup_unpol) *
			  						  (prop_taum_dec_unpol) * (1.0/(2.0*mA*Gamma_formula)) *
			  						  (prop_taup_dec_unpol) * (1.0/(2.0*mB*Gamma_formula)) ;

	// Sum of individual spin locked amiplitude squared contributions
	double result_spin_locked_sum = (1.0/(2.0*mA*Gamma_formula))*(1.0/(2.0*mB*Gamma_formula))
			  								  * (G_Fermi*G_Fermi*G_Fermi*G_Fermi)
											  * PSConst_AB
											  * PSConst_A123 * PSWeight_A123
											  * PSConst_B123 * PSWeight_B123 
											  *
			  							     ( 
												amp_h_tautau_11*amp_A_1*amp_B_1 + 
												amp_h_tautau_12*amp_A_1*amp_B_2 + 
												amp_h_tautau_21*amp_A_2*amp_B_1 + 
												amp_h_tautau_22*amp_A_2*amp_B_2
											  ) ;


//	double amplitude_piecewise = 

	double ratio = full_calculation/prod_x_BR;

	double ratio_spin_locked = full_calculation/result_spin_locked_sum;
	
	double ratio_raw = (amp_full_intact ) /
			  				 ((amp_h_taum_unpol+amp_h_taup_unpol)*(amp_taum_dec_unpol)*(amp_taup_dec_unpol));

   // Debug info for PSConstants and weights
	// Debuggg 
//	std::cerr << " PSConst " << PSConst << std::endl;
//	std::cerr << " PSWeight " << PSWeight << std::endl;
//
//	std::cerr << " PSConst_AB " << PSConst_AB << std::endl;
//	std::cerr << " PSConst_A123 " << PSConst_A123 << std::endl;
//	std::cerr << " PSConst_B123 " << PSConst_B123 << std::endl;
//	std::cerr << " PSConst_AB*PSConst_A123*PSConst_B123" << PSConst_AB*PSConst_A123*PSConst_A123 << std::endl;
//
//	std::cerr << " PSWeight_A123 " << PSWeight_A123 << std::endl;
//	std::cerr << " PSWeight_B123 " << PSWeight_B123 << std::endl;
//
//	std::cerr << " PSWeight_A123*PSWeight_B123 " << PSWeight_A123*PSWeight_B123 << std::endl;
//	std::cerr << " PSConst_AB*PSWeight_A123*PSWeight_B123 " << PSConst_AB*PSWeight_A123*PSWeight_B123 << std::endl;
//
//	std::cerr << "---- Check ----" << std::endl;
//
//	std::cerr << " PSConst*PSWeight " << PSConst*PSWeight << std::endl;
//
//	std::cerr << " PSConst_AB*PSConst_A123*PSWeight_A123*PSConst_B123*PSWeight_B123 " << PSConst_AB*PSConst_A123*PSWeight_A123*PSConst_B123*PSWeight_B123 << std::endl;



//	printf("\n");
//	printf("xAB1: %.4f\n", xAB1);
//	printf("xAB2: %.4f\n", xAB2);
//	printf("xA1:  %.4f\n", xA1);
//	printf("xA2:  %.4f\n", xA2);
//	printf("xA3:  %.4f\n", xA3);
//	printf("xA4:  %.4f\n", xA4);
//	printf("xA5:  %.4f\n", xA5);
//	printf("xB1:  %.4f\n", xB1);
//	printf("xB2:  %.4f\n", xB2);
//	printf("xB3:  %.4f\n", xB3);
//	printf("xB4:  %.4f\n", xB4);
//	printf("xB5:  %.4f\n", xB5);
//
//	printf("\n");
////	decay.DisplayMomenta();
//	printf("\n");
//	printf("Full calculation: %16.8e\n", full_calculation);
//	printf("(Prod)x(BR):      %16.8e\n", prod_x_BR);
//	printf("ratio:            %16.8e\n", ratio);
//	printf("\n");

	output <<
//	xAB1 << " " <<
//	xAB1 << " " <<
//	xAB2 << " " <<
//	xA1  << " " <<
//	xA2  << " " <<
//	xA3  << " " <<
//	xA4  << " " <<
//	xA5  << " " <<
//	xB1  << " " <<
//	xB2  << " " <<
//	xB3  << " " <<
//	xB4  << " " <<
//	xB5  << " " <<
	full_calculation << " " <<
	prod_x_BR << " " <<
	ratio << " " <<
	ratio_raw << " " <<
	ratio_spin_locked <<
	std::endl;

  }

  return 0;

}
