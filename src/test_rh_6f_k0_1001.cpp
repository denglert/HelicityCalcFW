////////////////////////////
// test_rh_6f_k0_1001.cpp //
////////////////////////////

// Auxiliary vector 'k0' is taken as:
//	TLorentzVector k0 (1.0, 0.0, 0.0, 1.0);

#include <iostream>
#include <cmath>
#include <complex>
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
#include "HelicityTools.h"

// Eucledian to Minkowski indeces mapping
int e2m[4] = {1, 2, 3, 0};

/////////////////////////////////////////////////
int main( int argc, const char *argv[] )
{

	////////////////////////
	// 					  //
	// -- Initializing -- //
	// 				      //
	////////////////////////
	

	int nEvents = 1;

	// Event parameters
	int nParticles = 2;

	// Auxiliary vector(s)
	TLorentzVector k0 (1.0, 0.0, 0.0, 1.0);

	// Fortran COMMON blocks
	masses_.rmtau             = m_tau;
	couplings_.wcl            = 1.0;
	couplings_.gh_tautau      = 1.0;
	amplitudes_.rh_6f_tautau  = 0.0;
	amplitudes_.rh_6f_taum    = 0.0;
	amplitudes_.rh_6f_taup    = 0.0;
	amplitudes_.rh_6f_res     = 0.0;
	amplitudes_.rh_6f_res_nwa = 0.0;

	// Ampltiudes
	// CMatrix_2_2 = cval cdec_taum[2][2];
	CMatrix_2_2 cdec_taum;
	CMatrix_2_2 cdec_taup;
	CMatrix_2_2 ch_tautau;

	// TauMatrix
	TauMatrix h_tautau;
	TauMatrix taum;
	TauMatrix taup;
	TauMatrix c7_568;
	h_tautau.SetName("H ---> tau- tau+");
	taum.SetName("taum");
	taup.SetName("taup");
	c7_568.SetName("c7_568");

	////////////////////////////////////////////////////
	// -- Generate HiggsDecays with TGenPhaseSpace -- //
	////////////////////////////////////////////////////
	
	// Higgs
	TLorentzVector Higgs(0.0, 0.0, 0.0, m_higgs);
	
	// Masses
	double TauMasses[2] = {m_tau, m_tau};
	double LeptonMasses1[3] = { m_nu_tau, m_ele, m_nu_ele };
	double LeptonMasses2[3] = { m_nu_tau, m_muo, m_nu_muo };

	// Phase spaces
	TGenPhaseSpace HiggsDecay;
	TGenPhaseSpace p568Decay;
	TGenPhaseSpace p734Decay;

	// Higgs decay
	HiggsDecay.SetDecay(Higgs, 2, TauMasses);
	
	// Momenta
	double p1_[4];
	double p2_[4];
	double p3_[4];
	double p4_[4];
	double p5_[4];
	double p6_[4];
	double p7_[4];
	double p8_[4];

	std::cout << "\n\n#########################" << std::endl;
	std::cout << "### Consistency check ###" << std::endl;
	std::cout << "#########################\n" << std::endl;

	std::cout << "Process:" << std::endl;
	std::cout << "H -->  (tau+) (tau-) --> (3 leptons) (3 leptons)\n" << std::endl;

	std::cout << Form("m_higgs (used in phase space gen): %8.4f [GeV]", m_higgs) << std::endl;
	std::cout << Form("masses_.rmtau (used in rh_6f):     %8.4f [GeV]", masses_.rmtau) << std::endl;
	std::cout << Form("m_tau (used in phase space gen):   %8.4f [GeV]", m_tau) << std::endl;
	std::cout << Form("m_ele (used in phase space gen):   %8.4f [GeV]", m_ele) << std::endl;
	std::cout << Form("m_muo (used in phase space gen):   %8.4f [GeV]", m_muo) << std::endl;
	std::cout << Form("and all neutrinos are massless") << std::endl;
	
	for (int iEv = 0; iEv < nEvents; iEv++)
	{

		std::cout << Form("\n################\n### iEv: %3d ###\n################", iEv) << std::endl;

		// Generate Higgs Decay
	   Double_t weight = HiggsDecay.Generate();

		std::cout << "weight: " << weight << std::endl;

		// Get tau+ and tau-
	   TLorentzVector *p734 = HiggsDecay.GetDecay(1);
	   TLorentzVector *p568 = HiggsDecay.GetDecay(0);


		// Get BoostVector
		TVector3 p734_BoostVector = p734->BoostVector();
		TVector3 p568_BoostVector = p568->BoostVector();

		// Make taus decay
		p734Decay.SetDecay( (*p734), 3, LeptonMasses1 );
		p568Decay.SetDecay( (*p568), 3, LeptonMasses2 );

	   Double_t weight1 = p568Decay.Generate();
	   Double_t weight2 = p734Decay.Generate();

		std::cout << "weight1: " << weight1 << std::endl;
		std::cout << "weight2: " << weight2 << std::endl;

		// Get TLorentzVectors momenta
	   TLorentzVector *p7 = p734Decay.GetDecay(0); // v_tau
	   TLorentzVector *p3 = p734Decay.GetDecay(1); // ele-
	   TLorentzVector *p4 = p734Decay.GetDecay(2); // v_ele_bar
	   TLorentzVector *p8 = p568Decay.GetDecay(0); // v_tau_bar
	   TLorentzVector *p6 = p568Decay.GetDecay(1); // mu+
	   TLorentzVector *p5 = p568Decay.GetDecay(2); // v_muo
 
		// Creating (p k0) scalar products
 		double p568k0 = k0 * (*p568);
 		double p734k0 = k0 * (*p734);
 		double p3k0   = k0 * (*p3);
 		double p4k0   = k0 * (*p4); 
 		double p5k0   = k0 * (*p5);
 		double p6k0   = k0 * (*p6);
		double p7k0   = k0 * (*p7);
		double p8k0   = k0 * (*p8);

		// Sum of momenta
	   TLorentzVector *sum = new TLorentzVector;
		(*sum) = (*p8) + (*p6) + (*p5) + (*p7) + (*p3) + (*p4);

		std::cout << "\n# --- Generated momenta --- #\n";
		std::cout << "- sum: "; displayTLorentzVector(sum);
		std::cout << "- p734 "; displayTLorentzVector(p734);
		std::cout << "- p568 "; displayTLorentzVector(p568);
		std::cout << "- p7   (v_tau )"; displayTLorentzVector(p7);
		std::cout << "- p3       (e-)"; displayTLorentzVector(p3);
		std::cout << "- p4   (v_ebar)"; displayTLorentzVector(p4);
		std::cout << "- p8 (v_taubar)"; displayTLorentzVector(p8);
		std::cout << "- p6    (muon+)"; displayTLorentzVector(p6);
		std::cout << "- p5     (v_mu)"; displayTLorentzVector(p5);


		////////////////////////////////
		// -- Standard calculation -- //
		////////////////////////////////

		// -- Unpolarized term -- //
		double taum_amp = p734->Dot((*p4)) * p3->Dot((*p7));
		double taup_amp = p568->Dot((*p5)) * p6->Dot((*p8));

		// -- Polarized term -- //

		// polarization vectors
		TLorentzVector polvec_taum;
		TLorentzVector polvec_taup;

		//  standard polarization vector with p and k0
		for(int nu = 0; nu<4; nu++)
		{
			polvec_taum[nu] = ( (*p734)[nu]/m_tau) - (m_tau/p734k0)*k0[nu];
			polvec_taup[nu] = ( (*p568)[nu]/m_tau) - (m_tau/p568k0)*k0[nu];
		}

		double taum_pol_amp = m_tau*( polvec_taum * (*p4) ) * ( (*p3) * (*p7) ); // (ek')(qk)
		double taup_pol_amp = m_tau*( polvec_taup * (*p5) ) * ( (*p6) * (*p8) ); // (ek')(qk)


		std::cout << Form("\n# --- Standard calculation --- #") << std::endl;
		std::cout << Form("Unpolarized case:") << std::endl;
		std::cout << Form("- tau- (pk')(qk) (unpolarized): %.4f", taum_amp ) << std::endl;
		std::cout << Form("- tau+ (pk')(qk) (unpolarized): %.4f", taup_amp ) << std::endl;
		std::cout << Form("- tau-/tau+ (unpolarized): %.4f", taum_amp/taup_amp ) << std::endl;
		std::cout << Form("Polarized case:") << std::endl;

		std::cout << Form("- tau- m_tau*(ek')(qk) (polarized term): %.4f", taum_pol_amp ) << std::endl;
		std::cout << Form("- tau+ m_tau*(ek')(qk) (polarized term): %.4f", taup_pol_amp ) << std::endl;

		double taum_standard_pol_1 = taum_amp - taum_pol_amp;
		double taum_standard_pol_2 = taum_amp + taum_pol_amp;
		double taup_standard_pol_1 = taup_amp - taup_pol_amp;
		double taup_standard_pol_2 = taup_amp + taup_pol_amp;

		std::cout << Form("- taum_standard_pol_1 ((pk')(qk)-m_tau*(ek')(qk)): %.4f", taum_standard_pol_1) << std::endl;
		std::cout << Form("- taum_standard_pol_2 ((pk')(qk)+m_tau*(ek')(qk)): %.4f", taum_standard_pol_2) << std::endl;
		std::cout << Form("- taup_standard_pol_1 ((pk')(qk)-m_tau*(ek')(qk)): %.4f", taup_standard_pol_1) << std::endl;
		std::cout << Form("- taup_standard_pol_2 ((pk')(qk)+m_tau*(ek')(qk)): %.4f", taup_standard_pol_2) << std::endl;

		for (int i = 0; i < 4; i++)
		{
			p1_[e2m[i]] = (*p568)[i];
			p2_[e2m[i]] = (*p734)[i];
			p3_[e2m[i]] = (*p3)[i];
			p4_[e2m[i]] = (*p4)[i];
			p5_[e2m[i]] = (*p5)[i];
			p6_[e2m[i]] = (*p6)[i];
			p7_[e2m[i]] = (*p7)[i];
			p8_[e2m[i]] = (*p8)[i];
		}

//		double rh_tautau_val = rh_tautau_(p1_,p2_);
//		std::cout << Form("rh_tautau: %.2f\n", rh_tautau_val) << std::endl;

		//////////////////////////////////////////
		// -- Helicity amplitude calculation -- //
		//////////////////////////////////////////
		
		std::cout << Form("\n# --- Helicity amplitude calculation --- #") << std::endl;

		double rh_6f_val = rh_6f_(p3_,p4_,p5_,p6_,p7_,p8_,cdec_taum,cdec_taup,ch_tautau);

		// Read in tau matrices
		h_tautau.ReadInCMatrix_2_2( ch_tautau );
		taum.ReadInCMatrix_2_2(cdec_taum);
		taup.ReadInCMatrix_2_2(cdec_taup);

		c7_568.ReadInCMatrix_2_2(taumatrices_.c7_568);

		// Show tau matrices
		h_tautau.Show();
		h_tautau.ShowSum();
		taum.Show();
		taum.ShowSum();
		taup.Show();
		taup.ShowSum();
		c7_568.Show();

		std::cout << std::endl;

		// Calculate polarized terms
		double taum_rh_6f_pol_1 = std::abs(taum.m[1][0])*std::abs(taum.m[1][0]);
		double taum_rh_6f_pol_2 = std::abs(taum.m[1][1])*std::abs(taum.m[1][1]);

		double taup_rh_6f_pol_1 = std::abs(taup.m[0][1])*std::abs(taup.m[0][1]);
		double taup_rh_6f_pol_2 = std::abs(taup.m[1][1])*std::abs(taup.m[1][1]);

		double c7_568_pol_1 = std::abs(c7_568.m[1][0])*std::abs(c7_568.m[1][0]);
		double c7_568_pol_2 = std::abs(c7_568.m[1][1])*std::abs(c7_568.m[1][1]);

		taum_rh_6f_pol_1 = taum_rh_6f_pol_1/p3k0/p4k0/p7k0/p734k0;
		taum_rh_6f_pol_2 = taum_rh_6f_pol_2/p3k0/p4k0/p7k0/p734k0;

		taup_rh_6f_pol_1 = taup_rh_6f_pol_1/p5k0/p6k0/p8k0/p568k0;
		taup_rh_6f_pol_2 = taup_rh_6f_pol_2/p5k0/p6k0/p8k0/p568k0;

		c7_568_pol_1 = c7_568_pol_1/p3k0/p4k0/p7k0/p568k0;
		c7_568_pol_2 = c7_568_pol_2/p3k0/p4k0/p7k0/p568k0;

//		for (int i=0; i<2; i++)
//		for (int j=0; j<2; j++)
//		{
//		  std::cout << Form("rh_6f_tau-[%d][%d].r: %.4f \n",i,j,cdec_taum[i][j].r);
//		  std::cout << Form("rh_6f_tau-[%d][%d].i: %.4f \n",i,j,cdec_taum[i][j].i);
//		}
//
//		for (int i=0; i<2; i++)
//		for (int j=0; j<2; j++)
//		{
//		  std::cout << Form("rh_6f_tau+[%d][%d].r: %.4f \n",i,j,cdec_taup[i][j].r);
//		  std::cout << Form("rh_6f_tau+[%d][%d].i: %.4f \n",i,j,cdec_taup[i][j].i);
//		}

		double taum_unpol = taum.CalcSum();
		double taup_unpol = taup.CalcSum();

		taum_unpol = taum_unpol/p3k0/p4k0/p7k0/p734k0;
		taup_unpol = taup_unpol/p5k0/p6k0/p8k0/p568k0;

		std::cout << Form("- taum_unpolarized (calc via sum of square of tau elements): %.4f", taum_unpol) << std::endl;
		std::cout << Form("- taup_unpolarized (calc via sum of square of tau elements): %.4f", taup_unpol) << std::endl;

		std::cout << Form("- taum_rh_6f_pol_1 (calc via square): %.4f", taum_rh_6f_pol_1) << std::endl;
		std::cout << Form("- taum_rh_6f_pol_2 (calc via square): %.4f", taum_rh_6f_pol_2) << std::endl;
		std::cout << Form("- taup_rh_6f_pol_1 (calc via square): %.4f", taup_rh_6f_pol_1) << std::endl;
		std::cout << Form("- taup_rh_6f_pol_2 (calc via square): %.4f", taup_rh_6f_pol_2) << std::endl;
		std::cout << Form("- c7_568_rh_6f_pol_1 (calc via square): %.4f", c7_568_pol_1) << std::endl;
		std::cout << Form("- c7_568_rh_6f_pol_2 (calc via square): %.4f", c7_568_pol_2) << std::endl;

		std::cout << Form("- (taum_rh_6f_pol_1-taum_rh_6f_pol_2) (calc via square): %.4f", taum_rh_6f_pol_1-taum_rh_6f_pol_2) << std::endl;
		std::cout << Form("- (taup_rh_6f_pol_1-taup_rh_6f_pol_2) (calc via square): %.4f", taup_rh_6f_pol_1-taup_rh_6f_pol_2) << std::endl;


		std::cout << Form("- rh_6f_taum (passed through FORTRAN common): %.4f", amplitudes_.rh_6f_taum) << std::endl;
		std::cout << Form("- rh_6f_taup (passed through FORTRAN common): %.4f", amplitudes_.rh_6f_taup ) << std::endl;
		std::cout << Form("- rh_6f tau-/tau+ (passed through FORTRAN common): %.4f", amplitudes_.rh_6f_taum/amplitudes_.rh_6f_taup	) << std::endl;
		std::cout << Form("- rh_6f_res (passed through FORTRAN common): %.4f", amplitudes_.rh_6f_res ) << std::endl;
		std::cout << Form("- rh_6f_res_nwa (passed through FORTRAN common): %.4f", amplitudes_.rh_6f_res_nwa ) << std::endl;

		std::cout << Form("\n# --- Comparison --- #") << std::endl;
		std::cout << Form("- rh_6f_taum/taum ratio (unpolarized): %.4f", amplitudes_.rh_6f_taum/taum_amp ) << std::endl;
		std::cout << Form("- rh_6f_taup/taup ratio (unpolarized): %.4f", amplitudes_.rh_6f_taup/taup_amp ) << std::endl;

		std::cout << Form("- taum_pol_1 (rh_6f/standard ratio) (polarized): %.4f", taum_rh_6f_pol_1/taum_standard_pol_1 ) << std::endl;
		std::cout << Form("- taum_pol_2 (rh_6f/standard ratio) (polarized): %.4f", taum_rh_6f_pol_2/taum_standard_pol_2 ) << std::endl;
		std::cout << Form("- taup_pol_1 (rh_6f/standard ratio) (polarized): %.4f", taup_rh_6f_pol_1/taup_standard_pol_1 ) << std::endl;
		std::cout << Form("- taup_pol_2 (rh_6f/standard ratio) (polarized): %.4f", taup_rh_6f_pol_2/taup_standard_pol_2 ) << std::endl;


		//std::cout << Form("rh_6f: %.2f", rh_6f_val) << std::endl;
		//	double LeptonMasses1[3] = {m_nu_tau, m_ele, m_nu_ele};
		//	double LeptonMasses2[3] = {m_nu_tau, m_muo, m_nu_muo};
		// H(p) -> e-(p3) vebar(p4) vmu(p5) mu+(p6) vtau(p7) vtaubar(p8)                 
	}

}
