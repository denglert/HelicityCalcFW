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
	

	int nEvents = 10;

	// Event parameters
	int nParticles = 2;
	double E_cm = 125.; 					  // GeV
	double masses[2] = {m_tau, m_tau}; // GeV

	// Fortran COMMON blocks
	masses_.rmtau             = m_tau;
	couplings_.wcl            = 1.0;
	couplings_.gh_tautau      = 1.0;
	amplitudes_.rh_6f_taup    = 0.0;
	amplitudes_.rh_6f_taum    = 0.0;
	amplitudes_.rh_6f_tautau  = 0.0;
	amplitudes_.rh_6f_res     = 0.0;
	amplitudes_.rh_6f_res_nwa = 0.0;

	///////////////////////////////////////////////
	// -- Generate HiggsDecays with TGenPhaseSpace -- //
	///////////////////////////////////////////////
	
	TGenPhaseSpace HiggsDecay;
	TLorentzVector Higgs(0.0, 0.0, 0.0, m_higgs);
	
	double TauMasses[2] = {m_tau, m_tau};
	TGenPhaseSpace TauPosDecay;
	TGenPhaseSpace TauNegDecay;

	double LeptonMasses1[3] = {m_nu_tau, m_ele, m_nu_ele};
	double LeptonMasses2[3] = {m_nu_tau, m_muo, m_nu_muo};

	HiggsDecay.SetDecay(Higgs, 2, TauMasses);
	
	TH2F *h2 = new TH2F("h2","h2", 50,1.1,1.8, 50,1.1,1.8);

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

	std::cout << "Process" << std::endl;
	std::cout << "H -->  (tau+) (tau-) --> (3 leptons) (3 leptons)\n" << std::endl;
	
	for (int iEv = 0; iEv < nEvents; iEv++)
	{

		std::cout << Form("################\n### iEv: %3d ###\n################", iEv) << std::endl;

		// Higgs Decay
	   Double_t weight = HiggsDecay.Generate();

		// Get tau+ and tau-
	   TLorentzVector *TauNeg = HiggsDecay.GetDecay(1);
	   TLorentzVector *TauPos = HiggsDecay.GetDecay(0);


		TVector3 TauNeg_BoostVector = TauNeg->BoostVector();
		TVector3 TauPos_BoostVector = TauPos->BoostVector();

		// Make taus decay
		TauNegDecay.SetDecay( (*TauNeg), 3, LeptonMasses1 );
		TauPosDecay.SetDecay( (*TauPos), 3, LeptonMasses2 );

	   Double_t weight2 = TauNegDecay.Generate();
	   Double_t weight1 = TauPosDecay.Generate();

	   TLorentzVector *TauNeg_Daughter_0 = TauNegDecay.GetDecay(0); // v_tau
	   TLorentzVector *TauNeg_Daughter_1 = TauNegDecay.GetDecay(1); // ele-
	   TLorentzVector *TauNeg_Daughter_2 = TauNegDecay.GetDecay(2); // v_ele
	   TLorentzVector *TauPos_Daughter_0 = TauPosDecay.GetDecay(0); // v_tau
	   TLorentzVector *TauPos_Daughter_1 = TauPosDecay.GetDecay(1); // mu+
	   TLorentzVector *TauPos_Daughter_2 = TauPosDecay.GetDecay(2); // v_muo

	   TLorentzVector *sum = new TLorentzVector;
		(*sum) = (*TauPos_Daughter_0) + (*TauPos_Daughter_1) + (*TauPos_Daughter_2) + (*TauNeg_Daughter_0) + (*TauNeg_Daughter_1) + (*TauNeg_Daughter_2);

		std::cout << "# Momenta #\n";
		std::cout << "- sum: "; displayTLorentzVector(sum);
		std::cout << "- TauNeg "; displayTLorentzVector(TauNeg);
		std::cout << "- TauPos "; displayTLorentzVector(TauPos);
		std::cout << "- TauNeg_Daughter_0   (v_tau )"; displayTLorentzVector(TauNeg_Daughter_0);
		std::cout << "- TauNeg_Daughter_1       (e-)"; displayTLorentzVector(TauNeg_Daughter_1);
		std::cout << "- TauNeg_Daughter_2   (v_ebar)"; displayTLorentzVector(TauNeg_Daughter_2);
		std::cout << "- TauPos_Daughter_0 (v_taubar)"; displayTLorentzVector(TauPos_Daughter_0);
		std::cout << "- TauPos_Daughter_1    (muon+)"; displayTLorentzVector(TauPos_Daughter_1);
		std::cout << "- TauPos_Daughter_2     (v_mu)"; displayTLorentzVector(TauPos_Daughter_2);


		double tauneg_amp = TauNeg->Dot((*TauNeg_Daughter_2)) * TauNeg_Daughter_1->Dot((*TauNeg_Daughter_0));
		double taupos_amp = TauPos->Dot((*TauPos_Daughter_2)) * TauPos_Daughter_1->Dot((*TauPos_Daughter_0));
		std::cout << Form("# Standard calculation #") << std::endl;
		std::cout << Form("- tau- (pk')(qk): %.4f", tauneg_amp ) << std::endl;
		std::cout << Form("- tau+ (pk')(qk): %.4f", taupos_amp ) << std::endl;
		std::cout << Form("- tau-/tau+: %.4f", tauneg_amp/taupos_amp ) << std::endl;


		for (int i = 0; i < 4; i++)
		{
			p1_[e2m[i]] = (*TauPos)[i];
			p2_[e2m[i]] = (*TauNeg)[i];
			p3_[e2m[i]] = (*TauNeg_Daughter_1)[i];
			p4_[e2m[i]] = (*TauNeg_Daughter_2)[i];
			p5_[e2m[i]] = (*TauPos_Daughter_2)[i];
			p6_[e2m[i]] = (*TauPos_Daughter_1)[i];
			p7_[e2m[i]] = (*TauNeg_Daughter_0)[i];
			p8_[e2m[i]] = (*TauPos_Daughter_0)[i];
		}

//		double rh_tautau_val = rh_tautau_(p1_,p2_);
//		std::cout << Form("rh_tautau: %.2f\n", rh_tautau_val) << std::endl;

		double rh_6f_val = rh_6f_(p3_,p4_,p5_,p6_,p7_,p8_);
		std::cout << Form("# Helicity amplitude calculation #") << std::endl;
		std::cout << Form("- rh_6f_tau-: %.4f", amplitudes_.rh_6f_taum ) << std::endl;
		std::cout << Form("- rh_6f_tau+: %.4f", amplitudes_.rh_6f_taup ) << std::endl;
		std::cout << Form("- rh_6f tau-/tau+: %.4f", amplitudes_.rh_6f_taum/amplitudes_.rh_6f_taup	) << std::endl;
		std::cout << Form("- rh_6f_res: %.4f", amplitudes_.rh_6f_res ) << std::endl;
		std::cout << Form("- rh_6f_res_nwa: %.4f", amplitudes_.rh_6f_res_nwa ) << std::endl;

		std::cout << Form("###") << std::endl;
		std::cout << Form("- rh_6f_taum/taum ratio: %.4f", amplitudes_.rh_6f_taum/tauneg_amp ) << std::endl;
		std::cout << Form("- rh_6f_taup/taup ratio: %.4f", amplitudes_.rh_6f_taup/taupos_amp ) << std::endl;
//		std::cout << Form("rh_6f: %.2f", rh_6f_val) << std::endl;
		//	double LeptonMasses1[3] = {m_nu_tau, m_ele, m_nu_ele};
		//	double LeptonMasses2[3] = {m_nu_tau, m_muo, m_nu_muo};
		// H(p) -> e-(p3) vebar(p4) vmu(p5) mu+(p6) vtau(p7) vtaubar(p8)                 
	}

}
