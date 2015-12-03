#include <iostream>
#include <cmath>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TGenPhaseSpace.h>
#include <TCanvas.h>
#include "FortranInterface.h"
#include "PhysConst.h"
#include "UtilFunctions.h"
#include "HelicityTools.h"

int main( int argc, const char *argv[] )
{

	////////////////////////
	// 					    //
	// -- Initializing -- //
	// 				       //
	////////////////////////
	
	int nEvents = 3;

	// Event parameters
	int nParticles = 2;
	double E_cm = 125.; 					  // GeV
	double masses[2] = {m_tau, m_tau}; // GeV
	masses_.rmtau = m_tau;

//	HelicityTools ht;

	////////////////////////////////////////////////////
	// -- Generate HiggsDecays with TGenPhaseSpace -- //
	////////////////////////////////////////////////////
	
	TGenPhaseSpace HiggsDecay_mtau_zero;
	TGenPhaseSpace HiggsDecay_mtau_mtau;

	TLorentzVector Higgs(0.0, 0.0, 0.0, m_higgs);
	
	double TauMasses_zero[2] = {0, 0};
	double TauMasses_mtau[2] = {m_tau, m_tau};

	HiggsDecay_mtau_zero.SetDecay(Higgs, 2, TauMasses_zero);
	HiggsDecay_mtau_mtau.SetDecay(Higgs, 2, TauMasses_mtau);
	
	double p1_[4];
	double p2_[4];
	
	// Eucledian to Minkowski indeces mapping
	int e2m[4] = {1, 2, 3, 0};

	std::cout << "\n\n#########################" << std::endl;
	std::cout << "### Consistency check ###" << std::endl;
	std::cout << "#########################\n" << std::endl;

	std::cout << "Process:" << std::endl;
	std::cout << "H -->  (tau+) (tau-)\n\n";

	std::cout << "#########################" << std::endl;
	std::cout << "### Massless scenario ###" << std::endl;
	std::cout << "#########################\n" << std::endl;


	masses_.rmtau = 0;
	std::cout << Form("masses_.rmtau = %.2f\n", masses_.rmtau);
	std::cout << Form("TauMasses_zero[%d] = %.2f \n", 0, TauMasses_zero[0]);
	std::cout << Form("TauMasses_zero[%d] = %.2f \n", 1, TauMasses_zero[1]);

	std::cout << Form("The unpolarized amplitude should be: \n|M|^2 = 2*m_higgs^2 = %.2f\n", 2*m_higgs*m_higgs);
	
	for (int iEv = 0; iEv < nEvents; iEv++)
	{

		std::cout << Form("\niEv: %d", iEv) << std::endl;

		// Higgs Decay
	   Double_t weight = HiggsDecay_mtau_zero.Generate();

		// Get tau+ and tau-
	   TLorentzVector *TauPos = HiggsDecay_mtau_zero.GetDecay(0);
	   TLorentzVector *TauNeg = HiggsDecay_mtau_zero.GetDecay(1);

		for (int i = 0; i < 4; i++)
		{
			p1_[e2m[i]] = (*TauPos)[i];
			p2_[e2m[i]] = (*TauNeg)[i];
		}

//		std::cout << "Polarized case" << std::endl;
//		int pol1, pol2;
//		pol1 = 1; pol2 = 1;
//		std::cout << Form("(1,1): %2.2f \n", ht.rh_tautau_polarized(p1_,p2_,pol1,pol2));
//		pol1 = 1; pol2 = 2;
//		std::cout << Form("(1,2): %2.2f \n", ht.rh_tautau_polarized(p1_,p2_,pol1,pol2));
//		pol1 = 2; pol2 = 1;
//		std::cout << Form("(2,1): %2.2f \n", ht.rh_tautau_polarized(p1_,p2_,pol1,pol2));
//		pol1 = 2; pol2 = 2;
//		std::cout << Form("(2,2): %2.2f \n", ht.rh_tautau_polarized(p1_,p2_,pol1,pol2));
//		std::cout << "Unpolarized case" << std::endl;
//		std::cout << Form("%2.2f \n", ht.rh_tautau_unpolarized(p1_,p2_));


	}


	std::cout << "\n########################" << std::endl;
	std::cout << "### Massive scenario ###" << std::endl;
	std::cout << "########################\n" << std::endl;

	masses_.rmtau = m_tau;
	std::cout << Form("masses_.rmtau = %.2f\n", masses_.rmtau);
	std::cout << Form("TauMasses_mtau[%d] = %.2f \n", 0, TauMasses_mtau[0]);
	std::cout << Form("TauMasses_mtau[%d] = %.2f \n", 1, TauMasses_mtau[1]);

	std::cout << Form("The unpolarized amplitude should be:\n");
	std::cout << Form("|M|^2 = 2*m_higgs^2 - 4*m_tau^2) = %.2f\n",
					 2*m_higgs*m_higgs-4*m_tau*m_tau );
	std::cout << Form("|M|^2 = 2*m_higgs^2 - 8*m_tau^2) = %.2f\n",
					 2*m_higgs*m_higgs-8*m_tau*m_tau );

	for (int iEv = 0; iEv < nEvents; iEv++)
	{

	std::cout << Form("\niEv: %d", iEv) << std::endl;

		// Higgs Decay
	   Double_t weight = HiggsDecay_mtau_mtau.Generate();

		// Get tau+ and tau-
	   TLorentzVector *TauPos = HiggsDecay_mtau_mtau.GetDecay(0);
	   TLorentzVector *TauNeg = HiggsDecay_mtau_mtau.GetDecay(1);

		for (int i = 0; i < 4; i++)
		{
			p1_[e2m[i]] = (*TauPos)[i];
			p2_[e2m[i]] = (*TauNeg)[i];
		}

//		std::cout << "Polarized case" << std::endl;
//		int pol1, pol2;
//		pol1 = 1; pol2 = 1;
//		std::cout << Form("(1,1): %2.2f \n", ht.rh_tautau_polarized(p1_,p2_,pol1,pol2));
//		pol1 = 1; pol2 = 2;
//		std::cout << Form("(1,2): %2.2f \n", ht.rh_tautau_polarized(p1_,p2_,pol1,pol2));
//		pol1 = 2; pol2 = 1;
//		std::cout << Form("(2,1): %2.2f \n", ht.rh_tautau_polarized(p1_,p2_,pol1,pol2));
//		pol1 = 2; pol2 = 2;
//		std::cout << Form("(2,2): %2.2f \n", ht.rh_tautau_polarized(p1_,p2_,pol1,pol2));
//		std::cout << "Unpolarized case" << std::endl;
//		std::cout << Form("%2.2f \n", ht.rh_tautau_unpolarized(p1_,p2_));
	}

}
