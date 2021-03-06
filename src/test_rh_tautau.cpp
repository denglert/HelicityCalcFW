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
	
	int nEvents = 1;

	// Event parameters
	int nParticles = 2;
	double E_cm = 125.; 					  // GeV

	double theta1 = M_PI/2;
	double theta2 = M_PI/2;

	double phi1 = M_PI/2;
	double phi2 = 0;

	//double masses[2] = {m_tau, m_tau}; // GeV
	//masses_.rmtau = m_tau;
	//
	double masses[2] = {m_tau, m_tau}; // GeV
	masses_.rmtau = m_tau;

	CMatrix_2_2 cth;
	TauMatrix taumatrix;

	taumatrix.SetName("H->tau+ tau- taumatrix");
	taumatrix.SetnParts(2);

	taumatrix.Setk0(1.0,1.0,0.0,0.0);
	TLorentzVector k0(1.0, 0.0, 0.0, 1.0);

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

		// Warning!!!
		// Setting momenta directions manually
//		TauPos->SetTheta(theta1);
//		TauPos->SetPhi(phi1);
//
//		TauNeg->SetTheta(theta2);
//		TauNeg->SetPhi(phi2);
		
		displayTLorentzVector(TauPos);
		displayTLorentzVector(TauNeg);

		for (int i = 0; i < 4; i++)
		{
			p1_[e2m[i]] = (*TauPos)[i];
			p2_[e2m[i]] = (*TauNeg)[i];
		}

		double unpolarizedamp = rh_tautau_(p1_,p2_,cth);
		taumatrix.ReadInCMatrix_2_2( cth );

		std::cout << "CPP TauMatrix class:" << std::endl;
		taumatrix.Setpi(0,p1_[0],p1_[1],p1_[2],p1_[3]);
		taumatrix.Setpi(1,p2_[0],p2_[1],p2_[2],p2_[3]);

		taumatrix.Show();
		taumatrix.ShowSumOfSquares();
		taumatrix.CalcUnpolarizedAmp();
		taumatrix.ShowUnpolarizedAmp();

		std::cout << "Fortran rh_tautau:" << std::endl;
		std::cout << Form("unpolarized amplitude (computed directly in rh_tautau): %.4f", unpolarizedamp) << std::endl;

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
//	std::cout << Form("|M|^2 = 2*m_higgs^2 - 4*m_tau^2) = %.2f\n",
//					 2*m_higgs*m_higgs-4*m_tau*m_tau );
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

		// Warning!
		// Setting momenta directions manually
//		TauPos->SetTheta(theta1);
//		TauPos->SetPhi(phi1);
//
//		TauNeg->SetTheta(theta2);
//		TauNeg->SetPhi(phi2);

		displayTLorentzVector(TauPos);
		displayTLorentzVector(TauNeg);

		for (int i = 0; i < 4; i++)
		{
			p1_[e2m[i]] = (*TauPos)[i];
			p2_[e2m[i]] = (*TauNeg)[i];
		}

		double unpolarizedamp = rh_tautau_(p1_,p2_,cth);
		taumatrix.ReadInCMatrix_2_2( cth );

		std::cout << "CPP TauMatrix class:" << std::endl;
		taumatrix.Setpi(0,p1_[0],p1_[1],p1_[2],p1_[3]);
		taumatrix.Setpi(1,p2_[0],p2_[1],p2_[2],p2_[3]);

		taumatrix.Show();
		taumatrix.Showk0();
		taumatrix.ShowSumOfSquares();
		taumatrix.CalcUnpolarizedAmp();
		taumatrix.ShowUnpolarizedAmp();

		std::cout << "Fortran rh_tautau:" << std::endl;
		std::cout << Form("unpolarized amplitude (computed directly in rh_tautau): %.4f", unpolarizedamp) << std::endl;

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
