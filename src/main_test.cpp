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

int main( int argc, const char *argv[] )
{

	std::cout << "Main program running.." << std::endl;

	////////////////////////
	// 					    //
	// -- Initializing -- //
	// 				       //
	////////////////////////
	

	int nEvents = 10;

	// Event parameters
	int nParticles = 2;
	double E_cm = 125.; 					  // GeV
	double masses[2] = {m_tau, m_tau}; // GeV
	masses_.rmtau = m_tau;

	const int pDistr_nBins = 50;
	const double pDistr_pmin = -70.;
	const double pDistr_pmax =  70.;

	TH1D *pDistr[3];
	
	for (int i = 0; i < 3; i++)
	{
		pDistr[i] = new TH1D(Form("pDistr[%d]", i),Form("%d;p_{}[GeV/c];", i), pDistr_nBins, pDistr_pmin, pDistr_pmax);
	}

	TH3D *ampl_sqr = new TH3D("ampl_sqr",";p_{x};p_{y};p_{z}",pDistr_nBins,pDistr_pmin,pDistr_pmax,pDistr_nBins,pDistr_pmin,pDistr_pmax,pDistr_nBins,pDistr_pmin,pDistr_pmax);

	///////////////////////////////////////
	// -- Generate HiggsDecays with RAMBO -- ///
	///////////////////////////////////////

	int lenv = 1;
	float rvec[1] = {-1.0};

	int nLoop = 5;
	for (int i = 0; i < nLoop; i++)
	{
		ranmar_(rvec,&lenv);
		double rand = rvec[0];
//		std::cout << "ranmar: " << rand << std::endl;
	}

	double **p = new double*[4];
	for (int i = 0; i < 4; i++)
	{
		p[i] = new double[nParticles];
		for (int iPart = 0; iPart < nParticles; iPart++)
		{
			p[i][iPart] = -1.0;
//			std::cout << "p[" << i << "][" << iPart << "]: " << p[i][iPart] << std::endl;
		}
	}

	double wt = 1.0;
	int lw = 1;
	
	for (int iEv = 0; iEv < nEvents; iEv++)
	{
//		std::cout << "RAMBO values" << std::endl;
		rambo_( &nParticles, &E_cm, masses, &p[0][0], &wt, &lw);

		for (int iPart = 0; iPart < nParticles; iPart++)
		{
			double E = masses[0]*masses[0];
			for (int i = 0; i < 3; i++)
			{
//				std::cout << "p[" << i << "][" << iPart << "]: " << p[iPart][i] << std::endl;
				E = E + p[iPart][i]*p[iPart][i];
			}
			E = sqrt(E);
//			std::cout << "p[" << 3 << "][" << iPart << "]: " << p[iPart][3] << std::endl;
//			std::cout << "E: " << E << std::endl;
		}
	}


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
	
	// Eucledian to Minkowski indeces mapping
	int e2m[4] = {1, 2, 3, 0};

	std::cout << "Process" << std::endl;
	std::cout << "H -->  (tau+) (tau-) --> (3 leptons) (3 leptons)" << std::endl;
	
	for (int iEv = 0; iEv < nEvents; iEv++)
	{

		std::cout << Form("iEv: %d", iEv) << std::endl;

		// Higgs Decay
	   Double_t weight = HiggsDecay.Generate();

		// Get tau+ and tau-
	   TLorentzVector *TauPos = HiggsDecay.GetDecay(0);
	   TLorentzVector *TauNeg = HiggsDecay.GetDecay(1);

		TVector3 TauPos_BoostVector = TauPos->BoostVector();
		TVector3 TauNeg_BoostVector = TauNeg->BoostVector();

		// Debuggg 
//		std::cerr << "TauPos_BoostVector" << std::endl;
//		std::cerr << Form("%4.2f %4.2f %4.2f", TauPos_BoostVector.x(), TauPos_BoostVector.y(), TauPos_BoostVector.z()) << std::endl;
//		std::cerr << Form("beta: %4.6f \n", sqrt(TauPos_BoostVector.x()*TauPos_BoostVector.x() + TauPos_BoostVector.y()*TauPos_BoostVector.y() + TauPos_BoostVector.z()*TauPos_BoostVector.z()) );
//		std::cerr << Form("beta: %4.6f \n" ,TauPos->Beta());

		// Make taus decay
		TauPosDecay.SetDecay( (*TauPos), 3, LeptonMasses1 );
		TauNegDecay.SetDecay( (*TauNeg), 3, LeptonMasses2 );

	   Double_t weight1 = TauPosDecay.Generate();
	   Double_t weight2 = TauNegDecay.Generate();

	   TLorentzVector *TauPos_Daughter_0 = TauPosDecay.GetDecay(0);
	   TLorentzVector *TauPos_Daughter_1 = TauPosDecay.GetDecay(1);
	   TLorentzVector *TauPos_Daughter_2 = TauPosDecay.GetDecay(2);
	   TLorentzVector *TauNeg_Daughter_0 = TauNegDecay.GetDecay(0);
	   TLorentzVector *TauNeg_Daughter_1 = TauNegDecay.GetDecay(1);
	   TLorentzVector *TauNeg_Daughter_2 = TauNegDecay.GetDecay(2);


	   TLorentzVector *sum = new TLorentzVector;
		(*sum) = (*TauPos_Daughter_0) + (*TauPos_Daughter_1) + (*TauPos_Daughter_2) + (*TauNeg_Daughter_0) + (*TauNeg_Daughter_1) + (*TauNeg_Daughter_2);

//		TauPos_Daughter_0->Boost(-TauPos_BoostVector);
//		TauPos_Daughter_1->Boost(-TauPos_BoostVector);
//		TauPos_Daughter_2->Boost(-TauPos_BoostVector);

	   TLorentzVector *TauPosDaughterRestFrameSum = new TLorentzVector;
		(*TauPosDaughterRestFrameSum) =  (*TauPos_Daughter_0) + (*TauPos_Daughter_1) + (*TauPos_Daughter_2);

		std::cout << "sum: "; displayTLorentzVector(sum);
		std::cout << "TauPosDaughterRestFrameSum "; displayTLorentzVector(TauPosDaughterRestFrameSum);
		std::cout << "TauPos "; displayTLorentzVector(TauPos);
		std::cout << "TauNeg "; displayTLorentzVector(TauNeg);
		std::cout << "TauPos_Daughter_0 "; displayTLorentzVector(TauPos_Daughter_0);
		std::cout << "TauPos_Daughter_1 "; displayTLorentzVector(TauPos_Daughter_1);
		std::cout << "TauPos_Daughter_2 "; displayTLorentzVector(TauPos_Daughter_2);
		std::cout << "TauNeg_Daughter_0 "; displayTLorentzVector(TauNeg_Daughter_0);
		std::cout << "TauNeg_Daughter_1 "; displayTLorentzVector(TauNeg_Daughter_1);
		std::cout << "TauNeg_Daughter_2 "; displayTLorentzVector(TauNeg_Daughter_2);

		pDistr[0]->Fill( TauPos->Px() );
		pDistr[1]->Fill( TauPos->Py() );
		pDistr[2]->Fill( TauPos->Pz() );

		for (int i = 0; i < 4; i++)
		{
			p1_[e2m[i]] = (*TauPos)[i];
			p2_[e2m[i]] = (*TauNeg)[i];
		}
//			std::cout << Form("p1_[0]: %+2.2f, p1_[1]: %.2f p1_[2]: %.2f p1_[3]: %.2f", p1_[0], p1_[1], p1_[2], p1_[3]) << std::endl;

		int pol1, pol2;
		pol1 = 1; pol2 = 1;
		double value = rh_tautau_(p1_,p2_,&pol1,&pol2);
		std::cout << Form("rh_tatau: %.2f", value) << std::endl;

		ampl_sqr->Fill( TauPos->Px(), TauPos->Py(), TauPos->Pz(), value);
	}

	std::string filebase;
	std::string figOUTpng;
	std::string figOUTpdf;

	TCanvas canvas ("canvas","", 800, 600);

	for (int i = 0; i < 3; i++)
	{
  		filebase = Form("./output/pDistr_%d", i);
		figOUTpng = filebase+".png";
		figOUTpdf = filebase+".pdf";
		pDistr[i]->Draw();
// 	canvas.SaveAs(figOUTpng.c_str());
		canvas.SaveAs(figOUTpdf.c_str());
		canvas.Clear();

	}


	ampl_sqr->Draw();
	filebase = "./output/ampl_sqr";
	figOUTpng = filebase+".png";
	figOUTpdf = filebase+".pdf";
	canvas.SaveAs(figOUTpdf.c_str());
	canvas.Clear();

}
