#include <iostream>
#include <cmath>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TLorentzVector.h>
#include <TGenPhaseSpace.h>
#include <TCanvas.h>
#include "HelicityFW.h"
#include "FortranInterface.h"
#include "PhysConst.h"

////////////////////////////////////////

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
	// -- Generate events with RAMBO -- ///
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
	// -- Generate events with TGenPhaseSpace -- //
	///////////////////////////////////////////////

	
	TLorentzVector Mother(0.0, 0.0, 0.0, 125.0);
	
	TGenPhaseSpace event;
	event.SetDecay(Mother, 2, masses);
	
	TH2F *h2 = new TH2F("h2","h2", 50,1.1,1.8, 50,1.1,1.8);

	double p1_[4];
	double p2_[4];
	
	// Eucledian to Minkowski indeces mapping
	int e2m[4] = {1, 2, 3, 0};
	
	for (int iEv = 0; iEv < nEvents; iEv++)
	{
	   Double_t weight = event.Generate();
	   TLorentzVector *p1 = event.GetDecay(0);
	   TLorentzVector *p2 = event.GetDecay(1);
//		std::cout << "TGenPhaseSpace values" << std::endl;
//		std::cout << Form("p1[0]: %+2.2f, p1[1]: %.2f p1[2]: %.2f p1[3]: %.2f", (*p1)[0], (*p1)[1], (*p1)[2], (*p1)[3]) << std::endl;
//		std::cout << Form("p2[0]: %+2.2f, p2[1]: %.2f p2[2]: %.2f p2[3]: %.2f", (*p2)[0], (*p2)[1], (*p2)[2], (*p2)[3]) << std::endl;

		pDistr[0]->Fill( p1->Px() );
		pDistr[1]->Fill( p1->Py() );
		pDistr[2]->Fill( p1->Pz() );

		for (int i = 0; i < 4; i++)
		{
			p1_[e2m[i]] = (*p1)[i];
			p2_[e2m[i]] = (*p2)[i];
		}
//			std::cout << Form("p1_[0]: %+2.2f, p1_[1]: %.2f p1_[2]: %.2f p1_[3]: %.2f", p1_[0], p1_[1], p1_[2], p1_[3]) << std::endl;

		double value = rh_tautau_(p1_,p2_);
		std::cout << Form("rh_tatau: %.2f", value) << std::endl;

		ampl_sqr->Fill( p1->Px(), p1->Py(), p1->Pz(), value);
	}

	std::string filebase;
	std::string figOUTpng;
	std::string figOUTpdf;

	TCanvas canvas ("canvas","", 800, 600);

	for (int i = 0; i < 3; i++)
	{
  		filebase = Form("./pDistr_%d", i);
		figOUTpng = filebase+".png";
		figOUTpdf = filebase+".pdf";
		pDistr[i]->Draw();
// 	canvas.SaveAs(figOUTpng.c_str());
		canvas.SaveAs(figOUTpdf.c_str());
		canvas.Clear();

	}


		ampl_sqr->Draw();
		filebase = "./ampl_sqr";
		figOUTpng = filebase+".png";
		figOUTpdf = filebase+".pdf";
		canvas.SaveAs(figOUTpdf.c_str());
		canvas.Clear();


}
