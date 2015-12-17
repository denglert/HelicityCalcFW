#ifndef PHASESPACETOOLS_H
#define PHASESPACETOOLS_H

#include <cmath>
#include <iostream>
#include <stdio.h>
#include "UtilFunctions.h"
#include <TLorentzVector.h>

#define DEBUG 0


// Two-body decay //
namespace TwoBodyFunc
{
	double beta(double s, double m1_sqr, double m2_sqr);
	double p   (double s, double m1_sqr, double m2_sqr);
	double E   (double s, double  m_sqr, double m_other_sqr);

}

//////////////////////
// Three-body decay //
namespace ThreeBody
{
	double x_to_sqrt_s (double x, double M, double m1, double m2, double m3 );
	double sqrt_s_to_x (double sqrt_s, double M, double m1, double m2, double m3 );

}

// General 
double x_to_cos (double x);
double cos_to_x (double cos);

double x_to_phi (double x);
double phi_to_x (double cos);


//////////////////////////
// class ThreeBodyDecay //
class ThreeBodyDecay
{

	public:
	ThreeBodyDecay(double M, double m1, double m2, double m3);
	~ThreeBodyDecay();


	// OBSOLETE
	//void SetMotherMass(double mass);
	//void SetDecayMass (int i, double mass);
	void SetMotherEPxPyPz(double E, double px, double py, double pz);
	void SetMotherMPThetaPhi(double M, double mom, double theta, double phi);
	void SetBitBoostBack(bool flag);
	double GetPhaseSpaceWeight(double x1, double x2, double x3, double x4, double x5);
	double GetJacobian();
	double GetPSConst();
	void SetPhaseSpace(double x1, double x2, double x3, double x4, double x5);

	void DisplayAll();
	void DisplayKinematics();
	void DisplayMomenta();
	void DisplayMasses();

	TLorentzVector *P;
	TLorentzVector **p;

	///////
	
	private:


	double M;
	double m[3];
	double M_sqr;
	double m_sqr[3];

	double Jacobian;
	double PSFactor;
	double PSWeight;
	double PSConst;

	bool BitBoostBack;

	double 	    betabar;        
	double 	     lambda;        
	double 	lambda_sqrt;        
	double 	  betabar23;        
	double 	   lambda23;        
	double lambda23_sqrt;        
	double 	        s23;        
	double 	   s23_sqrt;        
	double 	    s23_min;        
	double    s23_length;        
	double 	  costheta1;  
	double 	  sintheta1;  
	double 	       phi1;  
	double 	    cosphi1;  
	double 	    sinphi1;  
	double 	 costheta23; 
	double 	 sintheta23; 
	double 	      phi23;
	double 	   cosphi23;  
	double 	   sinphi23;  

	// In the rest frame of the mother particle
	double E1;
	double E23;
	double pmag;

	// In the rest frame of 23
	double p2mag_23;
	double E2_23;
	double E3_23;

	TLorentzVector sum;
	TVector3 Beta23;
	TVector3 BetaVec;

};



#endif
