#ifndef PHASESPACETOOLS_H
#define PHASESPACETOOLS_H

#include <cmath>
#include <iostream>
#include <stdio.h>
#include "UtilFunctions.h"
#include <TLorentzVector.h>

#define DEBUG 0


////////////////////
// Two-body decay //
////////////////////

namespace TwoBodyFunc
{
	double beta(double s, double m1_sqr, double m2_sqr);
	double lambda(double a, double b, double c);
	double p   (double s, double m1_sqr, double m2_sqr);
	double E   (double s, double  m_sqr, double m_other_sqr);

}

//////////////////////
// Three-body decay //
//////////////////////

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

////////////////////////////////
// --- TwoBodyDecay class --- //
////////////////////////////////

class TwoBodyDecay
{

	//////
	public:
	TwoBodyDecay(double sqrt_s_, double m1_, double m2_);
	~TwoBodyDecay();

	void SetTag( const char tag_[] );
	void SetTotalThreeMomentum_pxpypz( double px, double py, double pz );
	void SetTotalThreeMomentum_pthetaphi(double p_, double theta, double phi_);

	void SetBitBoostBack(bool flag);
	void SetPhaseSpace(double x1_, double x2_);
	double GetPhaseSpaceWeight();

	void DisplayAll();
	void DisplayTag();
	void DisplayKinematics();
	void DisplayConfiguration();
	void DisplayMomenta();
	void DisplayMasses();

	TLorentzVector *P;
	TLorentzVector *p1;
	TLorentzVector *p2;
	TLorentzVector *sum;

	///////
	private:

	std::string tag;
	bool BitBoostBack;
	TVector3 BetaVec;

	double PSWeight_DecayAB;
	double PSWeight_DecayA123;
	double PSWeight_DecayB123;
	double PSWeight;

	double x1, x2;

	// System kinematics
	double E;
	double s;
	double sqrt_s;

	// Daughther particles
	double m1;
	double m2;
	double m1_sqr;
	double m2_sqr;
	
	// CoM Frame energy
	double E1;
	double E2;

	double pmag;
	double lambda;        
	double sqrt_lambda;        

	double costheta;
	double sintheta;
	double      phi;
	double   cosphi;
	double   sinphi;

};

//////////////////////////////////
// --- ThreeBodyDecay class --- //
//////////////////////////////////

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
	void SetTag( const char tag_[] );
	double GetPhaseSpaceWeight(double x1_, double x2_, double x3_, double x4_, double x5_);
	double GetJacobian();
	double GetPSConst();
	void SetPhaseSpace(double x1_, double x2_, double x3_, double x4_, double x5_);

	void DisplayAll();
	void DisplayTag();
	void DisplayKinematics();
	void DisplayConfiguration();
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

	double PSWeight;
	double PSConst;

	bool BitBoostBack;

	std::string tag;

	double x1, x2, x3, x4, x5;

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

/////////////////////////////////
// --- DecayChain126 class --- //
/////////////////////////////////


class DecayChain126
{


	public:
	DecayChain126(double sqrt_s_, double mA_,  double mB_,
					  double mA1_,    double mA2_, double mA3_,
					  double mB1_,    double mB2_, double mB3_);
	~DecayChain126();

	TwoBodyDecay   *DecayAB;
	ThreeBodyDecay *DecayA123;
	ThreeBodyDecay *DecayB123;

	void SetTag( const char tag_[] );
	void DisplayTag();
	void SetTotalThreeMomentum_pxpypz( double px, double py, double pz );
	void SetTotalThreeMomentum_pthetaphi(double p_, double theta, double phi_);

	void SetBitBoostBack(bool flag);
	double GetPhaseSpaceWeight(        double x1_,  double x2_,  double x3_,
					 						   double x4_,  double x5_,  double x6_,	
					 						   double x7_,  double x8_,  double x9_,	
					 						  double x10_, double x11_, double x12_);

	double GetPSConst();


	void SetPhaseSpace(  double x1,  double x2,  double x3,
					 						   double x4,  double x5,  double x6,	
					 						   double x7,  double x8,  double x9,	
					 						  double x10, double x11, double x12);

	void DisplayAllInfo();
	void DisplayKinematics();
	void DisplayMomenta();
	void DisplayConfiguration();
	void DisplayMasses();

	TLorentzVector *P;
	TLorentzVector *pA;
	TLorentzVector *pB;
	TLorentzVector *pA1;
	TLorentzVector *pA2;
	TLorentzVector *pA3;
	TLorentzVector *pB1;
	TLorentzVector *pB2;
	TLorentzVector *pB3;
	TLorentzVector *sum;
	TLorentzVector *sumA;
	TLorentzVector *sumB;

	private:

	bool BitBoostBack;
	TVector3 BetaVec;

	double PSWeight;
	double PSConst;

	double x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12;

	std::string tag;

	// System kinematics
	double E;
	double s;
	double sqrt_s;

	// Daughther particles
	double mA;
	double mB;
	double mA1;
	double mA2;
	double mA3;
	double mB1;
	double mB2;
	double mB3;

	double mA_sqr;
	double mB_sqr;
	double mA1_sqr;
	double mA2_sqr;
	double mA3_sqr;
	double mB1_sqr;
	double mB2_sqr;
	double mB3_sqr;
	
	// CoM Frame energy
	//double E1;
	//double E2;

	//double pmag;
	//double lambda;        
	//double sqrt_lambda;        

	//double costheta;
	//double sintheta;
	//double      phi;
	//double   cosphi;
	//double   sinphi;

};

#endif
