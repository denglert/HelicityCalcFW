#ifndef PHASESPACETOOLS_H
#define PHASESPACETOOLS_H

#include <cmath>
#include <iostream>
#include "UtilFunctions.h"
#include <TLorentzVector.h>

////////////////////
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
	ThreeBodyDecay();
	ThreeBodyDecay(double M, double m1, double m2, double m3);
	~ThreeBodyDecay(){};

	void SetMotherMass(double mass);
	void SetDecayMass (int i, double mass);
	void SetMotherp(double E, double px, double py, double pz);
	double GetPhaseSpaceWeight(double x1, double x2, double x3, double x4, double x5);
	void SetPhaseSpace(double x1, double x2, double x3, double x4, double x5);

	void DisplayAll();
	void DisplayKinematics();
	void DisplayMomenta();
	void DisplayMasses();

	private:
	double M;
	double m[3];
	double M_sqr;
	double m_sqr[3];

	double Jacobian;
	double PSFactor;
	double PSWeight;

	double    betabar;        
	double  betabar23;        
	double   s23_sqrt;        
	double        s23;        
	double  costheta1;  
	double  sintheta1;  
	double       phi1;  
	double    cosphi1;  
	double    sinphi1;  
	double costheta23; 
	double sintheta23; 
	double      phi23;
	double   cosphi23;  
	double   sinphi23;  

	// In the rest frame of the mother particle
	double E1;
	double E23;
	double pmag;

	// In the rest frame of 23
	double p2mag_23;
	double E2_23;
	double E3_23;

	TLorentzVector P;
	TLorentzVector p[3];
	TLorentzVector sum;
	TVector3 Beta23;

};



#endif
