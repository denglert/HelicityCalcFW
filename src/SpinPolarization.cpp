#include <iostream>
#include <cmath>
#include <TH3D.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TGenPhaseSpace.h>
#include "HelicityFW.h"
#include "FortranInterface.h"
#include "PhysConst.h"
#include "UtilFunctions.h"
#include "PhaseSpaceTools.h"

int main( int argc, const char *argv[] )
{

	std::cout << "Testing spin polarization." << std::endl;

	double M = 5.00;
	double m1 = 1.0;
	double m2 = 1.0;

	double DaughterMasses[2] = { m1, m2 };
	TLorentzVector P (M, 0.0, 0.0, 0.0);

//	TGenPhaseSpace ScalarDecay;
//	ScalarDecay.SetDecay((*P), 2, DaughterMasses);
//
//	TLorentzVector *p1 = ScalarDecay.GetDecay(0);
//	TLorentzVector *p2 = ScalarDecay.GetDecay(1);

	double p1_3mag = sqrt( (M/2.0)*(M/2.0) - (m1*m1) );
	double p2_3mag = sqrt( (M/2.0)*(M/2.0) - (m2*m2) );

	TVector3 p1_3vec (0.0, 0.0,  p1_3mag);
	TVector3 p2_3vec (0.0, 0.0, -p1_3mag);

	TLorentzVector p1 ( p1_3vec, (M/2.0) );
	TLorentzVector p2 ( p2_3vec, (M/2.0) );

	TVector3 spinvec1 ( 0.0, 0.0, 1.0 );
	TVector3 spinvec2 ( 0.0, 0.0, 1.0 );


  	double theta1 = (0.0);
  	double phi1   = (0.0);
  	double theta2 = (0.0);
  	double phi2   = (0.0);
//
//	double theta1 = (M_PI/2.0);
//	double phi1   = (0.0);
//	double theta2 = (M_PI/2.0);
//	double phi2   = (M_PI);

	spinvec1.SetTheta(theta1);
	spinvec1.SetPhi(phi1);

	spinvec2.SetTheta(theta2);
	spinvec2.SetPhi(phi2);

	TLorentzVector polvec1( spinvec1, 0.0);
	TLorentzVector polvec2( spinvec2, 0.0);

	TVector3 BoostVec1 = p1.BoostVector();
	TVector3 BoostVec2 = p2.BoostVector();

	polvec1.Boost( BoostVec1 );
	polvec2.Boost( BoostVec2 );

	std::cout << "p1: " << std::endl;
	displayTLorentzVector(&p1);
	std::cout << "p2: " << std::endl;
	displayTLorentzVector(&p2);

	std::cout << "polvec1: " << std::endl;
	displayTLorentzVector(&polvec1);
	std::cout << "polvec2: " << std::endl;
	displayTLorentzVector(&polvec2);

	std::cout <<  std::endl;
	std::cout << "Consistency checks: " << std::endl;

	std::cout << "polvec1*polvec1: " << polvec1*polvec1 << std::endl;
	std::cout << "polvec2*polvec2: " << polvec2*polvec2 << std::endl;

	std::cout << "p1*polvec1: " << p1*polvec1 << std::endl;
	std::cout << "p2*polvec2: " << p2*polvec2 << std::endl;

	std::cout <<  std::endl;

	std::cout << "polvec1*polvec2: " << polvec1*polvec2 << std::endl;
	std::cout << "p1*polvec2: " << p1*polvec2 << std::endl;
	std::cout << "p2*polvec1: " << p2*polvec1 << std::endl;

	std::cout << "Result" << std::endl;

	double result = ((p1*p2) - m1*m1)*( 1 - (polvec1*polvec2)) - (p1*polvec2)*(p2*polvec1);
	std::cout << "[ (p1p2)-m^{2} ] [ 1 - (s1s2) ] - (p1s2)(p2s1)" << std::endl;
	std::cout << result << std::endl;


}
