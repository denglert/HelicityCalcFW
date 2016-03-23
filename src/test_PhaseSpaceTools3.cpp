#include <iostream>
#include <cmath>
#include <TCanvas.h>
#include "HelicityFW.h"
#include "FortranInterface.h"
#include "PhysConst.h"
#include "UtilFunctions.h"
#include "PhaseSpaceTools.h"

int main( int argc, const char *argv[] )
{

	std::cout <<  std::endl;
	std::cout << "##################################" << std::endl;
	std::cout << "### --- DecayChain126 test --- ###" << std::endl;
	std::cout << "##################################" << std::endl;
	std::cout <<  std::endl;


	const double mA     = 5.00;
	const double mB     = 5.00;
	const double mA1    = 1.00;
	const double mA2    = 0.00;
	const double mA3    = 0.00;
	const double mB1    = 1.00;
	const double mB2    = 0.00;
	const double mB3    = 0.00;

	const double x1   = 0.500; //    costheta = 2*x1 - 1
	const double x2   = 0.250; //         phi = 2*pi*x2
	const double x3   = 0.9999; //        sA23 = x3 * sA23_length + sA23_min
	const double x4   = 1.000; //   costhetaA = 2*x4 - 1
	const double x5   = 0.000; //        phiA = 2*pi*x5
	const double x6   = 0.000; // costhetaA23 = 2*x6 - 1
	const double x7   = 0.000; //      phiA23 = 2*pi*x7
	const double x8   = 0.050; //        sB23 = x8 * sB23_length + sB23_min
	const double x9   = 0.000; //   costhetaB = 2*x9 - 1
	const double x10  = 0.000; //        phiB = 2*pi*x10
	const double x11  = 0.000; // costhetaB23 = 2*x11 - 1
	const double x12  = 0.000; //      phiB23 = 2*pi*x12

	const double M = 100.00;
	const double Px = 5.00;
	const double Py = 0.00;
	const double Pz = 0.00;
	const double E  = sqrt( Px*Px + Py*Py + Pz*Pz + M*M );

	const bool BitBoostBack = true;

	// Initialize DecayChain126
	DecayChain126 decay(M, mA,  mB,
						 	     mA1, mA2, mA3,
								  mB1, mB2, mB3);

	decay.SetTag("DecayChain126");

	TLorentzVector *P = decay.P;
	P->SetPxPyPzE(Px,Py,Pz,E);
	decay.SetBitBoostBack( BitBoostBack );

	// Pointers to momenta
	TLorentzVector *pA;
	TLorentzVector *pB;
	TLorentzVector *pA1;
	TLorentzVector *pA2;
	TLorentzVector *pA3;
	TLorentzVector *pB1;
	TLorentzVector *pB2;
	TLorentzVector *pB3;
	
	// Link pointers
	pA  = decay.pA;
	pB  = decay.pB;
	pA1 = decay.pA1;
	pA2 = decay.pA2;
	pA3 = decay.pA3;
	pB1 = decay.pB1;
	pB2 = decay.pB2;
	pB3 = decay.pB3;

	// Set phase space to config {x1, x2, ..., x12}
	decay.SetPhaseSpace( x1, x2,
						 	   x3, x4,  x5,  x6,  x7,
								x8, x9, x10, x11, x12);

	// Display information on std output
	decay.DisplayAllInfo();

}
