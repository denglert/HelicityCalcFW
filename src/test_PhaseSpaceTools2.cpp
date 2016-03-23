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
	std::cout << "###################################" << std::endl;
	std::cout << "## test_PhaseSpaceTools2 running ##" << std::endl;
	std::cout << "###################################" << std::endl;

	std::cout <<  std::endl;
	std::cout << "costheta = 2*x1 - 1" << std::endl;
	std::cout << "phi      = 2*pi*x2" << std::endl;

	double sqrt_s = 1.00;
	double m1 	  = 0.10;
	double m2 	  = 0.10;

	double x1; 
	double x2; 

	bool BoostBack = false;

	double px = 0.50;
	double py = 0.00;
	double pz = 0.00;


	TwoBodyDecay twobod( sqrt_s, m1, m2);
	twobod.SetTotalThreeMomentum_pxpypz( px, py, pz );
	twobod.SetBitBoostBack( BoostBack );

	x1 = 1.00;
	x2 = 0.00;

	printf("\n");
	printf("x1 = %.2f\n", x1);
	printf("x2 = %.2f\n", x2);

	twobod.SetPhaseSpace(x1,x2);
	twobod.DisplayMomenta();	

	x1 =  0.00;
	x2 =  0.00;

	printf("\n");
	printf("x1 = %.2f\n", x1);
	printf("x2 = %.2f\n", x2);

	twobod.SetPhaseSpace(x1,x2);
	twobod.DisplayMomenta();	

	x1 = 0.50;
	x2 = 0.25;

	printf("\n");
	printf("x1 = %.2f\n", x1);
	printf("x2 = %.2f\n", x2);

	twobod.SetPhaseSpace(x1,x2);
	twobod.DisplayMomenta();	

	x1 = 0.10;
	x2 = 0.55;
	printf("\n");
	printf("x1 = %.2f\n", x1);
	printf("x2 = %.2f\n", x2);

	twobod.SetPhaseSpace(x1,x2);
	twobod.DisplayMomenta();	

}
