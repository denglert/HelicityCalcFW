#include <iostream>
#include "HelicityTools.h"

int main( int argc, const char *argv[] )
{
	std::cout << "HelicityTools test program" << std::endl;

	TauMatrix tau;

	tau.m[0][0] = 0+2I;
	tau.m[1][1] = 1+0I;

	std::cout << "tau.m[0][0].real(): " << tau.m[0][0].real() << std::endl;
	std::cout << "tau.m[0][0].real(): " << tau.m[0][0].imag() << std::endl;
	std::cout << "tau.m[1][1].real(): " << tau.m[1][1].real() << std::endl;
	std::cout << "tau.m[1][1].real(): " << tau.m[1][1].imag() << std::endl;

	std::cout << "Sum of the matrix elements squared: " << tau.CalcSumOfSquares() <<
			  std::endl;

	tau.Show();

}
