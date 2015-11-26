#ifndef HELICITYTOOLS_H
#define HELICITYTOOLS_H
#include "FortranInterface.h"
#include <complex>
#include <iomanip>
#include <iostream>

#define I compv(0.0,1.0)

typedef std::complex<double> compv;


class HelicityTools
{
	public:
	double rh_tautau_polarized(double *p1, double *p2, int &pol1, int &pol2);	double rh_tautau_unpolarized(double *p1, double *p2);	
};


class TauMatrix
{
	
	public:
	compv m[2][2];
	int indeces;
	double CalcSum();
	void Display();
	void ReadInCMatrix_2_2( CMatrix_2_2 &cmatrix );

};

#endif
