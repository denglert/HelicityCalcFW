#ifndef HELICITYTOOLS_H
#define HELICITYTOOLS_H
#include "FortranInterface.h"
#include <complex>
#include <iomanip>
#include <iostream>

#define I compv(0.0,1.0)

typedef std::complex<double> compv;


//class HelicityTools
//{
//	public:
//	double rh_tautau_polarized(double *p1, double *p2, int &pol1, int &pol2);	double rh_tautau_unpolarized(double *p1, double *p2);	
//};


class TauMatrix
{
	private:	
	std::string name;
	double k0[4];
	int nParticles;
	double *p[4];
	double unpolarizedamp;

	public:
	compv m[2][2];
	int indeces;
	double CalcSumOfSquares();
	double CalcSumOfMagnitudes();
	void Show();
	void ShowSumOfSquares();
	void ShowUnpolarizedAmp();
	void Showk0();
	void Showpi(int i);
	void SetName(const char str[]);
	void ReadInCMatrix_2_2( CMatrix_2_2 &cmatrix );
	void Setk0(double k0_0, double k0_1, double k0_2, double k0_3);
	void SetnParts(int nParts);
	void Setpi(int i, double p_0, double p_1, double p_2, double p_3);
	void CalcUnpolarizedAmp();
};

#endif
