#include "HelicityTools.h"

//double HelicityTools::rh_tautau_polarized(double *p1, double *p2, int &pol1, int &pol2)
//{
//	double res = rh_tautau_(p1,p2,&pol1,&pol2);
//};
//
//double HelicityTools::rh_tautau_unpolarized(double *p1, double *p2)
//{
//	double res = 0;
//	for(int pol1=1; pol1 < 3; pol1++)	
//	for(int pol2=1; pol2 < 3; pol2++)	
//	{ res = res + rh_tautau_(p1,p2,&pol1,&pol2); }
//	return res;
//};


/////////////////////////////
// --- TauMatrix class --- //
/////////////////////////////

double TauMatrix::CalcSum()
{
	double sum = 0;
	for (int i=0; i<2; i++)
	for (int j=0; j<2; j++)
	{
		sum = sum + m[i][j].real()*m[i][j].real()+m[i][j].imag()*m[i][j].imag();
	}
	return sum;
}

void TauMatrix::Display()
{
	std::cout << name << " matrix: " << std::endl;
	for (int i=0; i<2; i++)
	{
		for (int j=0; j<2; j++)
		{
			std::cout << std::setw(12) << std::setprecision(4) << m[i][j] << " " ;
		}
		std::cout << std::endl;
	}
}

void TauMatrix::DisplaySum()
{
	std::cout << name << " sum of squared matrix elements: " << std::setw(12) << std::setprecision(4) << TauMatrix::CalcSum() << std::endl ;
}


void TauMatrix::ReadInCMatrix_2_2( CMatrix_2_2 &cmatrix )
{

	for (int i=0; i<2; i++)
	{
		for (int j=0; j<2; j++)
		{
			m[i][j] = cmatrix[j][i].r + I*cmatrix[j][i].i;
		}
	}

}


void TauMatrix::SetName(const char str[])
{
	name = str;
}
