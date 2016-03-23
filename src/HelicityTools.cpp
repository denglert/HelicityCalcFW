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

void TauMatrix::Show()
{
	std::cout << std::endl << name << " matrix: " << std::endl;
	for (int i=0; i<2; i++)
	{
		for (int j=0; j<2; j++)
		{
			std::cout << std::setw(12) << std::setprecision(4) << m[i][j] << " " ;
		}
		std::cout << std::endl;
	}

	std::cout << std::endl;
};

void TauMatrix::ShowSum()
{
	std::cout << name << " sum of squared matrix elements: " << std::setw(12) << std::setprecision(4) << TauMatrix::CalcSum() << std::endl ;
};


void TauMatrix::ReadInCMatrix_2_2( CMatrix_2_2 &cmatrix )
{

	for (int i=0; i<2; i++)
	{
		for (int j=0; j<2; j++)
		{
			m[i][j] = cmatrix[j][i].r + I*cmatrix[j][i].i;
		}
	}

};


void TauMatrix::SetName(const char str[])
{
	name = str;
};


void TauMatrix::Setk0(double k0_0, double k0_1, double k0_2, double k0_3)
{
	k0[0] = k0_0;	
	k0[1] = k0_1;	
	k0[2] = k0_2;	
	k0[3] = k0_3;	
};

void TauMatrix::SetnParts(int nParts)
{
	nParticles=nParts;
	for(int mu=0; mu<4; mu++)
	{
		p[mu] = new double[nParticles];
	}
};


void TauMatrix::Setpi(int i, double p_0, double p_1, double p_2, double p_3)
{
	p[0][i] = p_0;
	p[1][i] = p_1;
	p[2][i] = p_2;
	p[3][i] = p_3;
}


void TauMatrix::CalcUnpolarizedAmp()
{
	double sum = TauMatrix::CalcSum();
	double factor = 1.0;

	for(int iPart = 0; iPart < nParticles; iPart++)
	{
		factor = factor*(p[0][iPart]*k0[0]-p[1][iPart]*k0[1]-p[2][iPart]*k0[2]-p[3][iPart]*k0[3]);
	}

	unpolarizedamp = sum/factor;

}

void TauMatrix::ShowUnpolarizedAmp()
{
	std::cout << name << " unpolarized amplitude: " << std::setw(12) << std::setprecision(7) << unpolarizedamp << std::endl ;
};

void TauMatrix::Showk0()
{
	
	for(int mu=0; mu<4; mu++)
	{
		std::cout << "k0[" << mu << "] " << std::setprecision(5) << k0[mu] << std::endl ;
	}
}

void TauMatrix::Showpi(int i)
{
	
	for(int mu=0; mu<4; mu++)
	{
		std::cout << "p[" << mu << "][" << i << "] " << std::setprecision(5) << p[mu][i] << std::endl ;
	}
}
