#include "PhaseSpaceTools.h"

////////////////////
// Two-body decay //
namespace TwoBodyFunc
{

	double beta(double s, double m1_sqr, double m2_sqr)
	{
		double beta;
		beta = sqrt( 1 - 2*(m1_sqr + m2_sqr)/s + (m1_sqr-m2_sqr)*(m1_sqr-m2_sqr)/s/s );
		return beta;
	};


	double p (double s, double m1_sqr, double m2_sqr)
	{
		double p;
		p = sqrt(s)*TwoBodyFunc::beta(s,m1_sqr,m2_sqr)/2;
		return p;
	};


	double E (double s, double m_sqr, double m_other_sqr)
	{
		double E;
		E = sqrt(s)*( 1 + m_sqr/s - m_other_sqr/s )/2;
		return E;
	}

}

//////////////////////////
// class ThreeBodyDecay //
ThreeBodyDecay::ThreeBodyDecay()
{
	Jacobian = 1;
}

ThreeBodyDecay::ThreeBodyDecay(double M, double m1, double m2, double m3)
{
	Jacobian = 1;

	M = M;
	m[0] = m1;
	m[1] = m2;
	m[2] = m3;
}

void ThreeBodyDecay::SetDecayMass (int i, double mass)
{
	m[i]     = mass;
	m_sqr[i] = mass*mass;
};

void ThreeBodyDecay::SetMotherMass (double mass)
{
	M     = mass;
	M_sqr = mass*mass;
};


void ThreeBodyDecay::SetPhaseSpace(double x1, double x2, double x3, double x4, double x5)
{
	s23_sqrt  = x1*(M-m[0]-m[1]-m[2]) + m[1] + m[2];
		  s23  = s23_sqrt*s23_sqrt;
	costheta1 = 2*x2-1;
	sintheta1 = sqrt(1-costheta1*costheta1);
	     phi1 = 2*M_PI*x3;
	  cosphi1 = cos(phi1);
	  sinphi1 = sin(phi1);

	costheta23 = 2*x4-1;
	sintheta23 = sqrt(1-costheta23*costheta23);
	     phi23 = 2*M_PI*x5;
	  cosphi23 = cos(phi23);
	  sinphi23 = sin(phi23);

	// In the rest frame of the mother particle
	betabar = TwoBodyFunc::beta ( M_sqr, m_sqr[0], s23);
	pmag    = M*betabar/2;
	E1   = TwoBodyFunc::E ( M_sqr, m_sqr[0], s23);
	E23  = M-E1;

	// In the rest frame of 23
	betabar23 = TwoBodyFunc::p    (s23, m_sqr[1], m_sqr[2]);
	p2mag_23  = s23_sqrt*betabar/2;
		E2_23 = TwoBodyFunc::E    (s23, m_sqr[1], m_sqr[2]);
		E3_23 = s23_sqrt - E2_23;

	p[0].SetPxPyPzE(pmag*sintheta1*cosphi1,pmag*sintheta1*sinphi1,pmag*costheta1,E1);
	Beta23 = -p[0].Vect();
	Beta23.SetMag(Beta23.Mag()/E23);

	p[1].SetPxPyPzE(p2mag_23*sintheta23*cosphi23,p2mag_23*sintheta1*sinphi23,p2mag_23*costheta23,E2_23);
	p[2].SetPxPyPzE(-(p[1])[0],-(p[1])[1],-(p[1])[2],E3_23);

	p[1].Boost(Beta23);
	p[2].Boost(Beta23);

	PSFactor = betabar*betabar23;
	PSWeight = PSFactor/Jacobian;
}


double ThreeBodyDecay::GetPhaseSpaceWeight(double x1, double x2, double x3, double x4, double x5)
{
	ThreeBodyDecay::SetPhaseSpace(x1,x2,x3,x4,x5);
	return PSWeight;
}

void ThreeBodyDecay::DisplayAll()
{
	std::cout << "Displaying all info..." << std::endl;
	ThreeBodyDecay::DisplayMasses();
	ThreeBodyDecay::DisplayMomenta();
	ThreeBodyDecay::DisplayKinematics();
}


void ThreeBodyDecay::DisplayKinematics()
{
	std::cout << "betabar: " << betabar << std::endl;
	std::cout << "pmag: " << pmag << std::endl;
	std::cout << "E1: "   << E1 << std::endl;
	std::cout << "E23: "   << E1 << std::endl;
	std::cout << "s23_sqrt: " << s23_sqrt << std::endl;
	std::cout << "betabar23: " << betabar23 << std::endl;
	std::cout << "p2mag_23: " << p2mag_23 << std::endl;
	std::cout << "E2_23 " << E2_23 << std::endl;
	std::cout << "E3_23 " << E3_23 << std::endl;
	std::cout << "PSFactor " << PSFactor << std::endl;
	std::cout << "Jacobian " << Jacobian << std::endl;
	std::cout << "PSWeight " << PSWeight << std::endl;
}

void ThreeBodyDecay::DisplayMasses()
{
	std::cout << "M: " << M << std::endl;
	std::cout << "m[0]: " << m[0] << std::endl;
	std::cout << "m[1]: " << m[1] << std::endl;
	std::cout << "m[2]: " << m[2] << std::endl;
}

void ThreeBodyDecay::DisplayMomenta()
{
	sum = p[0] + p[1] + p[2];
	std::cout << "In the mother particle rest frame:" << std::endl;

	std::cout << "p1+p2+p3:" << std::endl;
	displayTLorentzVector(&sum);

	std::cout << "p1:" << std::endl;
	displayTLorentzVector(&p[0]);

	std::cout << "p2:" << std::endl;
	displayTLorentzVector(&p[1]);

	std::cout << "p3:" << std::endl;
	displayTLorentzVector(&p[2]);
}
