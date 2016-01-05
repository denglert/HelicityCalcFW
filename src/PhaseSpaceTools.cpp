#include "PhaseSpaceTools.h"

////////////////////
// Two-body decay //
namespace TwoBodyFunc
{

	double beta(double s, double m1_sqr, double m2_sqr)
	{
		double beta;
		beta = sqrt( 1.0 - 2.0*(m1_sqr + m2_sqr)/s + ((m1_sqr-m2_sqr)*(m1_sqr-m2_sqr)/s/s) );
		return beta;
	};


	double lambda(double a, double b, double c)
	{
		double value;
		value = (a-b-c)*(a-b-c)-4*b*c;
		return value;
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
		E = sqrt(s)*( 1.0 + m_sqr/s - m_other_sqr/s )/2.0;
		return E;
	}

}

ThreeBodyDecay::ThreeBodyDecay(double M_, double m1_, double m2_, double m3_)
{

	P = new TLorentzVector();
	p = new TLorentzVector*[3];
	p[0] = new TLorentzVector();
	p[1] = new TLorentzVector();
	p[2] = new TLorentzVector();

	M = M_;
	m[0] = m1_;
	m[1] = m2_;
	m[2] = m3_;

	M_sqr = M_*M_;
	m_sqr[0] = m1_*m1_;
	m_sqr[1] = m2_*m2_;
	m_sqr[2] = m3_*m3_;

	s23_min    = (m[1]+m[2])*(m[1]+m[2]);
	s23_length = (M-m[0])*(M-m[0]) - s23_min;
	PSConst = s23_length/pow(M_PI,3.0)/128.0;

	BitBoostBack = false;
}

ThreeBodyDecay::~ThreeBodyDecay()
{
	delete P;
	delete p[0];
	delete p[1];
	delete p[2];
}


void ThreeBodyDecay::SetBitBoostBack(bool flag)
{
	BitBoostBack = flag;
}

// OBSOLETE
//void ThreeBodyDecay::SetDecayMass (int i, double mass)
//{
//	m[i]     = mass;
//	m_sqr[i] = mass*mass;
//};
//
//void ThreeBodyDecay::SetMotherMass (double mass)
//{
//	M     = mass;
//	M_sqr = mass*mass;
//};


void ThreeBodyDecay::SetMotherEPxPyPz(double E, double px, double py, double pz)
{
	P->SetPxPyPzE(px,py,pz,E);
	BetaVec = P->BoostVector();

	# if DEBUG
	printf("\n\nSTART DEBUG INFO\n");
	printf("BetaVec:\n");
	printf("bx = %12.6e\n", BetaVec.x());
	printf("by = %12.6e\n", BetaVec.y());
	printf("bz = %12.6e\n", BetaVec.z());
	printf("END DEBUG INFO\n\n");
	# endif
}

void ThreeBodyDecay::SetMotherMPThetaPhi(double M, double mom, double theta, double phi)
{

	P->SetPxPyPzE(mom,0.0,0.0,sqrt(mom*mom+M*M));
	P->SetTheta(theta);
	P->SetPhi(phi);
	BetaVec = P->BoostVector();

	# if DEBUG
	printf("\n\nSTART DEBUG INFO\n");
	printf("BetaVec:\n");
	printf("bx = %12.6e\n", BetaVec.x());
	printf("by = %12.6e\n", BetaVec.y());
	printf("bz = %12.6e\n", BetaVec.z());
	printf("END DEBUG INFO\n\n");
	# endif
}


void ThreeBodyDecay::SetPhaseSpace(double x1, double x2, double x3, double x4, double x5)
{

		  s23  = s23_min + s23_length*x1;
   s23_sqrt  = sqrt(s23);
	costheta1 = 1-2*x2;
	sintheta1 = sqrt(1-costheta1*costheta1);
	     phi1 = 2*M_PI*x3;
	  cosphi1 = cos(phi1);
	  sinphi1 = sin(phi1);

	costheta23 = 1-2*x4;
	sintheta23 = sqrt(1-costheta23*costheta23);
	     phi23 = 2*M_PI*x5;
	  cosphi23 = cos(phi23);
	  sinphi23 = sin(phi23);

	// In the rest frame of the mother particle
	betabar = TwoBodyFunc::beta ( M_sqr, m_sqr[0], s23);
	lambda  = TwoBodyFunc::lambda ( M_sqr, m_sqr[0], s23);
	lambda_sqrt  = sqrt(lambda);
	pmag    = M*betabar/2;
	E1   = TwoBodyFunc::E ( M_sqr, m_sqr[0], s23);
	E23  = M-E1;

	// In the rest frame of 23
	betabar23 		= TwoBodyFunc::beta   (s23, m_sqr[1], m_sqr[2]);
	// betabar23 		= 1.0;
	lambda23  		= TwoBodyFunc::lambda (s23, m_sqr[1], m_sqr[2]);
	lambda23_sqrt  = sqrt(lambda23);
	p2mag_23  = s23_sqrt*betabar23/2;
		E2_23 = TwoBodyFunc::E    (s23, m_sqr[1], m_sqr[2]);
		E3_23 = s23_sqrt - E2_23;

	p[0]->SetPxPyPzE(pmag*sintheta1*cosphi1,pmag*sintheta1*sinphi1,pmag*costheta1,E1);
	Beta23 = -p[0]->Vect();
	Beta23.SetMag(Beta23.Mag()/E23);

	p[1]->SetPxPyPzE(p2mag_23*sintheta23*cosphi23, p2mag_23*sintheta23*sinphi23, p2mag_23*costheta23, E2_23);
	p[2]->SetPxPyPzE(-(*p[1])[0], -(*p[1])[1], -(*p[1])[2], E3_23);

	# if DEBUG
	std::cerr << "\n\nSTART DEBUG INFO.\n";

	printf("--- INPUT ---\n" );
	printf("x1:   %10.5f\n", x1 );
	printf("x2:   %10.5f\n", x2 );
	printf("x3:   %10.5f\n", x3 );
	printf("x4:   %10.5f\n", x4 );
	printf("x5:   %10.5f\n", x5 );

	printf("M:    %10.5f\n", M );
	printf("m[0]: %10.5f\n", m[0] );
	printf("m[1]: %10.5f\n", m[1] );
	printf("m[2]: %10.5f\n", m[2] );

	printf("M_sqr:    %10.5f\n", M_sqr );
	printf("m_sqr[0]: %10.5f\n", m_sqr[0] );
	printf("m_sqr[1]: %10.5f\n", m_sqr[1] );
	printf("m_sqr[2]: %10.5f\n", m_sqr[2] );

	printf("--- Physical variables ---\n" );
	printf("s23_sqrt:     %10.5f\n", s23_sqrt );
	printf("s23:          %10.5f\n", s23 );
	printf("cos(theta1):  %10.5f\n", costheta1 );
	printf("phi1:         %10.5f\n", (phi1) );
	printf("cos(theta23): %10.5f\n", costheta23 );
	printf("phi23:        %10.5f\n", (phi23) );

	printf("--- Calculated variables ---\n" );
	printf("E1 = (M/2)*(1 + (m1_sqr/M_sqr) - (s23_sqr/M_sqr)) = %10.5f\n", E1 );
	printf("betabar               %10.5f\n", betabar );
	printf("lambda_sqrt/M/M    %10.5f\n", (lambda_sqrt/M/M));
	printf("pmag = (M/2)*betabar  %10.5f\n", pmag );
	printf("E2_23: %10.5f\n", E2_23 );
	printf("E3_23: %10.5f\n", E3_23 );
	printf("Beta23:\n");
	printf("b23x = %12.6e\n", Beta23.x());
	printf("b23y = %12.6e\n", Beta23.y());
	printf("b23z = %12.6e\n", Beta23.z());


	std::cerr << "--- 2-3 rest frame ---\n";
	printf("betabar23               %10.5f\n", betabar23 );
	printf("lambda23_sqrt/s23       %10.5f\n", (lambda23_sqrt/s23) );

	std::cerr << "p[1]:";
	displayTLorentzVector(p[1]);
	std::cerr << "p[2]:";
	displayTLorentzVector(p[2]);

	printf("p[1]->CosTheta(): %10.5f\n", p[1]->CosTheta() );
	printf("p[1]->Phi():      %10.5f\n", p[1]->Phi() );
	printf("p[2]->CosTheta(): %10.5f\n", p[2]->CosTheta() );
	printf("p[2]->Phi():      %10.5f\n", p[2]->Phi() );

	std::cerr << "Their masses:\n";
	printf("p[1]->M(): %10.5f\n", p[1]->M() );
	printf("p[2]->M(): %10.5f\n", p[2]->M() );

	# endif

	p[1]->Boost(Beta23);
	p[2]->Boost(Beta23);

	PSWeight = betabar*betabar23;

	# if DEBUG
	std::cerr << "-- M rest frame --\n";
	printf("p[0]->M(): %10.5f\n", p[0]->M() );
	printf("p[1]->M(): %10.5f\n", p[1]->M() );
	printf("p[2]->M(): %10.5f\n", p[2]->M() );

	printf("P:");
	displayTLorentzVector(P);
	printf("p[0]:");
	displayTLorentzVector(p[0]);
	printf("p[1]:");
	displayTLorentzVector(p[1]);
	printf("p[2]:");
	displayTLorentzVector(p[2]);

	printf("p[0]->CosTheta()): %10.5f\n", p[0]->CosTheta() );
	printf("p[0]->Phi()):      %10.5f\n", p[0]->Phi() );
	printf("p[1]->CosTheta()): %10.5f\n", p[1]->CosTheta() );
	printf("p[1]->Phi()):      %10.5f\n", p[1]->Phi() );
	printf("p[2]->CosTheta()): %10.5f\n", p[2]->CosTheta() );
	printf("p[2]->Phi()):      %10.5f\n", p[2]->Phi() );

	printf("PSWeight = betabar*betabar23 = %10.5f\n", PSWeight);
	std::cerr << "END DEBUG INFO.\n\n";
	# endif

	if( BitBoostBack )
	{
	# if DEBUG
	printf("BetaVec:\n");
	printf("bx = %12.6e\n", BetaVec.x());
	printf("by = %12.6e\n", BetaVec.y());
	printf("bz = %12.6e\n", BetaVec.z());
	printf("END DEBUG INFO\n\n");
	# endif

	p[0]->Boost(BetaVec);
	p[1]->Boost(BetaVec);
	p[2]->Boost(BetaVec);


	}



}


double ThreeBodyDecay::GetJacobian()
{
	return Jacobian;
}

double ThreeBodyDecay::GetPSConst()
{
	return PSConst;
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
	sum = (*p[1]) + (*p[2]);
	std::cout << "(*p[1]) + (*p[2]).M(): " << sum.M() << std::endl;
}

void ThreeBodyDecay::DisplayMasses()
{
	printf("M: %10.5f\n", M) ;
	printf("m[0]: %10.5f\n", m[0] );
	printf("m[1]: %10.5f\n", m[1] );
	printf("m[2]: %10.5f\n", m[2] );
	printf("P->M(): %10.5f\n", P->M());
	printf("p[0]->M(): %10.5f\n", p[0]->M());
	printf("p[1]->M(): %10.5f\n", p[1]->M());
	printf("p[2]->M(): %10.5f\n", p[2]->M());
}

void ThreeBodyDecay::DisplayMomenta()
{
	sum = (*p[0]) + (*p[1]) + (*p[2]);
	std::cout << "In the mother particle rest frame:" << std::endl;

	std::cout << "p1+p2+p3:" << std::endl;
	displayTLorentzVector(&sum);

	std::cout << "p1:" << std::endl;
	displayTLorentzVector(p[0]);

	std::cout << "p2:" << std::endl;
	displayTLorentzVector(p[1]);

	std::cout << "p3:" << std::endl;
	displayTLorentzVector(p[2]);
}
