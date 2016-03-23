#include "PhaseSpaceTools.h"

///////////////////////////////////
// --- TwoBodyFunc namespace --- //
///////////////////////////////////

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
		p = sqrt(s)*TwoBodyFunc::beta(s,m1_sqr,m2_sqr)/2.0;
		return p;
	};


	double E (double s, double m_sqr, double m_other_sqr)
	{
		double E;
		double sqrt_s = sqrt(s);
		E = sqrt_s*( 1 + (m_sqr/s) - (m_other_sqr/s))/2.0;
		return E;
	}

}

////////////////////////////////
// --- TwoBodyDecay class --- //
////////////////////////////////
// Note: The BetaVec boost vector is calculated when someone calls
// 		the SetBitBoostBack() function

/////////////////////////////////////////
TwoBodyDecay::TwoBodyDecay( double sqrt_s_, double m1_, double m2_ )
{
	P   = new TLorentzVector();
	p1  = new TLorentzVector();
	p2  = new TLorentzVector();
	sum = new TLorentzVector();

	sqrt_s  = sqrt_s_;
	m1 	  = m1_;
	m2 	  = m2_;

	s      = sqrt_s*sqrt_s;
	m1_sqr = m1_*m1_;
	m2_sqr = m2_*m2_;

	BitBoostBack = false;

	lambda      = TwoBodyFunc::lambda ( s, m1_sqr, m2_sqr );
	sqrt_lambda = sqrt(lambda);
	pmag        = sqrt_lambda/2.0/sqrt_s;

	E1 = TwoBodyFunc::E ( s, m1, m2);
	E2 = sqrt_s - E1;

	PSWeight = sqrt_lambda/s/8.0/M_PI;
}


/////////////////////////////////////////
TwoBodyDecay::~TwoBodyDecay()
{
	delete P;
	delete p1;
	delete p2;
	delete sum;
}

///////////////////////////////////////////////////////////////////////
void TwoBodyDecay::DisplayTag()
{
	printf("\n### %s ###\n", tag.c_str());
}

/////////////////////////////////////////
void TwoBodyDecay::SetTag(const char tag_[])
{
	tag = tag_;
}

/////////////////////////////////////////
void TwoBodyDecay::SetTotalThreeMomentum_pxpypz( double px, double py, double pz )
{
	// Total energy of the (boosted system)/(decaying mother particle) viewed from the LAB frame
	E = sqrt( s + px*px + py*py + pz*pz );
	P->SetPxPyPzE(px,py,pz,E);
	BetaVec = P->BoostVector();

	# if DEBUG
	printf("\n\nSTART DEBUG INFO\n");
	printf("BetaVec components:\n");
	printf("bx = %12.6e\n", BetaVec.x());
	printf("by = %12.6e\n", BetaVec.y());
	printf("bz = %12.6e\n", BetaVec.z());
	printf("END DEBUG INFO\n\n");
	# endif
}


/////////////////////////////////////////
void TwoBodyDecay::SetTotalThreeMomentum_pthetaphi(double p_, double theta_, double phi_)
{
	E = sqrt( s + p_*p_ );
	P->SetPxPyPzE(p_,0.0,0.0,E);
	P->SetTheta(theta_);
	P->SetPhi(phi_);
	BetaVec = P->BoostVector();

	# if DEBUG
	printf("\n\nSTART DEBUG INFO\n");
	printf("TwoBodyDecay::SetMotherMPThetaPhi(double mom, double theta, double phi)\n");
	printf("BetaVec:\n");
	printf("bx = %12.6e\n", BetaVec.x());
	printf("by = %12.6e\n", BetaVec.y());
	printf("bz = %12.6e\n", BetaVec.z());
	printf("END DEBUG INFO\n\n");
	# endif
}


/////////////////////////////////////////
void TwoBodyDecay::SetPhaseSpace( double x1_, double x2_ )
{

	x1 = x1_;
	x2 = x2_;

	costheta = 2*x1-1;
	sintheta = sqrt(1-costheta*costheta);
	     phi = 2*M_PI*x2;
	  cosphi = cos(phi);
	  sinphi = sqrt(1-cosphi*cosphi);

	p1->SetPxPyPzE( pmag*sintheta*cosphi, pmag*sintheta*sinphi, pmag*costheta, E1);
	p2->SetPxPyPzE(            -(*p1)[0],            -(*p1)[1],     -(*p1)[2], E2);

	if( BitBoostBack )
	{

		# if DEBUG
		printf("BetaVec:\n");
		printf("bx = %12.6e\n", BetaVec.x());
		printf("by = %12.6e\n", BetaVec.y());
		printf("bz = %12.6e\n", BetaVec.z());
		printf("END DEBUG INFO\n\n");
		# endif

		p1->Boost(BetaVec);
		p2->Boost(BetaVec);

	}

}


/////////////////////////////////////////
double TwoBodyDecay::GetPhaseSpaceWeight()
{
	return PSWeight;
}


/////////////////////////////////////////
void TwoBodyDecay::SetBitBoostBack(bool flag)
{
	BetaVec = P->BoostVector();
	BitBoostBack = flag;
}


// Display functions

void TwoBodyDecay::DisplayAll()
{
	TwoBodyDecay::DisplayMasses();
	TwoBodyDecay::DisplayMomenta();
	TwoBodyDecay::DisplayKinematics();
}


/////////////////////////////////////////
void TwoBodyDecay::DisplayConfiguration()
{
	TwoBodyDecay::DisplayTag();
	printf("x1 = %5.2f <->  cos(theta) = %5.2f \n", x1, costheta);
	printf("x2 = %5.2f <->         phi = %5.2f \n", x2, phi);
}

/////////////////////////////////////////
void TwoBodyDecay::DisplayMasses()
{
	printf("\n-- Mass Information --\n") ;
	printf("sqrt_s:  %10.5f\n", sqrt_s) ;
	printf("m1:      %10.5f\n", m1 );
	printf("m2:      %10.5f\n", m2 );
	printf("P->M():  %10.5f\n", P->M());
	printf("p1->M(): %10.5f\n", p1->M());
	printf("p2->M(): %10.5f\n", p2->M());
}


/////////////////////////////////////////
void TwoBodyDecay::DisplayMomenta()
{

	(*sum) = (*p1) + (*p2);

	printf("\n-- Momenta Information --\n") ;
	std::cout << "P:" << std::endl;
	displayTLorentzVector(P);

	std::cout << "p1+p2:" << std::endl;
	displayTLorentzVector(sum);

	std::cout << "p1:" << std::endl;
	displayTLorentzVector(p1);

	std::cout << "p2:" << std::endl;
	displayTLorentzVector(p2);

}


/////////////////////////////////////////
void TwoBodyDecay::DisplayKinematics()
{
	// System kinematics
	printf("\n-- Kinematics Information --\n") ;
	printf("E:           %.5f\n", E);
	printf("s:           %.5f\n", s);
	printf("sqrts_s:     %.5f\n", sqrt_s);
	printf("E1:          %.5f\n", E1);
	printf("E2:          %.5f\n", E2);
	printf("pmag:        %.5f\n", pmag);
	printf("lambda:      %.5f\n", lambda);        
	printf("sqrt_lambda: %.5f\n", sqrt_lambda);        

	printf("\n");
	printf("costheta:   %.5f\n",costheta);
	printf("sintheta:   %.5f\n",sintheta);
	printf("phi:        %.5f\n",     phi);
	printf("cosphi:     %.5f\n",  cosphi);
	printf("sinphi:     %.5f\n",  sinphi);

	printf("\n");
	printf("1/(8*PI):   %.4f\n", 1.0/(8.0*M_PI)) ;
	printf("PSWeight = sqrt_lambda/(s*8*PI)\n") ;
	printf("PSWeight:   %.4f\n", PSWeight) ;

}


//////////////////////////////////////
//// --- ThreeBodyDecay class --- ////
//////////////////////////////////////
// Note: The BetaVec boost vector is calculated when someone calls
// 		the SetBitBoostBack() function

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


///////////////////////////////////////////////////////////////////////
ThreeBodyDecay::~ThreeBodyDecay()
{
	delete P;
	delete p[0];
	delete p[1];
	delete p[2];
}


///////////////////////////////////////////////////////////////////////
void ThreeBodyDecay::DisplayTag()
{
	printf("\n### %s ###\n", tag.c_str());
}

///////////////////////////////////////////////////////////////////////
void ThreeBodyDecay::SetTag(const char tag_[])
{
	tag = tag_;
}


///////////////////////////////////////////////////////////////////////
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


///////////////////////////////////////////////////////////////////////
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


///////////////////////////////////////////////////////////////////////
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


///////////////////////////////////////////////////////////////////////
void ThreeBodyDecay::SetPhaseSpace(double x1_, double x2_, double x3_, double x4_, double x5_)
{

	x1 = x1_;
	x2 = x2_;
	x3 = x3_;
	x4 = x5_;
	x5 = x5_;

		  s23  = s23_min + s23_length*x1;
   s23_sqrt  = sqrt(s23);
 //costheta1 = 1-2*x2;
	costheta1 = 2*x2-1;
	sintheta1 = sqrt(1-costheta1*costheta1);
	     phi1 = 2*M_PI*x3;
	  cosphi1 = cos(phi1);
	  sinphi1 = sin(phi1);

 //costheta23 = 1-2*x4;
	costheta23 = 2*x4-1;
	sintheta23 = sqrt(1-costheta23*costheta23);
	     phi23 = 2*M_PI*x5;
	  cosphi23 = cos(phi23);
	  sinphi23 = sin(phi23);

	// In the rest frame of the mother particle
	betabar = TwoBodyFunc::beta ( M_sqr, m_sqr[0], s23);
	lambda  = TwoBodyFunc::lambda ( M_sqr, m_sqr[0], s23);
	lambda_sqrt  = sqrt(lambda);
	pmag    = M*betabar/2;
	E1   = TwoBodyFunc::E ( M*M, m_sqr[0], s23);
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
	BetaVec = P->BoostVector();

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


////////////////////////////////////////////
double ThreeBodyDecay::GetPSConst()
{
	return PSConst;
}


////////////////////////////////////////////
double ThreeBodyDecay::GetPhaseSpaceWeight(double x1_, double x2_, double x3_, double x4_, double x5_)
{
	ThreeBodyDecay::SetPhaseSpace(x1_,x2_,x3_,x4_,x5_);
	return PSWeight;
}


////////////////////////////////////////////
void ThreeBodyDecay::DisplayAll()
{
	printf("\nDisplaying all info...\n");
	ThreeBodyDecay::DisplayMasses();
	printf("\n");
	ThreeBodyDecay::DisplayMomenta();
	printf("\n");
	ThreeBodyDecay::DisplayKinematics();
}


////////////////////////////////////////////
void ThreeBodyDecay::DisplayConfiguration()
{
	ThreeBodyDecay::DisplayTag();
	printf("x1 = %5.2f <->          s23 = %5.2f \n", x1, s23);
	printf("x2 = %5.2f <->  cos(theta1) = %5.2f \n", x2, costheta1);
	printf("x3 = %5.2f <->         phi1 = %5.2f \n", x3, phi1);
	printf("x4 = %5.2f <-> cos(theta23) = %5.2f \n", x4, costheta23);
	printf("x5 = %5.2f <->        phi23 = %5.2f \n", x5, phi23);
}


////////////////////////////////////////////
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
	std::cout << "PSWeight " << PSWeight << std::endl;
	sum = (*p[1]) + (*p[2]);
	std::cout << "(*p[1]) + (*p[2]).M(): " << sum.M() << std::endl;
}


////////////////////////////////////////////
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


/////////////////////////////////////
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


/////////////////////////////////////
//// --- DecayChain123 class --- ////
/////////////////////////////////////


DecayChain126::DecayChain126(double sqrt_s_, double mA_,  double mB_,
					  				  double mA1_,    double mA2_, double mA3_,
					  				  double mB1_,    double mB2_, double mB3_)
{

	sqrt_s = sqrt_s_;
	mA     = mA_;
	mB     = mB_;
	mB1    = mB1_;
	mB2    = mB2_;
	mB3    = mB3_;
	mA1    = mA1_;
	mA2    = mA2_;
	mA3    = mA3_;

	s       = sqrt_s * sqrt_s;
	mA_sqr  = mA  * mA ;
	mB_sqr  = mB  * mB ;
	mB1_sqr = mB1 * mB1;
	mB2_sqr = mB2 * mB2;
	mB3_sqr = mB3 * mB3;
	mA1_sqr = mA1 * mA1;
	mA2_sqr = mA2 * mA2;
	mA3_sqr = mA3 * mA3;

	// Creating sub-decay classes
	DecayAB   = new   TwoBodyDecay( sqrt_s, mA,  mB      );
	DecayA123 = new ThreeBodyDecay(     mA, mA1, mA2, mA3);
	DecayB123 = new ThreeBodyDecay(     mB, mB1, mB2, mB3);

	// Set the name of the tags
	DecayAB->SetTag("DecayAB");
	DecayA123->SetTag("DecayA123");
	DecayB123->SetTag("DecayB123");

	// Boost back to the original frame
	DecayAB->SetBitBoostBack(true);
	DecayA123->SetBitBoostBack(true);
	DecayB123->SetBitBoostBack(true);

	// Linking pointers
	//
	P = DecayAB->P;
	
	// pA and pB are determined from DecayAB
	pA = DecayAB->p1;
	pB = DecayAB->p2;

	// DecayA123 and DecayB123 total momentum is passed from pA and pB
	DecayA123->P = pA;
	DecayB123->P = pB;

	pA1 = DecayA123->p[0];
	pA2 = DecayA123->p[1];
	pA3 = DecayA123->p[2];

	pB1 = DecayB123->p[0];
	pB2 = DecayB123->p[1];
	pB3 = DecayB123->p[2];

	BitBoostBack = false; // - default value

	// Calc DecayChain126 PSConst
	PSConst = DecayAB->GetPhaseSpaceWeight() * DecayA123->GetPSConst() * DecayB123->GetPSConst();

};


///////////////////////////////////
DecayChain126::~DecayChain126()
{
	delete P;
	delete pA;
	delete pB;
	delete pA1;
	delete pA2;
	delete pA3;
	delete pB1;
	delete pB2;
	delete pB3;
};


///////////////////////////////////////////////////////////////////////
void DecayChain126::SetTag(const char tag_[])
{
	tag = tag_;
}

///////////////////////////////////////////////////////////////////////
void DecayChain126::DisplayTag()
{
	printf("\n### %s ###\n", tag.c_str());
}

///////////////////////////////////
void DecayChain126::SetTotalThreeMomentum_pxpypz( double px, double py, double pz )
{
	// Total energy of the (boosted system)/(decaying mother particle) viewed from the LAB frame
	E = sqrt( s + px*px + py*py + pz*pz );
	P->SetPxPyPzE(px,py,pz,E);

	# if DEBUG
	printf("\n\nSTART DEBUG INFO\n");
	printf("BetaVec components:\n");
	printf("bx = %12.6e\n", BetaVec.x());
	printf("by = %12.6e\n", BetaVec.y());
	printf("bz = %12.6e\n", BetaVec.z());
	printf("END DEBUG INFO\n\n");
	# endif
}


///////////////////////////////////
void DecayChain126::SetTotalThreeMomentum_pthetaphi(double p_, double theta_, double phi_)
{
	E = sqrt( s + p_*p_ );
	P->SetPxPyPzE(p_,0.0,0.0,E);
	P->SetTheta(theta_);
	P->SetPhi(phi_);

	# if DEBUG
	printf("\n\nSTART DEBUG INFO\n");
	printf("DecayChain126::SetMotherMPThetaPhi(double mom, double theta, double phi)\n");
	printf("BetaVec:\n");
	printf("bx = %12.6e\n", BetaVec.x());
	printf("by = %12.6e\n", BetaVec.y());
	printf("bz = %12.6e\n", BetaVec.z());
	printf("END DEBUG INFO\n\n");
	# endif
}


///////////////////////////////////
void DecayChain126::SetBitBoostBack(bool flag)
{
	BitBoostBack = flag;

	if (BitBoostBack == true )
	{
		printf("\nBoost back: yes\n");
		DecayAB->SetBitBoostBack(flag);
	}
	else if ( BitBoostBack == false )
	{
		printf("\nBoost back: no\n");
		DecayAB->SetBitBoostBack(flag);
	}
//	BetaVec = P->BoostVector();
}


///////////////////////////////////
void DecayChain126::SetPhaseSpace(  double x1_,  double x2_,  double x3_,
					 						   double x4_,  double x5_,  double x6_,	
					 						   double x7_,  double x8_,  double x9_,	
					 						  double x10_, double x11_, double x12_)
{

	// Setting internal variables
	x1  = x1_;
	x2  = x2_;
	x3  = x3_;
	x4  = x4_;
	x5  = x5_;
	x6  = x6_;
	x7  = x7_;
	x8  = x8_;
	x9  = x9_;
	x10 = x10_;
	x11 = x11_;
	x12 = x12_;
	
//	DecayA123->SetPhaseSpace(x3, x4,  x5,  x6,  x7);
//	DecayB123->SetPhaseSpace(x8, x9, x10, x11, x12);

	DecayAB->SetPhaseSpace(x1, x2);
	double PSWeight_DecayA123 = DecayA123->GetPhaseSpaceWeight(x3, x4, x5, x6, x7);     // - sets sub phase space and return the weight
	double PSWeight_DecayB123 = DecayB123->GetPhaseSpaceWeight(x8, x9, x10, x11, x12);  // - sets sub phase space and return the weight

	PSWeight = PSWeight_DecayA123 * PSWeight_DecayB123; // - you need to multiply this by PSConst!!


// You don't need this as the BitBoostBack flags are already set 
// in the initialization and take care of boosting all the particles
//	if( BitBoostBack )
//	{
//
//		# if DEBUG
//		printf("BetaVec:\n");
//		printf("bx = %12.6e\n", BetaVec.x());
//		printf("by = %12.6e\n", BetaVec.y());
//		printf("bz = %12.6e\n", BetaVec.z());
//		printf("END DEBUG INFO\n\n");
//		# endif
//
//		pA->Boost(BetaVec);
//		pB->Boost(BetaVec);
//		pA1->Boost(BetaVec);
//		pA2->Boost(BetaVec);
//		pA3->Boost(BetaVec);
//		pB1->Boost(BetaVec);
//		pB2->Boost(BetaVec);
//		pB3->Boost(BetaVec);
//
//	}

}

///////////////////////////////////////////////////////////////////////////////
double DecayChain126::GetPhaseSpaceWeight(  double x1_,  double x2_,  double x3_,
					 						   double x4_,  double x5_,  double x6_,	
					 						   double x7_,  double x8_,  double x9_,	
					 						  double x10_, double x11_, double x12_)
{

	DecayChain126::SetPhaseSpace( x1, x2,
						 	   x3, x4,  x5,  x6,  x7,
								x8, x9, x10, x11, x12);


	return PSWeight;

}

////////////////////////////////////
double DecayChain126::GetPSConst()
{
	return PSConst;
}

///////////////////////////////////
void DecayChain126::DisplayAllInfo()
{
	DecayChain126::DisplayMasses();
	DecayChain126::DisplayConfiguration();
	DecayAB->DisplayConfiguration();
	DecayA123->DisplayConfiguration();
	DecayB123->DisplayConfiguration();
	DecayChain126::DisplayMomenta();
}

///////////////////////////////////
void DecayChain126::DisplayMomenta()
{

	TLorentzVector sum;
	TLorentzVector sumA;
	TLorentzVector sumB;

	sum = (*pA) + (*pB);
	sumA = (*pA1) + (*pA2) + (*pA3); 
	sumB = (*pB1) + (*pB2) + (*pB3); 

	printf("\n-- Momenta Information --\n") ;
	std::cout << "P:" << std::endl;
	displayTLorentzVector(P);

	std::cout << "pA:" << std::endl;
	displayTLorentzVector(pA);
	std::cout << "pB:" << std::endl;
	displayTLorentzVector(pB);

	std::cout << "pA1:" << std::endl;
	displayTLorentzVector(pA1);
	std::cout << "pA2:" << std::endl;
	displayTLorentzVector(pA2);
	std::cout << "pA3:" << std::endl;
	displayTLorentzVector(pA3);

	std::cout << "pB1:" << std::endl;
	displayTLorentzVector(pB1);
	std::cout << "pB2:" << std::endl;
	displayTLorentzVector(pB2);
	std::cout << "pB3:" << std::endl;
	displayTLorentzVector(pB3);

	printf("\nInternal test with sum of specific momenta\n") ;
	std::cout << "pA+pB:" << std::endl;
	displayTLorentzVector(&sum);

	std::cout << "pA1+pA2+pA3:" << std::endl;
	displayTLorentzVector(&sumA);

	std::cout << "pB1+pB2+pB3:" << std::endl;
	displayTLorentzVector(&sumB);
	printf("\n");

}

///////////////////////////////////////////
void DecayChain126::DisplayConfiguration()
{
	DecayChain126::DisplayTag();
	printf("x1  = %.2f = DecayAB   x1\n", x1);
	printf("x2  = %.2f = DecayAB   x2\n", x2);
	printf("x3  = %.2f = DecayA123 x1\n", x3);
	printf("x4  = %.2f = DecayA123 x2\n", x4);
	printf("x5  = %.2f = DecayA123 x3\n", x5);
	printf("x6  = %.2f = DecayA123 x4\n", x6);
	printf("x7  = %.2f = DecayA123 x5\n", x7);
	printf("x8  = %.2f = DecayB123 x1\n", x8);
	printf("x9  = %.2f = DecayB123 x2\n", x9);
	printf("x10 = %.2f = DecayB123 x3\n", x10);
	printf("x11 = %.2f = DecayB123 x4\n", x11);
	printf("x12 = %.2f = DecayB123 x5\n", x12);
}


///////////////////////////////////////////
void DecayChain126::DisplayMasses()
{
	DecayChain126::DisplayTag();
	printf("sqrt_s  = %.2f\n", sqrt_s);
	printf("mA      = %.2f\n", mA);
	printf("mB      = %.2f\n", mB);
	printf("mA1     = %.2f\n", mA1);
	printf("mA2     = %.2f\n", mA2);
	printf("mA3     = %.2f\n", mA3);
	printf("mB1     = %.2f\n", mB1);
	printf("mB2     = %.2f\n", mB2);
	printf("mB3     = %.2f\n", mB3);
}
