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

	double momentum  = 0.0;
	double theta 	  = 0.53;
	double phi  	  = 1.34;
	double 	     M  = 10.0;
	double 	Energy  = sqrt(momentum*momentum+M*M);
	double 	    m1  = 0.0000;
	double 	    m2  = 0.0000;
	double 	    m3  = 0.0000;

//	double x1 = 0.5;
//	double x2 = 0.2;
//	double x3 = 0.4;
//	double x4 = 0.3342;
//	double x5 = 0.75;
//

	double x1 = 0.50;
	double x2 = 0.00;
	double x3 = 0.00;
	double x4 = 0.00;
	double x5 = 0.00;

	{
	std::cout << "test_PhaseSpaceTools " << std::endl;

	TLorentzVector p1;
	TLorentzVector p2;
	TLorentzVector p3;

	TVector3 v1;

	double costheta1 = 2*x2-1;
	double sintheta1 = sqrt(1-costheta1*costheta1);
	double 	   phi1 = 2*M_PI*x3;
	double sqrt_s23  = x1*(M-m1-m2-m3) + m1+m2;

	double costheta23 = 2*x4-1;
	double sintheta23 = sqrt(1-costheta23*costheta23);
	double 	   phi23 = 2*M_PI*x5;

	double      s = M*M;
	double m1_sqr = m1*m1;
	double m2_sqr = m2*m2;
	double m3_sqr = m3*m3;
	double 	 s23 = sqrt_s23*sqrt_s23;

	double    pmag  = TwoBodyFunc::p    (s, m1_sqr, s23);
	double beta1    = TwoBodyFunc::beta (s, m1_sqr, s23);
	double 	 E1    = TwoBodyFunc::E    (M, m1_sqr, s23);
	double 	 E23   = TwoBodyFunc::E    (M, s23, m1_sqr);

	p1.SetPxPyPzE(pmag*sintheta1*cos(phi1),pmag*sintheta1*sin(phi1),pmag*costheta1,E1);

	TVector3 v;
	v = -p1.Vect();
	v.SetMag(v.Mag()/E23);

	double    p2mag  = TwoBodyFunc::p    (s23, m2_sqr, m3_sqr);
	double    	 E2  = TwoBodyFunc::E    (sqrt_s23, m2_sqr, m3_sqr);
	double    	 E3  = TwoBodyFunc::E    (sqrt_s23, m3_sqr, m2_sqr);

	p2.SetPxPyPzE(p2mag*sintheta23*cos(phi23),p2mag*sintheta23*sin(phi23),p2mag*costheta23,E2);
	p3.SetPxPyPzE(-p2mag*sintheta23*cos(phi23),-p2mag*sintheta23*sin(phi23),-p2mag*costheta23,E3);

	std::cout << "p1:" << std::endl;
	displayTLorentzVector(&p1);

	std::cout << "p2 and p3 in the 23 rest frame:" << std::endl;

	std::cout << "p2:" << std::endl;
	displayTLorentzVector(&p2);

	std::cout << "p3:" << std::endl;
	displayTLorentzVector(&p3);

	p2.Boost(v);
	p3.Boost(v);


	std::cout << "All in the lab frame:" << std::endl;

	std::cout << "p1:" << std::endl;
	displayTLorentzVector(&p1);

	std::cout << "p2:" << std::endl;
	displayTLorentzVector(&p2);

	std::cout << "p3:" << std::endl;
	displayTLorentzVector(&p3);

	TLorentzVector sum;
	sum = p1 + p2 + p3;

	std::cout << "sum in the lab frame:" << std::endl;
	displayTLorentzVector(&sum);

	}

	{

	printf("######################################\n");
	printf("### --- TwoBudyFunc utils TEST --- ###\n");
	printf("######################################\n");

	double M_test      = 5.0;
	double m1_sqr_test = 3.0;
	double m2_sqr_test = 1.0;

	double E = TwoBodyFunc::E (M_test, m1_sqr_test, m2_sqr_test);

	printf("M_test: %.2f\n", M_test);
	printf("m1_sqr_test: %.2f\n", m1_sqr_test);
	printf("m2_sqr_test: %.2f\n", m2_sqr_test);
	printf("TwoBodyFunc::E (M_test, m1_sqr_test, m2_sqr_test): %.2f\n", E);

	}

	{

	printf("\n");
	printf("########################################\n");
	printf("### --- ThreeBodyDecayClass TEST --- ###\n");
	printf("########################################\n");
	printf("\n");
	// ThreeBodyDecay class
   ThreeBodyDecay tau(M, m1, m2, m3);
	tau.SetMotherMPThetaPhi(M,momentum,theta,phi);
	tau.SetBitBoostBack(true);


	tau.SetPhaseSpace(x1, x2, x3, x4, x5);
	TLorentzVector *P = tau.P;
	TLorentzVector *p1 = tau.p[0];
	TLorentzVector *p2 = tau.p[1];
	TLorentzVector *p3 = tau.p[2];
  	double amp = P->Dot( (*p3) ) * p1->Dot( (*p2) );
	tau.DisplayAll();

//	printf("\ns23_min: %.2f", );

	double weight = tau.GetPhaseSpaceWeight(x1, x2, x3, x4, x5);
	std::cout << "PSWeight: " << weight << std::endl;
	std::cout << "amp: " << amp << std::endl;

	std::cout << "P:" << std::endl;
	displayTLorentzVector(P);

	std::cout << "p1:" << std::endl;
	displayTLorentzVector(p1);
	std::cout << "p2:" << std::endl;
	displayTLorentzVector(p2);
	std::cout << "p3:" << std::endl;
	displayTLorentzVector(p3);

	TLorentzVector sum = (*p1) + (*p2) + (*p3);
	std::cout << "p1+p2+p3:" << std::endl;
	displayTLorentzVector(&sum);
	}

}
