#include "HelicityTools.h"
#include "FortranInterface.h"

double HelicityTools::rh_tautau_polarized(double *p1, double *p2, int &pol1, int &pol2)
{
	double res = rh_tautau_(p1,p2,&pol1,&pol2);
};

double HelicityTools::rh_tautau_unpolarized(double *p1, double *p2)
{
	double res = 0;
	for(int pol1=1; pol1 < 3; pol1++)	
	for(int pol2=1; pol2 < 3; pol2++)	
	{ res = res + rh_tautau_(p1,p2,&pol1,&pol2); }
	return res;
};
