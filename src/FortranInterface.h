// Define Fortran functions //
extern "C"
{
	void ranmar_(float *rvec, int *lenv);
	void rambo_(int *nParticles, double *E_cm, double *masses, double *p, double *wt, int *lw);
	double rh_tautau_(double *p1, double *p2);
   extern struct
	{
		double rmtau;		
   } masses_;	
}
