// Define Fortran functions //

extern "C"
{
	void ranmar_(float *rvec, int *lenv);
	void rambo_(int *nParticles, double *E_cm, double *masses, double *p, double *wt, int *lw);
	void vegas_(double *region, int *ndim, double (*fxn)(double *x), int *init, int *ncall,
				 	int *itmx, int *nprn, double *tgral, double *sd, double *chi2a, double *acc,
					double **xi, int *it, int *ndo, double *si, double *swgt, double *schi );
	double rh_tautau_(double *p1, double *p2);

   extern struct
	{
		double rmtau;		
   } masses_;	
}
