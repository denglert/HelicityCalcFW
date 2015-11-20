// Define Fortran functions //

extern "C"
{
	void ranmar_(float *rvec, int *lenv);
	void rambo_(int *nParticles, double *E_cm, double *masses, double *p, double *wt, int *lw);
	void vegas_(double *region, int *ndim, double (*fxn)(double *x), int *init, int *ncall,
				 	int *itmx, int *nprn, double *tgral, double *sd, double *chi2a, double *acc,
					double **xi, int *it, int *ndo, double *si, double *swgt, double *schi );

	double rh_tautau_(double *p1, double *p2, int *i1, int *i2);
	double rh_6f_(double *p3, double *p4, double *p5, double *p6, double *p7, double *p8);

   extern struct
	{
		double rmtau;		
   } masses_;	

   extern struct
	{
		double wcl;
		double gh_tautau;		
   } couplings_;	

   extern struct
	{
		double taum;
		double taup;
	} amplitudes_;
}
