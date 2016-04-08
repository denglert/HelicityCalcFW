#ifndef FORTRANINTERFACE_H
#define FORTRANINTERFACE_H

// Define Fortran functions //

struct cval
{
	double r;
	double i;
};

typedef struct cval cval;

typedef cval CMatrix_2_2[2][2] ;

extern "C"
{
	void ranmar_(float *rvec, int *lenv);
	void rambo_(int *nParticles, double *E_cm, double *masses, double *p, double *wt, int *lw);
	void vegas_(double *region, int *ndim, double (*fxn)(double *x), int *init, int *ncall,
				 	int *itmx, int *nprn, double *tgral, double *sd, double *chi2a, double *acc,
					double **xi, int *it, int *ndo, double *si, double *swgt, double *schi );

	double rh_tautau_(double *p1, double *p2, CMatrix_2_2 &cth);
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
		double rh_6f_tautau;
		double rh_6f_taum;
		double rh_6f_taup;
		double rh_6f_res_nwa;
		double rh_6f_res;
		double rh_6f_res_test;
	} amplitudes_;

   extern struct
	{
		CMatrix_2_2 ch_tautau;
		CMatrix_2_2 cdec_taum;
		CMatrix_2_2 cdec_taup;
		CMatrix_2_2 c7_568;
   } taumatrices_;	

   extern struct
	{
		CMatrix_2_2 c_amp_res;
		CMatrix_2_2 c_amp_htautau;
		CMatrix_2_2 c_amp_dec_taum;
		CMatrix_2_2 c_amp_dec_taup;
		CMatrix_2_2 c_amp_7_568;
   } tau_amplitudes_;	

}

#endif
