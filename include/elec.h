#ifndef ELECH_INCLUDED
#define ELECH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

double coulomb_none(const PARAM *param,double *delec,const double qel,
		      const double r2,const double rt);

void coulomb_full(ENERGY *ener,PARAM *param,PBC *box,double *x,double *y,
		  double *z,double *fx, double *fy, double *fz,double *q,
		  int **exclList,int *exclPair);

double coulomb_shift1(const PARAM *param,double *delec,const double qel,
		      const double r2,const double rt);

double coulomb_shift2(const PARAM *param,double *delec,const double qel,
		      const double r2,const double rt);

double coulomb_switch(const PARAM *param,double *delec,const double qel,
		      const double r2,const double rt);

double coulomb14_none(const PARAM *param,double *delec,const double qel,
		      const double r2,const double rt);

double coulomb14_full(const PARAM *param,double *delec,const double qel,
		      const double r2,const double rt);

double coulomb14_shift1(const PARAM *param,double *delec,const double qel,
			const double r2,const double rt);

double coulomb14_shift2(const PARAM *param,double *delec,const double qel,
			const double r2,const double rt);

double coulomb14_switch(const PARAM *param,double *delec,const double qel,
			const double r2,const double rt);

/** Pointer to the output file. **/
extern FILE *outFile;

#ifdef	__cplusplus
}
#endif

#endif