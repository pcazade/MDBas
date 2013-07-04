#ifndef PARALLELH_INCLUDED
#define PARALLELH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

  void init_para(int *argc, char ***argv,PARAM *param);

  int my_proc();

  int num_proc();
  
  void sum_double_para(double *buf1,double *buf2,int size);
  
  void sum_int_para(int *buf1,int *buf2,int size);
  
  void update_double_para(PARAM *param,PARALLEL *parallel,double *buf1,double *buf2);
  
  void test_para(int *buf1);
  
  void setup_para(CTRL *ctrl,PARAM *param,PARALLEL *parallel,ENERGY *ener,
		  BATH *bath,NEIGH *neigh,EWALD *ewald,PBC *box,ATOM **atom,CONSTRAINT **constList,
		  BOND **bond,ANGLE **angle,DIHE **dihe,DIHE **impr,BOND **ub,double **x,
		  double **y, double **z,double **vx,double **vy,double **vz,double **fx,
		  double **fy, double **fz,double **mass,double **rmass,double **q,
		  double **eps,double **sig,double **eps14,double **sig14,int **frozen,
		  int **nAtConst,int **neighList,int **neighPair,int **neighList14,
		  int ***exclList,int **exclPair);
  
  void close_para();

  void mpi_error(int err, char file[],int line);
#ifdef	__cplusplus
}
#endif

#endif
