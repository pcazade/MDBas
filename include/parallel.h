#ifndef PARALLELH_INCLUDED
#define PARALLELH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

    void init_para(int *argc, char ***argv);

  int my_proc();

  int num_proc();
  
  void barrier_para();
  
  void sum_double_para(double *buf1,double *buf2,int size);
  
  void sum_int_para(int *buf1,int *buf2,int size);
  
  void update_double_para(PARAM *param,PARALLEL *parallel,double *buf1,double *buf2);
  
  void test_para(int *buf1);
  
  void parallel_allocate_buffers(PARAM *param,PARALLEL *parallel,double **dBuffer,
				 int **iBuffer);
  
  void parallel_reallocate_buffers(PARALLEL *parallel,EWALD *ewald,double **dBuffer);
  
  void bcast_int_para(int *buf1,int size,int iNode);
  
  void bcast_double_para(double *buf1,int size,int iNode);
  
  void bcast_param_para(PARAM *param,PARALLEL *parallel,double *dBuffer,int *iBuffer);
  
  void bcast_ctrl_para(CTRL *ctrl,PARALLEL *parallel,int *iBuffer);
  
  void bcast_bath_para(BATH *bath,PARALLEL *parallel,double *dBuffer);
  
  void bcast_pbc_para(PBC *box,PARALLEL *parallel,double *dBuffer);
  
  void bcast_neigh_para(NEIGH *neigh,PARALLEL *parallel,int *iBuffer);
  
  void bcast_ewald_para(EWALD *ewald,PARALLEL *parallel,int *iBuffer,double *dBuffer);
  
  void bcast_const_para(CONSTRAINT *constList,PARALLEL *parallel,int *iBuffer,double *dBuffer,int size);
  
  void bcast_bond_para(BOND *bond,PARALLEL *parallel,int *iBuffer,double *dBuffer,int size);
  
  void bcast_angle_para(ANGLE *angle,PARALLEL *parallel,int *iBuffer,double *dBuffer,int size);
  
  void bcast_dihe_para(DIHE *dihe,PARALLEL *parallel,int *iBuffer,double *dBuffer,int size);
  
  void setup_para(CTRL *ctrl,PARAM *param,PARALLEL *parallel,ENERGY *ener,
		  BATH *bath,NEIGH *neigh,EWALD *ewald,PBC *box,ATOM **atom,CONSTRAINT **constList,
		  BOND **bond,ANGLE **angle,DIHE **dihe,DIHE **impr,BOND **ub,double **x,
		  double **y, double **z,double **vx,double **vy,double **vz,double **fx,
		  double **fy, double **fz,double **mass,double **rmass,double **q,
		  double **eps,double **sig,double **eps14,double **sig14,int **frozen,
		  int **nAtConst,double **dBuffer,int **iBuffer);
  
  void close_para();

  void mpi_error(int err, char file[],int line);

#ifdef	__cplusplus
}
#endif

#endif
