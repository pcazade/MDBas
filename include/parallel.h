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
  
  void close_para();

  void mpi_error(int err, char file[],int line);
#ifdef	__cplusplus
}
#endif

#endif
