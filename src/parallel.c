#ifdef MDBAS_PARA_MPI

#include <mpi.h>

#define BUFSIZ 8192

int init_para(int *argc, char ***argv,PARAM *param)
{
  
  MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
  
  err=MPI_init(argc,argv);
  
  if(err!=MPI_SUCCES)
    call mpi_error(err,__FILE__,__LINE__);
  
  param->nProc=num_proc();
  param->nAtProc=(param->nAtoms+param->nProc-1)/param->nProc;
  
  
}

int myproc()
{
  
  int idProc,err;
  
  err=MPI_Comm_rank(MPI_COMM_WORLD,&idProc);
  
  if(err!=MPI_SUCCES)
    call mpi_error(err,__FILE__,__LINE__);
  
  return (idProc);
  
}


int num_proc()
{
  int numProc,err;
  
  err=MPI_Comm_size(MPI_COMM_WORLD, &numProc));
  
  if(err!=MPI_SUCCES)
    call mpi_error(err,__FILE__,__LINE__);
  
  return(numProc);
  
}

void mpi_error(int err, char file[],int line)
{
  char errString[BUFSIZ];
  int idProc,lenErrString, errClass;
  
  idProc=myproc();
  
  MPI_Error_class(err, &errClass);
  
  MPI_Error_string(errClass, errString, &lenErrString);
  fprintf(outFile, "%3d: %s\n", idProc, errString);
  
  MPI_Error_string(err, errString, &lenErrString);
  fprintf(outFile, "%3d: %s\n", idProc, errString);
  
  MPI_Abort(MPI_COMM_WORLD, err);
  
  my_error(MPI_ERROR,file,line,0);
    
}

#endif //MDBAS_PARA_MPI