#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>

#define PI (3.14159265358979)

#define X2(x) ((x)*(x))

extern void dsyev_(char *jobz, char *uplo, int *n, double *A, int *lda,
           double *W, double *WORK, int *lwork, int *info);


int main(int argc, char* argv[])
{

  double a,b,c,alpha,beta,gamma,aldeg,bedeg,gadeg;
  double egvec[9],hth[9],crystal[6];
  
  char jobz = 'V',uplo='L';
  int n = 3 ;
  int lda = n;
  double *W = malloc(n*sizeof(double));
  double *WORK = NULL;
  int lwork = -1;
  int info=-42;

  double siz = 0.0;
  
  if(argc>=7)
  {
    a=atof(argv[1]);
    b=atof(argv[2]);
    c=atof(argv[3]);
    
    aldeg=atof(argv[4]);
    bedeg=atof(argv[5]);
    gadeg=atof(argv[6]);
    
    alpha=aldeg*PI/180.;
    beta=bedeg*PI/180.;
    gamma=gadeg*PI/180.;
    
  }
  else
  {
    printf("%s A B C Alpha Beta Gamma\n",argv[0]);
    exit(0);
  }
  
  hth[0]=X2(a);
  hth[4]=X2(b);
  hth[8]=X2(c);
  
  if(fabs(aldeg-90.)<DBL_EPSILON)
  {
    hth[5]=0.;
    hth[7]=0.;
  }
  else
  {
    hth[5]=b*c*cos(alpha);
    hth[7]=hth[5];
  }
  
  if(fabs(bedeg-90.)<DBL_EPSILON)
  {
    hth[2]=0.;
    hth[6]=0.;
  }
  else
  {
    hth[2]=a*c*cos(beta);
    hth[6]=hth[2];
  }
  
  if(fabs(gadeg-90.)<DBL_EPSILON)
  {
    hth[1]=0.;
    hth[3]=0.;
  }
  else
  {
    hth[1]=a*b*cos(gamma);
    hth[3]=hth[1];
  }
  
  memcpy(egvec,hth,9*sizeof(double));
  
  
  /**
   * Diagonalization of the HTH matrix with lapack routine:
   * 
   * http://www.netlib.org/lapack/explore-html/d2/d8a/group__double_s_yeigen.html
   *
  subroutine DSYEV    (    CHARACTER     JOBZ,
                CHARACTER     UPLO,
             INTEGER     N,
             DOUBLE PRECISION, dimension( lda, * )     A,
                INTEGER     LDA,
             DOUBLE PRECISION, dimension( * )     W,
                DOUBLE PRECISION, dimension( * )     WORK,
                INTEGER     LWORK,
             INTEGER     INFO
  )
  **/
  
  dsyev_(&jobz, &uplo, &n, egvec, &lda, W, &siz, &lwork, &info);
  
  if(info==0)
    WORK = (double*)malloc((size_t)siz*sizeof(double));
  else
  {
    printf("Error with DSYEV when initialising work array : lapack error code %d .\n",info);
    exit(42);
  }

  if (WORK == NULL)
  {
    printf("Error when allocating a array of double of %d elements; FILE %s LINE %d",(int)siz,__FILE__,__LINE__);
    exit(22);
  }

  lwork = (int)siz;
  
  dsyev_(&jobz, &uplo, &n, egvec, &lda, W, WORK, &lwork, &info);
  
  if(info != 0)
  {
    printf("Error with DSYEV after return : lapack error code %d .\n",info);
    exit(2012);
  }
  
  W[0]=sqrt(W[0]);
  W[1]=sqrt(W[1]);
  W[2]=sqrt(W[2]);
  
//   printf("%lf\t%lf\t%lf\n\n",W[0],W[1],W[2]);

//   memcpy(hth,egvec,9*sizeof(double));
  
  crystal[0]=W[0]*X2(egvec[0])+W[1]*X2(egvec[3])+W[2]*X2(egvec[6]);
  crystal[2]=W[0]*X2(egvec[1])+W[1]*X2(egvec[4])+W[2]*X2(egvec[7]);
  crystal[5]=W[0]*X2(egvec[2])+W[1]*X2(egvec[5])+W[2]*X2(egvec[8]);
  
  crystal[1]=W[0]*egvec[0]*egvec[1]+W[1]*egvec[3]*egvec[4]+W[2]*egvec[6]*egvec[7];
  crystal[3]=W[0]*egvec[0]*egvec[2]+W[1]*egvec[3]*egvec[5]+W[2]*egvec[6]*egvec[8];
  crystal[4]=W[0]*egvec[1]*egvec[2]+W[1]*egvec[4]*egvec[5]+W[2]*egvec[7]*egvec[8];
  
  printf("%lf\t%lf\t%lf\n",crystal[0],crystal[1],crystal[3]);
  printf("%lf\t%lf\t%lf\n",crystal[1],crystal[2],crystal[4]);
  printf("%lf\t%lf\t%lf\n",crystal[3],crystal[4],crystal[5]);

  free(WORK);
  free(W);
}
