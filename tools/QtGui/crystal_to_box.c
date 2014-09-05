#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>

#include "crystal_to_box.h"

#define PI (3.14159265358979)
#define X2(x) ((x)*(x))

extern void dsyev_(char *jobz, char *uplo, int *n, real *A, int *lda, real *W, real *WORK, int *lwork, int *info);

int lattice_to_cryst(real lattice[6], real matrix[9])
{

    real a,b,c,alpha,beta,gamma,aldeg,bedeg,gadeg;
    real egvec[9],hth[9],crystal[6];

    char jobz = 'V',uplo='L';
    int n = 3 ;
    int lda = n;
    real *W = malloc(n*sizeof(real));
    real *WORK = NULL;
    int lwork = -1;
    int info=-42;

    real siz = 0.0;

    a=lattice[0];
    b=lattice[1];
    c=lattice[2];

    aldeg=lattice[3];
    bedeg=lattice[4];
    gadeg=lattice[5];

    alpha=aldeg*PI/180.;
    beta=bedeg*PI/180.;
    gamma=gadeg*PI/180.;

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

    memcpy(egvec,hth,9*sizeof(real));


    /*
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
    */

    dsyev_(&jobz, &uplo, &n, egvec, &lda, W, &siz, &lwork, &info);

    if(info==0)
        WORK = (real*)malloc((size_t)siz*sizeof(real));
    else
        return info;

    lwork = (int)siz;

    dsyev_(&jobz, &uplo, &n, egvec, &lda, W, WORK, &lwork, &info);

    if(info != 0)
        return info;

    W[0]=sqrt(W[0]);
    W[1]=sqrt(W[1]);
    W[2]=sqrt(W[2]);

    crystal[0]=W[0]*X2(egvec[0])+W[1]*X2(egvec[3])+W[2]*X2(egvec[6]);
    crystal[2]=W[0]*X2(egvec[1])+W[1]*X2(egvec[4])+W[2]*X2(egvec[7]);
    crystal[5]=W[0]*X2(egvec[2])+W[1]*X2(egvec[5])+W[2]*X2(egvec[8]);

    crystal[1]=W[0]*egvec[0]*egvec[1]+W[1]*egvec[3]*egvec[4]+W[2]*egvec[6]*egvec[7];
    crystal[3]=W[0]*egvec[0]*egvec[2]+W[1]*egvec[3]*egvec[5]+W[2]*egvec[6]*egvec[8];
    crystal[4]=W[0]*egvec[1]*egvec[2]+W[1]*egvec[4]*egvec[5]+W[2]*egvec[7]*egvec[8];

    matrix[0] = crystal[0];
    matrix[1] = crystal[1];
    matrix[2] = crystal[3];
    matrix[3] = crystal[1];
    matrix[4] = crystal[2];
    matrix[5] = crystal[4];
    matrix[6] = crystal[3];
    matrix[7] = crystal[4];
    matrix[8] = crystal[5];

    free(WORK);
    free(W);

    return info;
}
