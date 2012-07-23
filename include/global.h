#ifndef GLOBALH_INCLUDED
#define GLOBALH_INCLUDED

#define clight  (299792458.)
#define elemchg (1.602176565e-19)
#define angstr  (1.e-10)
#define calory  (4.184)
#define kcaltoiu (418.4)
#define NA	(6.02214129e+23)
#define kboltz  (1.3806488e-23)
#define rboltz  (8.3144621)
#define rboltzui  (0.83144621)
#define mu0     (1.e-7)
#define chgcharmm (332.0716)
#define chgnamd   (332.0636)
#define sq6rt2  (1.12246204830937)
#define PI      (3.14159265358979)

#define X2(x) ((x)*(x))
#define X3(x) (X2(x)*(x))
#define X4(x) (X2(x)*X2(x))
#define X6(x) (X3(x)*X3(x))
#define X12(x) (X6(x)*X6(x))
#define MAX(x,y) ((x)>=(y)?(x):(y))
#define MIN(x,y) ((x)<=(y)?(x):(y))

#define NOELEC 0
#define FULL   1
#define SHIFT1 2
#define SHIFT2 3
#define SWITCH 4

#define NOVDW   0
#define VFULL   1
#define VSWITCH 2

#define COSDIH  1
#define HARMDIH 2

typedef struct
{
  int natom;
  char **atomLabel,**segi,**resi;
  int *atomType,*ires,*resn,*inconst;
  double *x,*y,*z,*m;
  double *vx,*vy,*vz;
  double *fx,*fy,*fz;
  double *numFx,*numFy,*numFz;
}ATOM;

typedef struct
{
//   FILE *cntrFile,*topFile,*psfFile,*parmFile,*corFile;
  char **types;
  int *typesNum;
  int nTypes,nBondTypes,nAngTypes,nUbTypes,nDiheTypes,nImprTypes,nNonBonded;
  int **bondTypes,**angTypes,**ubTypes,*nDiheTypesParm,**diheTypes,**imprTypes;
  double **bondTypesParm,**angTypesParm,**ubTypesParm;
  double **diheTypesParm,**imprTypesParm,**nonBondedTypesParm;
}INPUTS;

typedef struct
{
  int nBond,nAngle,nDihedral,nImproper,nUb;
  double *q,*dq,*ddq,*qlpoly,**qpoly;
  int *verList,*verPair,npr,**ver14,npr14,*nParmDihe;
  double **parmVdw,scal14;
  double **parmBond,**parmUb,**parmAngle,**parmDihe;
  double **parmImpr;
}FORCEFIELD;

typedef struct
{
  int lqpoly,nb14,step,nsteps,degfree,firstener;
  int printo,printtr,integrator,ens,enstime;
  int keyrand,seed,keytraj,keyener,keyforf,keymd;
  int *excludeNum,*excludeAtom;
  int **iBond,**iUb,**iAngle,**iDihedral,**iImproper;
  int elecType,vdwType,periodicType,mdNature,numDeriv;
  int *bondType,*ubType,*angleType,*diheType,*imprType;
  double chargeConst,cutoff,cuton,delr;
  double temp,timeStep,periodicBox[3][3];
}SIMULPARAMS;

typedef struct
{
  double energyTot,energyPot,energyKin;
  double energyElec,energyVdw;
  double energyBond,energyAng,energyUb,energyDih,energyImpr;
  double **hessian;
}ENERGYFORCE;

typedef struct
{
  int i,j;
  double r0;
}CONSTRAINT;

typedef struct
{
  double x,y,z
}DELTA;

#endif
