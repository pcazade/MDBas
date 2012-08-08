/***************************************************************
 * cmenerve is designed to compute the main energy terms for MM
 * simulations. The goal is to provide a toolkit of routines that
 * can be easily modified and combined.
 * 
 * P.-A. Cazade 02 May 2012
 * ************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "global.h"
#include "energy.h"
#include "io.h"
#include "rand.h"
#include "utils.h"
#include "integrate.h"
#include "numderiv.h"
#include "list.h"
#include "minim.h"


#if ( defined __linux__ || defined __FreeBSD__ )
#include <sys/resource.h>
#endif /* for some unixes */

#define _ENERFORMAT_ "%13d\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\n"
#define _ENERLABELS_ "%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\n"

int main(int argc, char* argv[])
{
  
//   void (*test)(int *par1,double *par2,...);
//   test=&(elec);

  FILE *fener=NULL,*traj=NULL,*frc=NULL,*numfrc=NULL;
  INPUTS inp;
  ATOM *atom=NULL;
  FORCEFIELD ff;
  ENERGY ener;
  SIMULPARAMS simulCond;
  CONSTRAINT *constList=NULL;
  DELTA *nForce=NULL;
  PBC box;
  
  char enerLabel[11][7]={ {"Step"}  , {"Etot"}  , {"Ekin"}  , {"Epot"} ,
			  {"Ecoul"} , {"Evdw"}  , {"Ebond"} , {"Eangle"} ,
			  {"Eub"}   , {"Edihe"} , {"Eimpr"} };
  int i;
  
  init_rand(time(NULL));
  
  read_SIMU(&simulCond,&ff,&box);
  
  printf("SIMU file read\n");
  
  if(simulCond.keyrand)
    init_rand(simulCond.seed);
  
  if(simulCond.mdNature==1)
    simulCond.chargeConst=chgcharmm*kcaltoiu;
  else if(simulCond.mdNature==2)
    simulCond.chargeConst=chgnamd*kcaltoiu;
  else
    simulCond.chargeConst=mu0*X2(clight)*X2(elemchg)*NA*0.1/(angstr);
  
  read_PSF(&inp,&atom,&ff,&simulCond,&constList);
  
  printf("PSF file read\n");
  
  read_TOP(&inp);
  
  printf("TOP file read\n");
  
  read_PAR(&inp);
  
  printf("PAR file read\n");
  
  read_CONF(atom,&simulCond);
  
  printf("CONF file read\n");
  
  setup(&inp,atom,&ff,&simulCond,constList);
  
  printf("Setup done\n");
  
  if(simulCond.keyforf)
    write_FORF(&inp,atom,&ff,&simulCond);
  
  free_temp_array(&inp);
  
  printf("Free temp array done\n");
  
  get_kinfromtemp(atom,&simulCond,&box);
  
  //init_vel(atom,&simulCond,constList,&box);
  
  init_box(&box);

  makelist(&simulCond,atom,&ff,constList,&box);
  
  if(simulCond.keyminim)
  {
    steepestDescent(atom,&ff,&ener,&simulCond,&box);
    makelist(&simulCond,atom,&ff,constList,&box);
  }
  
  init_vel(atom,&simulCond,constList,&box);
  
  if(simulCond.keyener)
  {
    fener=fopen("ener.dat","w");
    
    fprintf(fener,_ENERLABELS_,enerLabel[0],enerLabel[1],enerLabel[2],enerLabel[3],
			       enerLabel[4],enerLabel[5],enerLabel[6],enerLabel[7],
			       enerLabel[8],enerLabel[9],enerLabel[10]);
  }
  
  if(simulCond.keytraj)
  {
    traj=fopen("traj.xyz","w");
    fprintf(traj,"%d\n",simulCond.natom);
    fprintf(traj,"initial config (minimised)\n");
    for(i=0;i<simulCond.natom;i++)
      fprintf(traj,"%s\t%7.3lf\t%7.3lf\t%7.3lf\n",atom[i].label,atom[i].x,atom[i].y,atom[i].z);
  }
  
  if(simulCond.integrator==1 || !simulCond.keymd)
  {

//   Computes kinetic energy at time=0

    ener.kin=kinetic(atom,&simulCond);
  
//   Computes potential energies and forces at time=0
  
    energy(atom,&ff,&ener,&simulCond,&box);
    
    ener.tot=ener.kin+ener.pot;
    
    if(simulCond.keyener)
    { 
      fprintf(fener,_ENERFORMAT_,
	simulCond.step,ener.tot/kcaltoiu,ener.kin/kcaltoiu,ener.pot/kcaltoiu,
	ener.elec/kcaltoiu,ener.vdw/kcaltoiu,ener.bond/kcaltoiu,
	ener.ang/kcaltoiu,ener.ub/kcaltoiu,ener.dihe/kcaltoiu,
	ener.impr/kcaltoiu);
    }
    
//     if(simulCond.keytraj)
//     {
//       fprintf(traj,"%d\n",simulCond.natom);
//       fprintf(traj,"step %d\n",simulCond.step);
//       for(i=0;i<simulCond.natom;i++)
// 	  fprintf(traj,"%s\t%7.3lf\t%7.3lf\t%7.3lf\n",atom.atomLabel[i],atom[i].x,atom[i].y,atom[i].z);
//     }

//   Numerical derivatives to estimate forces for initial configuration.
      
    if(simulCond.numDeriv==1)
    {
      nForce=(DELTA*)malloc(simulCond.natom*sizeof(*nForce));
      
      frc=fopen("force.dat","w");
      numfrc=fopen("numforce.dat","w");
    
//     Write analytical forces into frc file.
    
      fprintf(frc,"#step=%d\n",simulCond.step);
      for(i=0;i<simulCond.natom;i++)
	fprintf(frc,"%.15lf\t%.15lf\t%.15lf\n",atom[i].fx,atom[i].fy,atom[i].fz);
      
      numforce(atom,nForce,&ff,&ener,&simulCond,&box,4,1.e-4);
      
//     Write numerical forces into numfrc file.
      
      fprintf(numfrc,"#step=%d\n",simulCond.step);
      for(i=0;i<simulCond.natom;i++)
      {
	fprintf(numfrc,"%.15lf\t%.15lf\t%.15lf\n",nForce[i].x,nForce[i].y,nForce[i].z);
	
	atom[i].fx=nForce[i].x;
	atom[i].fy=nForce[i].y;
	atom[i].fz=nForce[i].z;
      }      
    }
  }
  
  if(simulCond.keymd)
  {
//   Molecular dynamics loop starts here.

    for(simulCond.step=1;simulCond.step<=simulCond.nsteps;simulCond.step++)
    {
      
//     Integration of the Newtonian equations. First stage
//     of Velocity Verlet algorithm.

      if(simulCond.integrator==1)
      {
	vv_integrate(atom,&ener,&simulCond,constList,&box,1);
// 	printf("Velocity verlet fisrt stage done for step %d\n",simulCond.step);
      }
      
//     List update if needed.
      
      makelist(&simulCond,atom,&ff,constList,&box);
      
//     Energies calculation.

      energy(atom,&ff,&ener,&simulCond,&box);
//       printf("Energy done for step %d\n",simulCond.step);
      
//     Numerical derivatives to estimate forces.

      if(simulCond.numDeriv==1)
      {
	
	numforce(atom,nForce,&ff,&ener,&simulCond,&box,4,1.e-4);
	
//       Write analytical forces into frc file.
	
	fprintf(frc,"#step=%d\n",simulCond.step);
	for(i=0;i<simulCond.natom;i++)
	  fprintf(frc,"%.15lf\t%.15lf\t%.15lf\n",atom[i].fx,atom[i].fy,atom[i].fz);
	  
//       Write numerical forces into numfrc file.
	
	fprintf(numfrc,"#step=%d\n",simulCond.step);
      
	for(i=0;i<simulCond.natom;i++)
	{
	  fprintf(numfrc,"%.15lf\t%.15lf\t%.15lf\n",nForce[i].x,nForce[i].y,nForce[i].z);
	  
// 	Overwrite analytical forces with numerical forces.
	  atom[i].fx=nForce[i].x;
	  atom[i].fy=nForce[i].y;
	  atom[i].fz=nForce[i].z;
	}
	
      }
      
//     Integration of the Newtonian equations. Leapfrog or
//     second stage of Velocity Verlet algorithm.
      
      if(simulCond.integrator==0)
      {
	lf_integrate(atom,&ener,&simulCond,constList,&box);
// 	printf("Leap frog done for step %d\n",simulCond.step);
      }
      else if(simulCond.integrator==1)
      {
	vv_integrate(atom,&ener,&simulCond,constList,&box,2);
// 	printf("Velocity verlet second stage done for step %d\n",simulCond.step);
      }
      
      ener.tot=ener.kin+ener.pot;
      
//     Write thermodynamics properties into ener file.
      
      if( (simulCond.keyener) && (simulCond.step%simulCond.printo==0) )
	fprintf(fener,_ENERFORMAT_,
	  simulCond.step,ener.tot/kcaltoiu,ener.kin/kcaltoiu,ener.pot/kcaltoiu,
	  ener.elec/kcaltoiu,ener.vdw/kcaltoiu,ener.bond/kcaltoiu,
	  ener.ang/kcaltoiu,ener.ub/kcaltoiu,ener.dihe/kcaltoiu,
	  ener.impr/kcaltoiu);
	
//     Write coordinates into traj file.
      
      if( (simulCond.keytraj) && (simulCond.step%simulCond.printtr==0) )
      {
	fprintf(traj,"%d\n",simulCond.natom);
	fprintf(traj,"step %d\n",simulCond.step);
	for(i=0;i<simulCond.natom;i++)
	  fprintf(traj,"%s\t%7.3lf\t%7.3lf\t%7.3lf\n",atom[i].label,atom[i].x,atom[i].y,atom[i].z);
      }
      
    }
  
//   Molecular dynamics loop ends here.
  }

  if(simulCond.keyener)
    fclose(fener);
  
  if(simulCond.keytraj)
    fclose(traj);
  
  printf("Normal Termination of MDBas.\n");
  printf("Veni, Vedi, Vici.\n");
  
#if ( defined __linux__ || defined __FreeBSD__ )
//compatible with some unixes-like OS: the struct rusage communicates with the kernel directly.
  struct rusage infos_usage;
  getrusage(RUSAGE_SELF,&infos_usage);
  fprintf(stdout,"Memory used in kBytes is : %ld\n",infos_usage.ru_maxrss);
  fprintf(stdout,"Execution time in Seconds : %lf\n",(double)infos_usage.ru_utime.tv_sec+infos_usage.ru_utime.tv_usec/1000000.0);
#endif /* unixes part */
  
  return 0;
}
