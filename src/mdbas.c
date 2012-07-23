/***************************************************************
 * cmenerve is designed to compute the main energy terms for MM
 * simulations. The goal is to provide a toolkit of routines that
 * can be easily modified and combined.
 * 
 * P.-A. Cazade 02 May 2012
 * ************************************************************/

#include <stdio.h>
#include <time.h>
#include "global.h"
#include "energy.h"
#include "io.h"
#include "rand.h"
#include "utils.h"
#include "integrate.h"
#include "numderiv.h"


#if ( defined __linux__ || defined __FreeBSD__ )
#include <sys/resource.h>
#endif /* for some unixes */

int main(int argc, char* argv[])
{
  
//   void (*test)(int *par1,double *par2,...);
//   test=&(elec);

  FILE *ener,*traj,*frc,*numfrc;
  INPUTS inp;
  ATOM atom;
  FORCEFIELD ff;
  ENERGYFORCE enerFor;
  SIMULPARAMS simulCond;
  CONSTRAINT *constList;
  int i;
  
  init_rand(time(NULL));
  
  read_SIMU(&simulCond,&ff);
  
  printf("SIMU file read\n");
  
  if(simulCond.keyrand)
    init_rand(simulCond.seed);
  
  if(simulCond.mdNature==1)
    simulCond.chargeConst=chgcharmm*kcaltoiu;
  else if(simulCond.mdNature==2)
    simulCond.chargeConst=chgnamd*kcaltoiu;
  else
    simulCond.chargeConst=mu0*X2(clight)*X2(elemchg)*NA*0.1/(angstr);
  
  printf("pointer adress=%p pointer adress=%p\n",constList,&constList);
  read_PSF(&inp,&atom,&ff,&simulCond,constList);
  printf("test pointer %d %d\n",constList[1].a,constList[1].b);
  
  printf("PSF file read\n");
  
  read_TOP(&inp);
  
  printf("TOP file read\n");
  
  read_PAR(&inp);
  
  printf("PAR file read\n");
  
  read_CONF(&atom);
  
  printf("CONF file read\n");
  
  setup(&inp,&atom,&ff,&simulCond,constList);
  
  printf("Setup done\n");
  
  if(simulCond.keyforf)
    write_FORF(&inp,&atom,&ff,&simulCond);
  
  free_temp_array(&inp);
  
  printf("Free temp array done\n");
  
  init_vel(&atom,&simulCond);
  
  //test debug
  for(i=0;i<atom.natom;i++)
  {
    atom.vx[i]=0.;
    atom.vy[i]=0.;
    atom.vz[i]=0.;
  }
  
  if(simulCond.keyener)
  {
    ener=fopen("ener.dat","w");
    
    fprintf(ener,"Step\tEtot\tEkin\tEpot\tEcoul\tEvdw\tEbond\tEangle\tEub\tEdihe\tEimpr\n");
  }
  
  if(simulCond.keytraj)
  {
    traj=fopen("traj.xyz","w");
  }
  
  if(simulCond.integrator==1)
  {

//   Computes kinetic energy at time=0

    enerFor.energyKin=kinetic(&atom);
  
//   Computes potential energies and forces at time=0
  
    energy(&atom,&ff,&enerFor,&simulCond);
    
    enerFor.energyTot=enerFor.energyKin+enerFor.energyPot;
    
    if(simulCond.keyener)
    { 
      fprintf(ener,"%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
	simulCond.step,enerFor.energyTot/kcaltoiu,enerFor.energyKin/kcaltoiu,enerFor.energyPot/kcaltoiu,
	enerFor.energyElec/kcaltoiu,enerFor.energyVdw/kcaltoiu,enerFor.energyBond/kcaltoiu,
	enerFor.energyAng/kcaltoiu,enerFor.energyUb/kcaltoiu,enerFor.energyDih/kcaltoiu,
	enerFor.energyImpr/kcaltoiu);
    }
    
    if(simulCond.keytraj)
    {
      fprintf(traj,"%d\n",atom.natom);
      fprintf(traj,"step %d\n",simulCond.step);
      for(i=0;i<atom.natom;i++)
	  fprintf(traj,"%s\t%7.3lf\t%7.3lf\t%7.3lf\n",atom.atomLabel[i],atom.x[i],atom.y[i],atom.z[i]);
    }
//   Numerical derivatives to estimate forces for initial configuration.
      
    if(simulCond.numDeriv==1)
    {
      frc=fopen("force.dat","w");
      numfrc=fopen("numforce.dat","w");
    
//     Write analytical forces into frc file.
    
      fprintf(frc,"#step=%d\n",simulCond.step);
      for(i=0;i<atom.natom;i++)
	fprintf(frc,"%.15lf\t%.15lf\t%.15lf\n",atom.fx[i],atom.fy[i],atom.fz[i]);
      
      numforce(&atom,&ff,&enerFor,&simulCond,4,1.e-4);
      
//     Write numerical forces into numfrc file.
      
      fprintf(numfrc,"#step=%d\n",simulCond.step);
      for(i=0;i<atom.natom;i++)
      {
	fprintf(numfrc,"%.15lf\t%.15lf\t%.15lf\n",atom.numFx[i],atom.numFy[i],atom.numFz[i]);
	
	atom.fx[i]=atom.numFx[i];
	atom.fy[i]=atom.numFy[i];
	atom.fz[i]=atom.numFz[i];
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
	vv_nve(&atom,&enerFor,&simulCond,1);
	printf("Velocity verlet fisrt stage done for step %d\n",simulCond.step);
      }
      
//     Energies calculation.

      energy(&atom,&ff,&enerFor,&simulCond);
      printf("Energy done for step %d\n",simulCond.step);
      
//     Numerical derivatives to estimate forces.

      if(simulCond.numDeriv==1)
      {
	
	numforce(&atom,&ff,&enerFor,&simulCond,4,1.e-4);
	
//       Write analytical forces into frc file.
	
	fprintf(frc,"#step=%d\n",simulCond.step);
	for(i=0;i<atom.natom;i++)
	  fprintf(frc,"%.15lf\t%.15lf\t%.15lf\n",atom.fx[i],atom.fy[i],atom.fz[i]);
	  
//       Write numerical forces into numfrc file.
	
	fprintf(numfrc,"#step=%d\n",simulCond.step);
      
	for(i=0;i<atom.natom;i++)
	{
	  fprintf(numfrc,"%.15lf\t%.15lf\t%.15lf\n",atom.numFx[i],atom.numFy[i],atom.numFz[i]);
	  
// 	Overwrite analytical forces with numerical forces.
	  
	  atom.fx[i]=atom.numFx[i];
	  atom.fy[i]=atom.numFy[i];
	  atom.fz[i]=atom.numFz[i];
	}
	
      }
      
//     Integration of the Newtonian equations. Leapfrog or
//     second stage of Velocity Verlet algorithm.
      
      if(simulCond.integrator==0)
      {
	lf_nve(&atom,&enerFor,&simulCond,constList);
	printf("Leap frog done for step %d\n",simulCond.step);
      }
      else if(simulCond.integrator==1)
      {
	vv_nve(&atom,&enerFor,&simulCond,2);
	printf("Velocity verlet second stage done for step %d\n",simulCond.step);
      }
      
      enerFor.energyTot=enerFor.energyKin+enerFor.energyPot;
      
//     Write thermodynamics properties into ener file.
      
      if( (simulCond.keyener) && (simulCond.step%simulCond.printo==0) )
	fprintf(ener,"%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
	  simulCond.step,enerFor.energyTot/kcaltoiu,enerFor.energyKin/kcaltoiu,enerFor.energyPot/kcaltoiu,
	  enerFor.energyElec/kcaltoiu,enerFor.energyVdw/kcaltoiu,enerFor.energyBond/kcaltoiu,
	  enerFor.energyAng/kcaltoiu,enerFor.energyUb/kcaltoiu,enerFor.energyDih/kcaltoiu,
	  enerFor.energyImpr/kcaltoiu);
	
//     Write coordinates into traj file.
      
      if( (simulCond.keytraj) && (simulCond.step%simulCond.printtr==0) )
      {
	fprintf(traj,"%d\n",atom.natom);
	fprintf(traj,"step %d\n",simulCond.step);
	for(i=0;i<atom.natom;i++)
	  fprintf(traj,"%s\t%7.3lf\t%7.3lf\t%7.3lf\n",atom.atomLabel[i],atom.x[i],atom.y[i],atom.z[i]);
      }
      
    }
  
//   Molecular dynamics loop ends here.
  }

  if(simulCond.keyener)
    fclose(ener);
  
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
