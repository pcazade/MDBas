/*
 * Copyright (c) 2013 Pierre-Andre Cazade
 * Copyright (c) 2013 Florent hedin
 * 
 * This file is part of MDBas.
 *
 * MDBas is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MDBas is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MDBas.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * \file mdbas.c
 * \brief Main entry file for the program MDBas.
 * \author Pierre-Andre Cazade and Florent Hedin
 * \version alpha-branch
 * \date 2012
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "global.h"
#include "energy.h"
#include "init.h"
#include "io.h"
#include "rand.h"
#include "utils.h"
#include "integrate.h"
#include "numderiv.h"
#include "list.h"
#include "minim.h"
#include "memory.h"
#include "user.h"

#if (defined TIMING && defined __unix__ && !defined __STRICT_ANSI__)
#define TIMER
#include "timing.h"
#endif

#ifdef __unix__
/* Informations about resources usage */
#include <sys/resource.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#define _ENERFORMAT_ "%13d\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\n%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\n%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\n"
#define _ENERLABELS_ "%13s\t%13s\t%13s\t%13s\t%13s\n%13s\t%13s\t%13s\t%13s\t%13s\n%13s\t%13s\t%13s\t%13s\t%13s\n"

FILE *outFile=NULL;

/**
 * \param argc Number of arguments of the command line including the name of the program, so argc >= 1.
 * \param argv Array of char strings containing parameters of the command line ; contains at least the name of the program at index 0.
 *
 * \return On return EXIT_SUCCESS if the program terminates properly.
 */
int main(int argc, char* argv[])
{
  /** Beginning of structures declaration. */
  
#ifdef TIMER
  init_timers();
  create_new_timer(TIMER_ALL);
  update_timer_begin(TIMER_ALL,__func__);
#endif
  
  IO inout;

  NEIGH neigh;
  
  CTRL ctrl;
  PARAM param;
  
  ENERGY ener;
  
  PBC box;
  
  BATH bath;
  
  EWALD ewald;
  
  ATOM *atom=NULL;
  
  BOND *bond=NULL,*ub=NULL;
  ANGLE *angle=NULL;
  DIHE *dihe=NULL,*impr=NULL;
  
  CONSTRAINT *constList=NULL;
  
  DELTA *nForce=NULL;
  
  PARALLEL *parallel=NULL;
  
  double *x,*y,*z;
  double *vx,*vy,*vz;
  double *fx,*fy,*fz;
  double *q,*mass,*rmass;
  double *eps,*sig,*eps14,*sig14;
  
  int *frozen,*nAtConst;
  int *neighOrder,*neighList,*neighPair,*neighList14;
  int **exclList,*exclPair;
  
  /** End of structures declarations. */
  
  /** Beginning of variables declaration. */
  
  char enerLabel[15][7]={ {"Step"} , {"Time"} , {"Temp"} , {"Press"} , {"Vol"} ,
			  {"Etot"} , {"Ekin"}  , {"Epot"} , {"Ecoul"} , {"Evdw"} ,
			  {"Ebond"} , {"Eangle"} , {"Eub"} , {"Edihe"} , {"Eimpr"} };
			  
  char dashes[81]={0};
  
  {
    int k;
    for(k=0;k<80;k++)
      dashes[k]='-';
    
    dashes[80]='\0';
  }
  
  int i,nPrint=0;
  double temp,press;
  
  #ifdef _OPENMP
  int num_threads=1;
  #pragma omp parallel
  {num_threads = omp_get_num_threads();}
  #endif
  
  /** End of variables declaration. */
  
  /** Initialization of the simulation starts here. */
  
  init_system(argc,argv,&inout,&ctrl,&param,&ener,&bath,&neigh,&ewald,&box,
	      &atom,&constList,&bond,&angle,&dihe,&impr,&ub,&x,&y,&z,
	      &vx,&vy,&vz,&fx,&fy,&fz,&mass,&rmass,&q,&eps,&sig,&eps14,
	      &sig14,&frozen,&nAtConst,&neighList,&neighPair,&neighOrder,
	      &neighList14,&exclList,&exclPair);
  
  #ifdef _OPENMP
  fprintf(outFile,"Multi-threading enabled : number of threads = %d\n\n",num_threads);
  #else
  fprintf(outFile,"Multi-threading disabled.\n\n");
  #endif

//  UserEnergyPtr userPtr = NULL;
//  userPtr = loadUserPlugin("user_functions.so","MyEnergyFunction");
//  userPtr();
  
  /** Initialization of the simulation ends here. */
  
  /** Minimization procedure starts here. */
  
  if(ctrl.keyMinim)
  {
    steepestDescent(&ctrl,&param,&ener,&box,&neigh,atom,bond,ub,angle,dihe,impr,x,y,z,fx,fy,fz);
    
    makelist(&ctrl,&param,&box,&neigh,constList,bond,angle,dihe,impr,x,y,z,frozen,
	     &neighList,&neighPair,&neighOrder,&neighList14,&exclList,&exclPair);
    
    init_vel(&param,&box,constList,x,y,z,vx,vy,vz,mass,rmass,frozen,nAtConst);
  }
  
  /** Minimization procedure ends here. */
  
  /** First calculation of the enregy starts here.
   * 
   *  Applies if VV integrator is used or if a single energy is required.
   * 
   */
  
  if(ctrl.integrator==VELOCITY || !ctrl.keyMd)
  {

   /** Computes kinetic energy at time=0. */

    ener.kin=kinetic(&param,vx,vy,vz,mass);
  
  /** Computes potential energies and forces at time=0. */
    
    energy(&ctrl,&param,&ener,&ewald,&box,&neigh,bond,ub,angle,dihe,impr,
	   x,y,z,vx,vy,vz,fx,fy,fz,q,eps,sig,eps14,sig14,frozen,
	   neighList,neighPair,neighOrder,neighList14,exclList,exclPair);
    
    ener.tot=ener.kin+ener.pot;
    
    ener.virtot=ener.virbond+ener.virub+ener.virelec+ener.virvdw+ener.virshake;
    
    if( (ctrl.keyProp) && (param.step%ctrl.printProp==0) )
    {
      write_prop(&inout,&param,&ener,&box);
    }
    
    if( (param.step%ctrl.printOut==0) )
    {
      temp=2.*ener.kin/((double)param.nDegFree*rboltzui);
      if(box.type>0)
	press=(2.*ener.kin-ener.virtot)/(3.*box.vol*bartoiu);
      else
	press=0.;
      
      if( (nPrint%5==0) )
      {
	fprintf(outFile,"\n");
	fprintf(outFile,"%s\n\n",dashes);
	
	fprintf(outFile,_ENERLABELS_,enerLabel[0],enerLabel[1],enerLabel[2],enerLabel[3],enerLabel[4],
			       enerLabel[5],enerLabel[6],enerLabel[7],enerLabel[8],enerLabel[9],
			       enerLabel[10],enerLabel[11],enerLabel[12],enerLabel[13],enerLabel[14]);
	
	fprintf(outFile,"\n%s\n",dashes);
	
	nPrint=0;
      }
      
      fprintf(outFile,_ENERFORMAT_,
	param.step,param.step*param.timeStep,temp,press,box.vol,
	ener.tot/kcaltoiu,ener.kin/kcaltoiu,ener.pot/kcaltoiu,ener.elec/kcaltoiu,ener.vdw/kcaltoiu,
	ener.bond/kcaltoiu,ener.ang/kcaltoiu,ener.ub/kcaltoiu,ener.dihe/kcaltoiu,ener.impr/kcaltoiu);
      
      fprintf(outFile,"%s\n",dashes);
      
      nPrint++;
    }

  /** Numerical derivatives to estimate forces for initial configuration. */
      
    if(ctrl.keyNumForce==1)
    {
      nForce=(DELTA*)my_malloc(param.nAtom*sizeof(*nForce));
      
      numforce(&ctrl,&param,&ener,&box,&neigh,bond,ub,angle,dihe,impr,nForce,x,y,z,4,1.e-4);
      
      for(i=0;i<param.nAtom;i++)
      {
	fx[i]=nForce[i].x;
	fy[i]=nForce[i].y;
	fz[i]=nForce[i].z;
      } 
    }
  }
  /** First calculation of the enregy ends here. */
  
  
  /** MD starts here if required. */
  
  if(ctrl.keyMd)
  {
  /** Molecular dynamics loop starts here.*/

    for(++(param.step);param.step<=param.nSteps;param.step++)
    {
      
      ener.consv=0.;
      
    /** Integration of the Newtonian equations.
     * First stage of Velocity Verlet algorithm.
     */

      if(ctrl.integrator==VELOCITY)
      {
	vv_integrate(&ctrl,&param,&ener,&box,&bath,constList,
		     x,y,z,vx,vy,vz,fx,fy,fz,mass,rmass,nAtConst,1);
      }
      
    /** List update if needed. */
    
      makelist(&ctrl,&param,&box,&neigh,constList,bond,angle,dihe,impr,x,y,z,frozen,
	       &neighList,&neighPair,&neighOrder,&neighList14,&exclList,&exclPair);
      
      
    /** Energies calculation. */
    
      energy(&ctrl,&param,&ener,&ewald,&box,&neigh,bond,ub,angle,dihe,impr,
	     x,y,z,vx,vy,vz,fx,fy,fz,q,eps,sig,eps14,sig14,frozen,
	     neighList,neighPair,neighOrder,neighList14,exclList,exclPair);
      
    /** Numerical derivatives to estimate forces. */

      if(ctrl.keyNumForce==1)
      {
	
	numforce(&ctrl,&param,&ener,&box,&neigh,bond,ub,angle,dihe,impr,nForce,x,y,z,4,1.e-4);
	  
	for(i=0;i<param.nAtom;i++)
	{
	/** Overwrite analytical forces with numerical forces. */
	  fx[i]=nForce[i].x;
	  fy[i]=nForce[i].y;
	  fz[i]=nForce[i].z;
	}
	
      }
      
    /** Integration of the Newtonian equations.
     * Leapfrog or second stage of Velocity Verlet algorithm.
     */
      
      if(ctrl.integrator==LEAPFROG)
      {
	lf_integrate(&ctrl,&param,&ener,&box,&bath,constList,
		     x,y,z,vx,vy,vz,fx,fy,fz,mass,rmass,nAtConst);
      }
      else if(ctrl.integrator==VELOCITY)
      {
	vv_integrate(&ctrl,&param,&ener,&box,&bath,constList,
		     x,y,z,vx,vy,vz,fx,fy,fz,mass,rmass,nAtConst,2);
      }
      
      ener.tot=ener.kin+ener.pot;
      
      ener.virtot=ener.virbond+ener.virub+ener.virelec+ener.virvdw+ener.virshake;
      
    /** Writes system properties with full precision into PROP file. */
      
      if( (ctrl.keyProp) && (param.step%ctrl.printProp==0) )
      {
	write_prop(&inout,&param,&ener,&box);
      }
      
      if( (param.step%ctrl.printOut==0) )
      {
	temp=2.*ener.kin/((double)param.nDegFree*rboltzui);
	if(box.type>0)
	  press=(2.*ener.kin-ener.virtot)/(3.*box.vol*bartoiu);
	else
	  press=0.;
	
	if( (nPrint%5==0) )
	{
	  fprintf(outFile,"\n");
	  fprintf(outFile,"%s\n\n",dashes);
	  
	  fprintf(outFile,_ENERLABELS_,enerLabel[0],enerLabel[1],enerLabel[2],enerLabel[3],enerLabel[4],
				enerLabel[5],enerLabel[6],enerLabel[7],enerLabel[8],enerLabel[9],
				enerLabel[10],enerLabel[11],enerLabel[12],enerLabel[13],enerLabel[14]);
	  
	  fprintf(outFile,"\n%s\n",dashes);
	  
	  nPrint=0;
	}
	
	fprintf(outFile,_ENERFORMAT_,param.step,param.step*param.timeStep,temp,press,box.vol,
	ener.tot/kcaltoiu,ener.kin/kcaltoiu,ener.pot/kcaltoiu,ener.elec/kcaltoiu,ener.vdw/kcaltoiu,
	ener.bond/kcaltoiu,ener.ang/kcaltoiu,ener.ub/kcaltoiu,ener.dihe/kcaltoiu,ener.impr/kcaltoiu);
	
	fprintf(outFile,"%s\n",dashes);
	
	nPrint++;
      }
	
    /** Writes coordinates into DCD file. */
      
      if( (ctrl.keyTraj) && (param.step%ctrl.printTraj==0) )
      {
	write_DCD_traj(&inout,&param,&box,x,y,z,frozen);
      }
      
    /** Writes restart files. */
      
      if(param.step%ctrl.printRest==0)
      {
	write_CONF(&inout,&param,atom,x,y,z);
// 	fprintf(outFile,"\n%s file written\n",inout.rconName);
	
	write_rest(&inout,&param,&ener,&bath,atom,x,y,z,vx,vy,vz,fx,fy,fz);
// 	fprintf(outFile,"\n%s file written\n",inout.restName);
      }
      
    }
  
  /** Molecular dynamics loop ends here. */
  
  }
  
  /** MD ends here. */
  
  /** Writes restart files. */
    
  write_CONF(&inout,&param,atom,x,y,z);
  fprintf(outFile,"\n%s file written\n",inout.rconName);
  
  write_rest(&inout,&param,&ener,&bath,atom,x,y,z,vx,vy,vz,fx,fy,fz);
  fprintf(outFile,"\n%s file written\n",inout.restName);
  
  fprintf(outFile,"\nNormal Termination of MDBas.\n");
  
#ifdef TIMER
  update_timer_end(TIMER_ALL,__func__);
  /** Write timings **/
  print_timers();
#endif
  
  /** Freeing arrays **/
  free_all(&ctrl,&param,&ewald,&atom,&constList,&bond,&angle,&dihe,&impr,&ub,
	   &x,&y,&z,&vx,&vy,&vz,&fx,&fy,&fz,&mass,&rmass,&q,&eps,&sig,&eps14,
	   &sig14,&frozen,&nAtConst,&neighList,&neighPair,&neighOrder,
	   &neighList14,&exclList,&exclPair);
  
#ifdef __unix__
  struct rusage infos_usage;
  getrusage(RUSAGE_SELF,&infos_usage);
  /** when using omp this time is not correct **/
//  fprintf(outFile,"Execution time in Seconds : %lf\n",(double)(infos_usage.ru_utime.tv_sec+infos_usage.ru_utime.tv_usec/1000000.0));
  fprintf(outFile,"Max amount of physical memory used (kBytes) : %ld\n",infos_usage.ru_maxrss);
#endif
  
  fclose(outFile);
  
  return EXIT_SUCCESS;
}
