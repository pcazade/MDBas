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
 * \file io.c
 * \brief Contains functions in charge of I/O and parsing.
 * \author Pierre-Andre Cazade and Florent Hedin
 * \version alpha-branch
 * \date 2012
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <float.h>

#include "global.h"
#include "io.h"
#include "memory.h"
#include "utils.h"
#include "errors.h"

#ifdef MPI_VERSION
#include "parallel.h"
#else
#include "serial.h"
#endif

/** Pointer to the output file. **/
extern FILE *outFile;

/** First write of traj in dcd, or not*/
static int first=1;

int read_command_line(int *argc, char ***argv,IO *inout,PARALLEL *parallel)
{
  char outName[FINAMELEN];
  
  int i=1;
  
  int err=0;
  
  outFile=NULL;
  
  if(parallel->idProc==0)
  {
    
    strcpy(outName,"OUTPUT");
    
    strcpy(inout->simuName,"SIMU");
    
    while(i<*argc)
    {
      if(!strcmp((*argv)[i],"-i"))
      {
	strcpy(inout->simuName,(*argv)[++i]);
      }
      else if (!strcmp((*argv)[i],"-o"))
      {
	strcpy(outName,(*argv)[++i]);
      }
      else if (!strcmp((*argv)[i],"--help"))
      {
	printf("%s [-i input_file] [-o output_file] [--help]\n",(*argv)[0]);
	err=1;
	break;
      }
      else
      {
	err=2;
	break;
      }
      i++;
    }
    
    if(err==0)
    {
      outFile=fopen(outName,"w");
      if(outFile==NULL)
      {
	outFile=stdout;
	my_error(UNKNOWN_GENERAL_ERROR,__FILE__,__LINE__,0);
	
      } 
    }
    else
    {
      outFile=stdout;
    }
  
  }
  
  bcast_int_para(&err,1,0);
  
  return(err);
}

void read_SIMU(IO *inout,CTRL *ctrl,PARAM *param,BATH *bath,NEIGH *neigh,EWALD *ewald,PBC *box)
{
    char buff1[1024]="", *buff2=NULL, *buff3=NULL, *buff4=NULL;

    inout->simuFile=fopen(inout->simuName,"r");

    if (inout->simuFile==NULL)
    {
        my_error(SIMU_FILE_ERROR,__FILE__,__LINE__,0);
    }

    strcpy(inout->confName,"CONF");
    strcpy(inout->forfName,"FORF");
    strcpy(inout->propName,"PROP");
    strcpy(inout->restName,"REST");
    strcpy(inout->rconName,"RCON");
    strcpy(inout->trajName,"TRAJ");

    while(fgets(buff1,1024,inout->simuFile)!=NULL)
    {
        buff2=strtok(buff1," \n\t");

        if(buff2==NULL)
            continue;

        nocase(buff2);

        if(!strcmp(buff2,"mdbas"))
            ctrl->mdType=0;

        else if(!strcmp(buff2,"charmm"))
            ctrl->mdType=1;

        else if(!strcmp(buff2,"namd"))
            ctrl->mdType=2;

        else if(!strcmp(buff2,"dlpoly"))
            ctrl->mdType=3;

        else if(!strcmp(buff2,"nomd"))
            ctrl->keyMd=0;

        else if(!strcmp(buff2,"restart"))
            ctrl->keyRest=1;

        else if(!strcmp(buff2,"file"))
        {
            buff3=strtok(NULL," \n\t");
            if(buff3==NULL)
                my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,2,buff2,"type of file");

            nocase(buff3);

            if(!strcmp(buff3,"conf"))
            {
                buff4=strtok(NULL," \n\t");
                if(buff4==NULL)
                    my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,2,buff3,"name of file");

                strcpy(inout->confName,buff4);
            }
            else if(!strcmp(buff3,"forf"))
            {
                buff4=strtok(NULL," \n\t");
                if(buff4==NULL)
                    my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,2,buff3,"name of file");

                strcpy(inout->forfName,buff4);
            }
            else if(!strcmp(buff3,"prop"))
            {
                buff4=strtok(NULL," \n\t");
                if(buff4==NULL)
                    my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,2,buff3,"name of file");

                strcpy(inout->propName,buff4);
            }
            else if(!strcmp(buff3,"rest"))
            {
                buff4=strtok(NULL," \n\t");
                if(buff4==NULL)
                    my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,2,buff3,"name of file");

                strcpy(inout->restName,buff4);
            }
            else if(!strcmp(buff3,"rcon"))
            {
                buff4=strtok(NULL," \n\t");
                if(buff4==NULL)
                    my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,2,buff3,"name of file");

                strcpy(inout->rconName,buff4);
            }
            else if(!strcmp(buff3,"traj"))
            {
                buff4=strtok(NULL," \n\t");
                if(buff4==NULL)
                    my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,2,buff3,"name of file");

                strcpy(inout->trajName,buff4);
            }
            else
                my_error(SIMU_PARAM_ERROR,__FILE__,__LINE__,2,buff2,buff3);
        }
        else if(!strcmp(buff2,"minim"))
        {
            ctrl->keyMinim=1;
            buff3=strtok(NULL," \n\t");
            if(buff3==NULL)
                my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

            nocase(buff3);

            if(!strcmp(buff3,"tol"))
            {
                buff4=strtok(NULL," \n\t");
                if(buff4==NULL)
                    my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

                param->tolMinim=atof(buff4);
            }
            else if(!strcmp(buff3,"maxcycle"))
            {
                buff4=strtok(NULL," \n\t");
                if(buff4==NULL)
                    my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

                param->maxminst=atoi(buff4);
            }
            else if(!strcmp(buff3,"maxsize"))
            {
                buff4=strtok(NULL," \n\t");
                if(buff4==NULL)
                    my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

                param->maxminsiz=atof(buff4);
            }
            else
                my_error(SIMU_PARAM_ERROR,__FILE__,__LINE__,2,buff2,buff3);
        }
        else if(!strcmp(buff2,"timestep"))
        {
            buff3=strtok(NULL," \n\t");
            if(buff3==NULL)
                my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

            param->timeStep=atof(buff3);
        }
        else if(!strcmp(buff2,"nsteps"))
        {
            buff3=strtok(NULL," \n\t");
            if(buff3==NULL)
                my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

            param->nSteps=atoi(buff3);
        }
        else if(!strcmp(buff2,"cutoff"))
        {
            buff3=strtok(NULL," \n\t");
            if(buff3==NULL)
                my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

            param->cutOff=atof(buff3);
        }
        else if(!strcmp(buff2,"cuton"))
        {
            buff3=strtok(NULL," \n\t");
            if(buff3==NULL)
                my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

            param->cutOn=atof(buff3);
        }
        else if(!strcmp(buff2,"delr"))
        {
            buff3=strtok(NULL," \n\t");
            if(buff3==NULL)
                my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

            param->delr=atof(buff3);
        }
        else if(!strcmp(buff2,"ewald"))
        {
            buff3=strtok(NULL," \n\t");
            if(buff3==NULL)
                my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

            nocase(buff3);

            if(!strcmp(buff3,"prec"))
            {
                buff4=strtok(NULL," \n\t");
                if(buff4==NULL)
                    my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

                ewald->prec=atof(buff4);
            }
            else if(!strcmp(buff3,"alpha"))
            {
                buff4=strtok(NULL," \n\t");
                if(buff4==NULL)
                    my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

                ctrl->keyAlpha=1;
                ewald->alpha=atof(buff4);
            }
            else if(!strcmp(buff3,"nbsp"))
            {
                buff4=strtok(NULL," \n\t");
                if(buff4==NULL)
                    my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

                ewald->nbsp=atoi(buff4);
            }
            else if(!strcmp(buff3,"recvec"))
            {
                buff4=strtok(NULL," \n\t");
                if(buff4==NULL)
                    my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

                ctrl->keyMmax=1;
                ewald->m1max=atof(buff4);

                buff4=strtok(NULL," \n\t");
                if(buff4==NULL)
                    my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

                ewald->m2max=atof(buff4);

                buff4=strtok(NULL," \n\t");
                if(buff4==NULL)
                    my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

                ewald->m3max=atof(buff4);
            }
            else
                my_error(SIMU_PARAM_ERROR,__FILE__,__LINE__,2,buff2,buff3);
        }
        else if(!strcmp(buff2,"elec"))
        {
            buff3=strtok(NULL," \n\t");
            if(buff3==NULL)
                my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

            nocase(buff3);

            if(!strcmp(buff3,"noelec"))
                ctrl->elecType=NOELEC;
            else if(!strcmp(buff3,"full"))
                ctrl->elecType=FULL;
            else if(!strcmp(buff3,"shift1"))
                ctrl->elecType=SHIFT1;
            else if(!strcmp(buff3,"shift2"))
                ctrl->elecType=SHIFT2;
            else if(!strcmp(buff3,"switch"))
                ctrl->elecType=SWITCH;
            else if(!strcmp(buff3,"ewald"))
                ctrl->keyEwald=1;
            else if(!strcmp(buff3,"spme"))
                ctrl->keyEwald=2;
            else
                my_error(SIMU_PARAM_ERROR,__FILE__,__LINE__,2,buff2,buff3);
        }
        else if(!strcmp(buff2,"vdw"))
        {
            buff3=strtok(NULL," \n\t");
            if(buff3==NULL)
                my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

            nocase(buff3);

            if(!strcmp(buff3,"novdw"))
                ctrl->vdwType=NOVDW;
            else if(!strcmp(buff3,"full"))
                ctrl->vdwType=VFULL;
            else if(!strcmp(buff3,"switch"))
                ctrl->vdwType=VSWITCH;
            else
                my_error(SIMU_PARAM_ERROR,__FILE__,__LINE__,2,buff2,buff3);
        }
        else if(!strcmp(buff2,"nb14"))
        {
            ctrl->keyNb14=1;
            buff3=strtok(NULL," \n\t");
            if(buff3==NULL)
                my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

            param->scal14=atof(buff3);
        }
        else if(!strcmp(buff2,"numforce"))
        {
            ctrl->keyNumForce=1;
        }
        else if(!strcmp(buff2,"list"))
        {
            buff3=strtok(NULL," \n\t");
            if(buff3==NULL)
                my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

            neigh->update=atoi(buff3);
        }
        else if(!strcmp(buff2,"link"))
        {
            buff3=strtok(NULL," \n\t");
            if(buff3==NULL)
                my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

            neigh->linkRatio=atoi(buff3);
        }
        else if(!strcmp(buff2,"nolink"))
        {
            ctrl->noLink=1;
        }
        else if(!strcmp(buff2,"integrator"))
        {
            buff3=strtok(NULL," \n\t");
            if(buff3==NULL)
                my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

            nocase(buff3);

            if(!strcmp(buff3,"leapfrog"))
                ctrl->integrator=LEAPFROG;
            else if(!strcmp(buff3,"velocity"))
                ctrl->integrator=VELOCITY;
            else
                my_error(SIMU_PARAM_ERROR,__FILE__,__LINE__,2,buff2,buff3);
        }
        else if(!strcmp(buff2,"ensemble"))
        {
            buff3=strtok(NULL," \n\t");
            if(buff3==NULL)
                my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

            nocase(buff3);

            if(!strcmp(buff3,"nve"))
                ctrl->ens=NVE;
            else if(!strcmp(buff3,"nvt"))
            {

                buff4=strtok(NULL," \n\t");
                if(buff4==NULL)
                    my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

                nocase(buff4);

                if(!strcmp(buff4,"beren"))
                {
                    ctrl->ens=NVT_B;
                }
                else if(!strcmp(buff4,"hoover"))
                {
                    ctrl->ens=NVT_H;
                }
                else
                    my_error(SIMU_PARAM_ERROR,__FILE__,__LINE__,2,buff3,buff4);
            }
            else if(!strcmp(buff3,"npt"))
            {

                buff4=strtok(NULL," \n\t");
                if(buff4==NULL)
                    my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

                nocase(buff4);

                if(!strcmp(buff4,"beren"))
                {
                    ctrl->ens=NPT_B;
                }
                else if(!strcmp(buff4,"hoover"))
                {
                    ctrl->ens=NPT_H;
                }
                else
                    my_error(SIMU_PARAM_ERROR,__FILE__,__LINE__,2,buff3,buff4);

            }
            else
                my_error(SIMU_PARAM_ERROR,__FILE__,__LINE__,2,buff2,buff3);
        }
        else if(!strcmp(buff2,"taut"))
        {

            buff3=strtok(NULL," \n\t");
            if(buff3==NULL)
                my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

            bath->tauT=atof(buff3);

        }
        else if(!strcmp(buff2,"taup"))
        {

            buff3=strtok(NULL," \n\t");
            if(buff3==NULL)
                my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

            bath->tauP=atof(buff3);

        }
        else if(!strcmp(buff2,"temperature"))
        {
            buff3=strtok(NULL," \n\t");
            if(buff3==NULL)
                my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

            param->temp0=atof(buff3);
        }
        else if(!strcmp(buff2,"pressure"))
        {
            buff3=strtok(NULL," \n\t");
            if(buff3==NULL)
                my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

            param->press0=atof(buff3)*bartoiu;
        }
        else if(!strcmp(buff2,"compressibility"))
        {
            buff3=strtok(NULL," \n\t");
            if(buff3==NULL)
                my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

            bath->compress=atof(buff3)/bartoiu;
        }
        else if(!strcmp(buff2,"consth"))
        {
            ctrl->keyConstH=1;
        }
        else if(!strcmp(buff2,"shake"))
        {
            buff3=strtok(NULL," \n\t");
            if(buff3==NULL)
                my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

            nocase(buff3);

            if(!strcmp(buff3,"tol"))
            {
                buff4=strtok(NULL," \n\t");
                if(buff4==NULL)
                    my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

                param->tolShake=atof(buff4);
            }
            else if(!strcmp(buff3,"maxcycle"))
            {
                buff4=strtok(NULL," \n\t");
                if(buff4==NULL)
                    my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

                param->maxCycle=atoi(buff4);
            }
            else
                my_error(SIMU_PARAM_ERROR,__FILE__,__LINE__,2,buff2,buff3);
        }
        else if(!strcmp(buff2,"seed"))
        {
            buff3=strtok(NULL," \n\t");
            if(buff3==NULL)
                my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

            ctrl->keyRand=1;
            ctrl->seed=atoi(buff3);
        }
        else if(!strcmp(buff2,"print"))
        {
            buff3=strtok(NULL," \n\t");
            if(buff3==NULL)
                my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

            ctrl->printOut=atoi(buff3);
        }
        else if(!strcmp(buff2,"prop"))
        {
            buff3=strtok(NULL," \n\t");
            if(buff3==NULL)
                my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

            ctrl->keyProp=1;
            ctrl->printProp=atoi(buff3);
        }
        else if(!strcmp(buff2,"traj"))
        {
            buff3=strtok(NULL," \n\t");
            if(buff3==NULL)
                my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

            ctrl->keyTraj=1;
            ctrl->printTraj=atoi(buff3);
        }
        else if(!strcmp(buff2,"resconf"))
        {
            buff3=strtok(NULL," \n\t");
            if(buff3==NULL)
                my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

            ctrl->printRest=atoi(buff3);
        }
        else if(!strcmp(buff2,"write"))
        {
            buff3=strtok(NULL," \n\t");
            if(buff3==NULL)
                my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

            nocase(buff3);

            if(!strcmp(buff3,"field"))
                ctrl->keyForF=1;
            else
                my_error(SIMU_PARAM_ERROR,__FILE__,__LINE__,2,buff2,buff3);
        }
        else if(!strcmp(buff2,"pbc"))
        {
            buff3=strtok(NULL," \n\t");
            if(buff3==NULL)
                my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

            box->type=atoi(buff3);

            if(fgets(buff1,1024,inout->simuFile)!=NULL)
            {

                buff2=strtok(buff1," \n\t");
                if(buff2==NULL)
                    my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);
                else if(isdigit(buff2[0])==0)
                    my_error(SIMU_PARAM_ERROR,__FILE__,__LINE__,0);

                box->a1=atof(buff2);

                buff2=strtok(NULL," \n\t");
                if(buff2==NULL)
                    my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);
                else if(isdigit(buff2[0])==0)
                    my_error(SIMU_PARAM_ERROR,__FILE__,__LINE__,0);

                box->a2=atof(buff2);

                buff2=strtok(NULL," \n\t");
                if(buff2==NULL)
                    my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);
                else if(isdigit(buff2[0])==0)
                    my_error(SIMU_PARAM_ERROR,__FILE__,__LINE__,0);

                box->a3=atof(buff2);

            }
            else
                my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

            if(fgets(buff1,1024,inout->simuFile)!=NULL)
            {

                buff2=strtok(buff1," \n\t");
                if(buff2==NULL)
                    my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);
                else if(isdigit(buff2[0])==0)
                    my_error(SIMU_PARAM_ERROR,__FILE__,__LINE__,0);

                box->b1=atof(buff2);

                buff2=strtok(NULL," \n\t");
                if(buff2==NULL)
                    my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);
                else if(isdigit(buff2[0])==0)
                    my_error(SIMU_PARAM_ERROR,__FILE__,__LINE__,0);

                box->b2=atof(buff2);

                buff2=strtok(NULL," \n\t");
                if(buff2==NULL)
                    my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);
                else if(isdigit(buff2[0])==0)
                    my_error(SIMU_PARAM_ERROR,__FILE__,__LINE__,0);

                box->b3=atof(buff2);

            }
            else
                my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);

            if(fgets(buff1,1024,inout->simuFile)!=NULL)
            {

                buff2=strtok(buff1," \n\t");
                if(buff2==NULL)
                    my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);
                else if(isdigit(buff2[0])==0)
                    my_error(SIMU_PARAM_ERROR,__FILE__,__LINE__,0);

                box->c1=atof(buff2);

                buff2=strtok(NULL," \n\t");
                if(buff2==NULL)
                    my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);
                else if(isdigit(buff2[0])==0)
                    my_error(SIMU_PARAM_ERROR,__FILE__,__LINE__,0);

                box->c2=atof(buff2);

                buff2=strtok(NULL," \n\t");
                if(buff2==NULL)
                    my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);
                else if(isdigit(buff2[0])==0)
                    my_error(SIMU_PARAM_ERROR,__FILE__,__LINE__,0);

                box->c3=atof(buff2);

            }
            else
                my_error(SIMU_NOTFOUND_ERROR,__FILE__,__LINE__,0);
        }
        else if(!strcmp(buff2,"end"))
            break;
        else
            my_error(SIMU_KEYWORD_ERROR,__FILE__,__LINE__,1,buff2);
    }

    param->rcutOff=1./param->cutOff;

    param->cutOff2=X2(param->cutOff);
    param->rcutOff2=1./param->cutOff2;

    param->cutOn2=X2(param->cutOn);

    param->switch2=1./X3(param->cutOff2-param->cutOn2);

    param->rTimeStep=1./param->timeStep;


}

void read_CONF(IO *inout,PARAM *param,ATOM **atom,double **x,double **y, double **z)
{

    char buff1[1024]="", *buff2=NULL;

    char ren[5]="",atl[5]="",sen[5]="";
    int i,/*j,*/atn,res,ire;
    double wei,xx,yy,zz;

    inout->confFile=fopen(inout->confName,"r");

    if (inout->confFile==NULL)
    {
        my_error(CONF_FILE_ERROR,__FILE__,__LINE__,0);
    }

    while(fgets(buff1,1024,inout->confFile)!=NULL)
    {

        if(buff1[0]!='*')
            break;

    }

    buff2=strtok(buff1," \n\t");
    param->nAtom=atoi(buff2);

    *atom=(ATOM*)my_malloc(param->nAtom*sizeof(ATOM));

    *x=(double*)my_malloc(param->nAtom*sizeof(double));
    *y=(double*)my_malloc(param->nAtom*sizeof(double));
    *z=(double*)my_malloc(param->nAtom*sizeof(double));

    for(i=0; i<param->nAtom; i++)
    {
        fscanf(inout->confFile,"%d %d %4s %4s %lf %lf %lf %4s %d %lf",&atn,&ire,ren,atl,&xx,&yy,&zz,sen,&res,&wei);

        (*x)[i]=xx;
        (*y)[i]=yy;
        (*z)[i]=zz;

        (*atom)[i].ires=ire;
        (*atom)[i].resi=res;

        strcpy((*atom)[i].label,atl);
        strcpy((*atom)[i].resn,ren);
        strcpy((*atom)[i].segn,sen);

    }

    fclose(inout->confFile);
}

void read_rest(IO *inout,PARAM *param,ENERGY *ener,BATH *bath,ATOM **atom,
               double **x,double **y, double **z,double **vx,double **vy,
               double **vz,double **fx,double **fy, double **fz)
{

    size_t ret;

    inout->restFile=fopen(inout->restName,"rb");

    ret=fread(&(param->step),sizeof(int),1,inout->restFile);
    if(ret!=1)
        my_error(UNKNOWN_GENERAL_ERROR,__FILE__,__LINE__,0);

    ret=fread(&(param->nAtom),sizeof(int),1,inout->restFile);
    if(ret!=1)
        my_error(UNKNOWN_GENERAL_ERROR,__FILE__,__LINE__,0);

    *atom=(ATOM*)my_malloc(param->nAtom*sizeof(ATOM));

    *x=(double*)my_malloc(param->nAtom*sizeof(double));
    *y=(double*)my_malloc(param->nAtom*sizeof(double));
    *z=(double*)my_malloc(param->nAtom*sizeof(double));

    *vx=(double*)my_malloc(param->nAtom*sizeof(double));
    *vy=(double*)my_malloc(param->nAtom*sizeof(double));
    *vz=(double*)my_malloc(param->nAtom*sizeof(double));

    *fx=(double*)my_malloc(param->nAtom*sizeof(double));
    *fy=(double*)my_malloc(param->nAtom*sizeof(double));
    *fz=(double*)my_malloc(param->nAtom*sizeof(double));

    ret=fread(*atom,sizeof(ATOM),param->nAtom,inout->restFile);
    if(ret!=(size_t)param->nAtom)
        my_error(UNKNOWN_GENERAL_ERROR,__FILE__,__LINE__,0);

    ret=fread(*x,sizeof(double),param->nAtom,inout->restFile);
    if(ret!=(size_t)param->nAtom)
        my_error(UNKNOWN_GENERAL_ERROR,__FILE__,__LINE__,0);

    ret=fread(*y,sizeof(double),param->nAtom,inout->restFile);
    if(ret!=(size_t)param->nAtom)
        my_error(UNKNOWN_GENERAL_ERROR,__FILE__,__LINE__,0);

    ret=fread(*z,sizeof(double),param->nAtom,inout->restFile);
    if(ret!=(size_t)param->nAtom)
        my_error(UNKNOWN_GENERAL_ERROR,__FILE__,__LINE__,0);

    ret=fread(*vx,sizeof(double),param->nAtom,inout->restFile);
    if(ret!=(size_t)param->nAtom)
        my_error(UNKNOWN_GENERAL_ERROR,__FILE__,__LINE__,0);

    ret=fread(*vy,sizeof(double),param->nAtom,inout->restFile);
    if(ret!=(size_t)param->nAtom)
        my_error(UNKNOWN_GENERAL_ERROR,__FILE__,__LINE__,0);

    ret=fread(*vz,sizeof(double),param->nAtom,inout->restFile);
    if(ret!=(size_t)param->nAtom)
        my_error(UNKNOWN_GENERAL_ERROR,__FILE__,__LINE__,0);

    ret=fread(*fx,sizeof(double),param->nAtom,inout->restFile);
    if(ret!=(size_t)param->nAtom)
        my_error(UNKNOWN_GENERAL_ERROR,__FILE__,__LINE__,0);

    ret=fread(*fy,sizeof(double),param->nAtom,inout->restFile);
    if(ret!=(size_t)param->nAtom)
        my_error(UNKNOWN_GENERAL_ERROR,__FILE__,__LINE__,0);

    ret=fread(*fz,sizeof(double),param->nAtom,inout->restFile);
    if(ret!=(size_t)param->nAtom)
        my_error(UNKNOWN_GENERAL_ERROR,__FILE__,__LINE__,0);

    ret=fread(ener,sizeof(ENERGY),1,inout->restFile);
    if(ret!=1)
        my_error(UNKNOWN_GENERAL_ERROR,__FILE__,__LINE__,0);

    ret=fread(&(bath->chiT),sizeof(double),1,inout->restFile);
    if(ret!=1)
        my_error(UNKNOWN_GENERAL_ERROR,__FILE__,__LINE__,0);

    ret=fread(&(bath->chiP),sizeof(double),1,inout->restFile);
    if(ret!=1)
        my_error(UNKNOWN_GENERAL_ERROR,__FILE__,__LINE__,0);

    fclose(inout->restFile);

}

void read_FORF(IO *inout,PARAM *param,ATOM atom[],CONSTRAINT **constList,BOND **bond,
               ANGLE **angle,DIHE **dihe,DIHE **impr,BOND **ub,double *eps,double *sig,
               double *eps14, double *sig14,double *mass,double *q,int *frozen,int *nAtConst)
{

    char buff1[1024]="", *buff2=NULL, buff3[1024]="", *buff4=NULL;
    int i,k,ia,ib;
    int nAtomCheck;

    inout->forfFile=fopen(inout->forfName,"r");

    param->nBond=0;
    param->nConst=0;
    param->nUb=0;
    param->nAngle=0;
    param->nDihedral=0;
    param->nImproper=0;

    while(fgets(buff1,1024,inout->forfFile)!=NULL)
    {
        nocase(buff1);

        if(buff1[0]=='#')
            continue;

        buff2=strtok(buff1," \n\t");

        if(!strcmp(buff2,"atoms"))
        {

            nAtomCheck=atoi(strtok(NULL," \n\t"));
            if(nAtomCheck!=param->nAtom)
                my_error(CONF_ATNUM_ERROR,__FILE__,__LINE__,2,&(param->nAtom),&nAtomCheck);

            k=0;
            while(k<param->nAtom)
            {
                if(fgets(buff3,1024,inout->forfFile)!=NULL)
                {
                    if(buff3[0]=='#')
                        continue;

                    i=atoi(strtok(buff3," \n\t"))-1;
                    if(i!=k)
                        my_error(PSF_FILE_ERROR,__FILE__,__LINE__,0);

                    buff4=strtok(NULL," \n\t");

                    atom[i].type=atoi(strtok(NULL," \n\t"));
                    q[i]=atof(strtok(NULL," \n\t"));
                    mass[i]=atof(strtok(NULL," \n\t"));
                    frozen[i]=atoi(strtok(NULL," \n\t"));

                    nAtConst[i]=0;

                    if( ( frozen[i]!=0 ) && ( frozen[i]!=1 ) )
                        my_error(UNKNOWN_GENERAL_ERROR,__FILE__,__LINE__,0);

                    k++;
                }
                else
                    my_error(PSF_BADLINE_ERROR,__FILE__,__LINE__,0);
            }
        }
        else if(!strcmp(buff2,"bonds"))
        {
            param->nBond=atoi(strtok(NULL," \n\t"));

            if(param->nBond==0)
                continue;

            *bond=(BOND*)my_malloc(param->nBond*sizeof(BOND));

            k=0;
            while(k<param->nBond)
            {
                if(fgets(buff3,1024,inout->forfFile)!=NULL)
                {
                    if(buff3[0]=='#')
                        continue;

                    (*bond)[k].a=atoi(strtok(buff3," \n\t"))-1;
                    (*bond)[k].b=atoi(strtok(NULL," \n\t"))-1;

                    (*bond)[k].type=atoi(strtok(NULL," \n\t"));

                    (*bond)[k].k=atof(strtok(NULL," \n\t"))*kcaltoiu;
                    (*bond)[k].r0=atof(strtok(NULL," \n\t"));
                    (*bond)[k].beta=atof(strtok(NULL," \n\t"));

                    k++;

                }
                else
                    my_error(PSF_BOND_SEQ_ERROR,__FILE__,__LINE__,0);
            }

        }
        else if(!strcmp(buff2,"constraints"))
        {

            param->nConst=atoi(strtok(NULL," \n\t"));

            if(param->nConst==0)
                continue;

            *constList=(CONSTRAINT*)malloc(param->nConst*sizeof(CONSTRAINT));

            k=0;
            while(k<param->nConst)
            {
                if(fgets(buff3,1024,inout->forfFile)!=NULL)
                {
                    if(buff3[0]=='#')
                        continue;

                    ia=atoi(strtok(buff3," \n\t"))-1;
                    ib=atoi(strtok(NULL," \n\t"))-1;

                    (*constList)[k].a=ia;
                    (*constList)[k].b=ib;

                    buff4=strtok(NULL," \n\t");
                    (*constList)[k].rc2=X2(atof(buff4));

                    nAtConst[ia]++;
                    nAtConst[ib]++;

                    k++;

                }
                else
                    my_error(UNKNOWN_GENERAL_ERROR,__FILE__,__LINE__,0);
            }

        }
        else if(!strcmp(buff2,"urey-bradley"))
        {
            param->nUb=atoi(strtok(NULL," \n\t"));

            if(param->nUb==0)
                continue;

            *ub=(BOND*)my_malloc(param->nUb*sizeof(BOND));

            k=0;
            while(k<param->nUb)
            {
                if(fgets(buff3,1024,inout->forfFile)!=NULL)
                {
                    if(buff3[0]=='#')
                        continue;

                    (*ub)[k].a=atoi(strtok(buff3," \n\t"))-1;
                    (*ub)[k].b=atoi(strtok(NULL," \n\t"))-1;

                    (*ub)[k].type=0;

                    (*ub)[k].k=atof(strtok(NULL," \n\t"))*kcaltoiu;
                    (*ub)[k].r0=atof(strtok(NULL," \n\t"));
                    (*ub)[k].beta=atof(strtok(NULL," \n\t"));

                    k++;

                }
                else
                    my_error(UNKNOWN_GENERAL_ERROR,__FILE__,__LINE__,0);
            }

        }
        else if(!strcmp(buff2,"angles"))
        {
            param->nAngle=atoi(strtok(NULL," \n\t"));

            if(param->nAngle==0)
                continue;

            *angle=(ANGLE*)my_malloc(param->nAngle*sizeof(ANGLE));

            k=0;
            while(k<param->nAngle)
            {
                if(fgets(buff3,1024,inout->forfFile)!=NULL)
                {
                    if(buff3[0]=='#')
                        continue;

                    (*angle)[k].a=atoi(strtok(buff3," \n\t"))-1;
                    (*angle)[k].b=atoi(strtok(NULL," \n\t"))-1;
                    (*angle)[k].c=atoi(strtok(NULL," \n\t"))-1;

                    (*angle)[k].type=atoi(strtok(NULL," \n\t"));

                    (*angle)[k].k=atof(strtok(NULL," \n\t"))*kcaltoiu;
                    (*angle)[k].theta0=atof(strtok(NULL," \n\t"))*PI/180.;

                    k++;

                }
                else
                    my_error(PSF_ANGL_SEQ_ERROR,__FILE__,__LINE__,0);
            }

        }
        else if(!strcmp(buff2,"dihedrals"))
        {
            param->nDihedral=atoi(strtok(NULL," \n\t"));

            if(param->nDihedral==0)
                continue;

            *dihe=(DIHE*)my_malloc(param->nDihedral*sizeof(DIHE));

            k=0;
            while(k<param->nDihedral)
            {
                if(fgets(buff3,1024,inout->forfFile)!=NULL)
                {
                    if(buff3[0]=='#')
                        continue;

                    (*dihe)[k].a=atoi(strtok(buff3," \n\t"))-1;
                    (*dihe)[k].b=atoi(strtok(NULL," \n\t"))-1;
                    (*dihe)[k].c=atoi(strtok(NULL," \n\t"))-1;
                    (*dihe)[k].d=atoi(strtok(NULL," \n\t"))-1;

                    (*dihe)[k].type=atoi(strtok(NULL," \n\t"));
                    (*dihe)[k].order=atoi(strtok(NULL," \n\t"));

                    (*dihe)[k].k=atof(strtok(NULL," \n\t"))*kcaltoiu;
                    (*dihe)[k].phi0=atof(strtok(NULL," \n\t"))*PI/180.;
                    (*dihe)[k].mult=atof(strtok(NULL," \n\t"));

                    k++;

                }
                else
                    my_error(PSF_DIHE_SEQ_ERROR,__FILE__,__LINE__,0);
            }

        }
        else if(!strcmp(buff2,"impropers"))
        {
            param->nImproper=atoi(strtok(NULL," \n\t"));

            if(param->nImproper==0)
                continue;

            *impr=(DIHE*)my_malloc(param->nImproper*sizeof(DIHE));

            k=0;
            while(k<param->nImproper)
            {
                if(fgets(buff3,1024,inout->forfFile)!=NULL)
                {
                    if(buff3[0]=='#')
                        continue;

                    (*impr)[k].a=atoi(strtok(buff3," \n\t"))-1;
                    (*impr)[k].b=atoi(strtok(NULL," \n\t"))-1;
                    (*impr)[k].c=atoi(strtok(NULL," \n\t"))-1;
                    (*impr)[k].d=atoi(strtok(NULL," \n\t"))-1;

                    (*impr)[k].type=atoi(strtok(NULL," \n\t"));
                    (*impr)[k].order=atoi(strtok(NULL," \n\t"));

                    (*impr)[k].k=atof(strtok(NULL," \n\t"))*kcaltoiu;
                    (*impr)[k].phi0=atof(strtok(NULL," \n\t"))*PI/180.;
                    (*impr)[k].mult=atof(strtok(NULL," \n\t"));

                    k++;

                }
                else
                    my_error(PSF_IMPR_SEQ_ERROR,__FILE__,__LINE__,0);
            }

        }
        else if(!strcmp(buff2,"vdw"))
        {

            nAtomCheck=atoi(strtok(NULL," \n\t"));
            if(nAtomCheck!=param->nAtom)
                my_error(CONF_ATNUM_ERROR,__FILE__,__LINE__,2,&(param->nAtom),&nAtomCheck);

            int type;
            double bet;

            k=0;
            while(k<param->nAtom)
            {
                if(fgets(buff3,1024,inout->forfFile)!=NULL)
                {
                    if(buff3[0]=='#')
                        continue;

                    i=atoi(strtok(buff3," \n\t"))-1;
                    if(i!=k)
                        my_error(PSF_FILE_ERROR,__FILE__,__LINE__,0);

                    type=atoi(strtok(NULL," \n\t"));

                    eps[i]=sqrt(atof(strtok(NULL," \n\t"))*kcaltoiu);
                    sig[i]=atof(strtok(NULL," \n\t"));
                    bet=atof(strtok(NULL," \n\t"));

                    eps14[i]=sqrt(atof(strtok(NULL," \n\t"))*kcaltoiu);
                    sig14[i]=atof(strtok(NULL," \n\t"));
                    bet=atof(strtok(NULL," \n\t"));

                    k++;
                }
                else
                    my_error(PSF_BADLINE_ERROR,__FILE__,__LINE__,0);
            }
        }
        else if(!strcmp(buff2,"end"))
        {
            break;
        }
        else
        {
            my_error(UNKNOWN_GENERAL_ERROR,__FILE__,__LINE__,0);
        }
    }

    fclose(inout->forfFile);
}

void write_CONF(IO *inout,PARAM *param,ATOM atom[],double *x,double *y, double *z)
{

    int i;
    double wei=0.;

    inout->rconFile=fopen(inout->rconName,"w");

    if (inout->rconFile==NULL)
    {
        my_error(RESCONF_FILE_ERROR,__FILE__,__LINE__,0);
    }

    fprintf(inout->rconFile,"*Restart Configuration\n");
    fprintf(inout->rconFile,"* Written by MDBas\n");
    fprintf(inout->rconFile,"*\n");
    fprintf(inout->rconFile,"%5d\n",param->nAtom);

    for(i=0; i<param->nAtom; i++)
    {
        fprintf(inout->rconFile,"%5d%5d %-4s %-4s%10.5lf%10.5lf%10.5lf %-4s %-4d%10.5lf\n",
                i+1,atom[i].ires,atom[i].resn,atom[i].label,x[i],y[i],
                z[i],atom[i].segn,atom[i].resi,wei);

    }

    fclose(inout->rconFile);
}

void write_prop(IO *inout,PARAM *param,ENERGY *ener,PBC *box)
{

    double buffer[29]= {0.};
    double temp,press;

    inout->propFile=fopen(inout->propName,"ab");

    temp=2.*ener->kin/((double)param->nDegFree*rboltzui);
    if(box->type>0)
        press=(2.*ener->kin-ener->virtot)/(3.*box->vol*bartoiu);
    else
        press=0.;

    buffer[0]=(double)param->step;
    buffer[1]=(double)param->step*param->timeStep;
    buffer[2]=temp;
    buffer[3]=press;
    buffer[4]=box->vol;
    buffer[5]=ener->tot/kcaltoiu;
    buffer[6]=ener->kin/kcaltoiu;
    buffer[7]=ener->pot/kcaltoiu;
    buffer[8]=ener->elec/kcaltoiu;
    buffer[9]=ener->vdw/kcaltoiu;
    buffer[10]=ener->bond/kcaltoiu;
    buffer[11]=ener->ang/kcaltoiu;
    buffer[12]=ener->ub/kcaltoiu;
    buffer[13]=ener->dihe/kcaltoiu;
    buffer[14]=ener->impr/kcaltoiu;
    buffer[15]=ener->virtot/kcaltoiu;
    buffer[16]=ener->virpot/kcaltoiu;
    buffer[17]=ener->virelec/kcaltoiu;
    buffer[18]=ener->virvdw/kcaltoiu;
    buffer[19]=ener->virbond/kcaltoiu;
    buffer[20]=ener->virub/kcaltoiu;
    buffer[21]=ener->virshake/kcaltoiu;
    buffer[22]=ener->consv/kcaltoiu;

    box_to_lattice(box,&(buffer[23]));

    fwrite(buffer,sizeof(double),29,inout->propFile);

    fclose(inout->propFile);

}

void write_rest(IO *inout,PARAM *param,ENERGY *ener,BATH *bath,ATOM atom[],
                double *x,double *y, double *z,double *vx,double *vy,double *vz,
                double *fx,double *fy, double *fz)
{

    size_t ret;

    inout->restFile=fopen(inout->restName,"wb");

    ret=fwrite(&(param->step),sizeof(int),1,inout->restFile);
    if(ret!=1)
        my_error(UNKNOWN_GENERAL_ERROR,__FILE__,__LINE__,0);

    ret=fwrite(&(param->nAtom),sizeof(int),1,inout->restFile);
    if(ret!=1)
        my_error(UNKNOWN_GENERAL_ERROR,__FILE__,__LINE__,0);

    ret=fwrite(atom,sizeof(ATOM),param->nAtom,inout->restFile);
    if(ret!=(size_t)param->nAtom)
        my_error(UNKNOWN_GENERAL_ERROR,__FILE__,__LINE__,0);

    ret=fwrite(x,sizeof(double),param->nAtom,inout->restFile);
    if(ret!=(size_t)param->nAtom)
        my_error(UNKNOWN_GENERAL_ERROR,__FILE__,__LINE__,0);

    ret=fwrite(y,sizeof(double),param->nAtom,inout->restFile);
    if(ret!=(size_t)param->nAtom)
        my_error(UNKNOWN_GENERAL_ERROR,__FILE__,__LINE__,0);

    ret=fwrite(z,sizeof(double),param->nAtom,inout->restFile);
    if(ret!=(size_t)param->nAtom)
        my_error(UNKNOWN_GENERAL_ERROR,__FILE__,__LINE__,0);

    ret=fwrite(vx,sizeof(double),param->nAtom,inout->restFile);
    if(ret!=(size_t)param->nAtom)
        my_error(UNKNOWN_GENERAL_ERROR,__FILE__,__LINE__,0);

    ret=fwrite(vy,sizeof(double),param->nAtom,inout->restFile);
    if(ret!=(size_t)param->nAtom)
        my_error(UNKNOWN_GENERAL_ERROR,__FILE__,__LINE__,0);

    ret=fwrite(vz,sizeof(double),param->nAtom,inout->restFile);
    if(ret!=(size_t)param->nAtom)
        my_error(UNKNOWN_GENERAL_ERROR,__FILE__,__LINE__,0);

    ret=fwrite(fx,sizeof(double),param->nAtom,inout->restFile);
    if(ret!=(size_t)param->nAtom)
        my_error(UNKNOWN_GENERAL_ERROR,__FILE__,__LINE__,0);

    ret=fwrite(fy,sizeof(double),param->nAtom,inout->restFile);
    if(ret!=(size_t)param->nAtom)
        my_error(UNKNOWN_GENERAL_ERROR,__FILE__,__LINE__,0);

    ret=fwrite(fz,sizeof(double),param->nAtom,inout->restFile);
    if(ret!=(size_t)param->nAtom)
        my_error(UNKNOWN_GENERAL_ERROR,__FILE__,__LINE__,0);

    ret=fwrite(ener,sizeof(ENERGY),1,inout->restFile);
    if(ret!=1)
        my_error(UNKNOWN_GENERAL_ERROR,__FILE__,__LINE__,0);

    ret=fwrite(&(bath->chiT),sizeof(double),1,inout->restFile);
    if(ret!=1)
        my_error(UNKNOWN_GENERAL_ERROR,__FILE__,__LINE__,0);

    ret=fwrite(&(bath->chiP),sizeof(double),1,inout->restFile);
    if(ret!=1)
        my_error(UNKNOWN_GENERAL_ERROR,__FILE__,__LINE__,0);

    fclose(inout->restFile);

}

void write_DCD_header(IO *inout,CTRL *ctrl,PARAM *param,PBC *box,int frozen[])
{
    inout->trajFile=fopen(inout->trajName,"wb");

    char HDR[4]= {'C','O','R','D'};

    int ICNTRL[20]= {0};
    ICNTRL[0]= param->nSteps/ctrl->printTraj;
    ICNTRL[1]= ctrl->printTraj;
    ICNTRL[2]= ctrl->printTraj;
    ICNTRL[3]= param->nSteps;
    ICNTRL[7]= param->nDegFree;
    ICNTRL[8]= param->nFrozen;
    //ICNTRL[9]= timestep in akma but in 32 bits mode
    ICNTRL[10]= (box->type == NOBOX)?0:1;
    ICNTRL[19]= 37;

    int NATOM=param->nAtom;
    int LNFREAT=NATOM-param->nFrozen;

    int NTITLE=3;
    char TITLE1[80]="";
    char TITLE2[80]="";
    char TITLE3[80]="";
    char tmp[4096]="";
    time_t rawtime;
    struct tm *timeinfo;

    time(&rawtime);
    timeinfo = localtime(&rawtime);
    char *user=NULL , *host=NULL , *pwd=NULL;
    user=getenv("USER");
    host=getenv("HOSTNAME");
    pwd=getenv("PWD");

    sprintf(tmp,"* CREATION DATE : %s",asctime(timeinfo));
    strncat(TITLE1,tmp,79);
    sprintf(tmp,"* USER : %s HOSTNAME : %s",(user!=NULL)?user:"UNKNOWN",(host!=NULL)?host:"UNKNOWN");
    strncat(TITLE2,tmp,79);
    sprintf(tmp,"* PWD : %s",(pwd!=NULL)?pwd:"UNKNOWN");
    strncat(TITLE3,tmp,79);

    unsigned int size;

    // first : HDR + ICNTRL
    size = 4*sizeof(char)+20*sizeof(int);
    fwrite(&size,sizeof(unsigned int),1,inout->trajFile);
    fwrite(HDR,sizeof(char),4,inout->trajFile);
    fwrite(ICNTRL,sizeof(int),20,inout->trajFile);
    fwrite(&size,sizeof(unsigned int),1,inout->trajFile);

    // second : NTITLE + TITLE
    size = sizeof(int) + NTITLE*80*sizeof(char);
    fwrite(&size,sizeof(unsigned int),1,inout->trajFile);
    fwrite(&NTITLE,sizeof(int),1,inout->trajFile);
    fwrite(TITLE1,sizeof(char),80,inout->trajFile);
    fwrite(TITLE2,sizeof(char),80,inout->trajFile);
    fwrite(TITLE3,sizeof(char),80,inout->trajFile);
    fwrite(&size,sizeof(unsigned int),1,inout->trajFile);

    // third : NATOM
    size = sizeof(int);
    fwrite(&size,sizeof(unsigned int),1,inout->trajFile);
    fwrite(&NATOM,sizeof(int),1,inout->trajFile);
    fwrite(&size,sizeof(unsigned int),1,inout->trajFile);

    // fourth : if some atoms are frozen, write here the list of free atoms
    if(NATOM != LNFREAT)
    {
        int * freeat=NULL;
        freeat=(int*)my_malloc(LNFREAT*sizeof(int));

        int i=0,j=0;
        for(i=0; i<NATOM; i++)
        {
            if (frozen[i])
                continue;
            else
                freeat[j++]=i;
        }
        size=LNFREAT*sizeof(int);
        fwrite(&size,sizeof(unsigned int),1,inout->trajFile);
        fwrite(freeat,sizeof(int),LNFREAT,inout->trajFile);
        fwrite(&size,sizeof(unsigned int),1,inout->trajFile);

        free(freeat);
    }

    fclose(inout->trajFile);
}

void write_DCD_traj(IO *inout,PARAM *param,PBC *box,double x[],double y[], double z[],int frozen[])
{
    inout->trajFile=fopen(inout->trajName,"ab");

    unsigned int size;

    // if some PBC write a part of the matrix
    if (box->type != NOBOX)
    {
        double box_matrix[6];
        box_to_crystal(box,box_matrix);

        size = 6*sizeof(double);
        fwrite(&size,sizeof(unsigned int),1,inout->trajFile);
        fwrite(box_matrix,sizeof(double),6,inout->trajFile);
        fwrite(&size,sizeof(unsigned int),1,inout->trajFile);
    }

    int siz = (first) ? (param->nAtom) : (param->nAtom - param->nFrozen);

    // alloc of crdinates array + loading of coordinates ; but in float mode !
    float *X = (float*)my_malloc(siz*sizeof(float));
    float *Y = (float*)my_malloc(siz*sizeof(float));
    float *Z = (float*)my_malloc(siz*sizeof(float));

    int i=0,j=0;
    for(i=0; i<param->nAtom; i++)
    {
        if (frozen[i] && !first)
            continue;
        else
        {
            X[j] = (float) x[i];
            Y[j] = (float) y[i];
            Z[j] = (float) z[i];
            j++;
        }
    }

    size = siz*sizeof(float);

    // writing X coordinates
    fwrite(&size,sizeof(unsigned int),1,inout->trajFile);
    fwrite(X,sizeof(float),siz,inout->trajFile);
    fwrite(&size,sizeof(unsigned int),1,inout->trajFile);

    // writing Y coordinates
    fwrite(&size,sizeof(unsigned int),1,inout->trajFile);
    fwrite(Y,sizeof(float),siz,inout->trajFile);
    fwrite(&size,sizeof(unsigned int),1,inout->trajFile);

    // writing Z coordinates
    fwrite(&size,sizeof(unsigned int),1,inout->trajFile);
    fwrite(Z,sizeof(float),siz,inout->trajFile);
    fwrite(&size,sizeof(unsigned int),1,inout->trajFile);

    if(first) first=0;

    fclose(inout->trajFile);
}
