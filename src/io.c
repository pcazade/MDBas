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

#include "global.h"
#include "io.h"
#include "memory.h"
#include "utils.h"

/** Pointer to the output file. **/
extern FILE *outFile;

void read_SIMU(SIMULPARAMS *simulCond,FORCEFIELD *ff,PBC *box)
{
  char buff1[1024]="", *buff2=NULL, *buff3=NULL, *buff4=NULL;
  FILE *simuFile;
  
  simuFile=fopen("SIMU","r");
  
  if (simuFile==NULL)
  {
    error(60);
  }
  
  simulCond->keymd=1;
  
  simulCond->keyrestart=0;
  
  simulCond->keyminim=0;
  simulCond->tolminim=1.e-3;
  simulCond->maxminst=1000;
  simulCond->maxminsiz=0.15;
  
  simulCond->step=0;
  simulCond->firstener=1;
  
  simulCond->mdNature=1;
  simulCond->timeStep=0.001;
  simulCond->nsteps=0;
 
  simulCond->cutoff=12.0;
  simulCond->cuton=10.0;
  simulCond->delr=2.0;

  simulCond->ens=0;
  simulCond->temp=300.0;
  simulCond->press=1.0;
  simulCond->taut=0.1;
  simulCond->taup=0.5;
  simulCond->compres=watercomp;
  
  simulCond->keyrand=0;
  simulCond->seed=12345;
  
  simulCond->elecType=FULL;
  simulCond->vdwType=VFULL;
  simulCond->nb14=0;
  ff->scal14=1.0;
  simulCond->numDeriv=0;
  simulCond->listupdate=20;
  simulCond->linkRatio=1;
  simulCond->nolink=0;

  simulCond->integrator=1;
  simulCond->tolshake=1.e-8;
  simulCond->maxcycle=150;
  simulCond->keyconsth=0;
  simulCond->nconst=0;
  
  simulCond->keyprop=0;
  simulCond->keytraj=0;
  simulCond->keyforf=0;
  simulCond->printo=1000;
  simulCond->printpr=1000;
  simulCond->printtr=1000;
  simulCond->fresconf=1000;
  
  box->type=0;
  box->a1=0.;
  box->b2=0.;
  box->c3=0.;
  
  while(fgets(buff1,1024,simuFile)!=NULL)
  {
    buff2=strtok(buff1," \n\t");
    
    if(buff2==NULL)
      continue;
    
    nocase(buff2);
    
    if(!strcmp(buff2,"mdbas"))
      simulCond->mdNature=0;
    
    else if(!strcmp(buff2,"charmm"))
      simulCond->mdNature=1;
    
    else if(!strcmp(buff2,"namd"))
      simulCond->mdNature=2;
    
    else if(!strcmp(buff2,"nomd"))
      simulCond->keymd=0;
    
    else if(!strcmp(buff2,"restart"))
      simulCond->keyrestart=1;
    
    else if(!strcmp(buff2,"minim"))
    {
      simulCond->keyminim=1;
      buff3=strtok(NULL," \n\t");
      if(buff3==NULL)
	error(63);
      
      nocase(buff3);
      
      if(!strcmp(buff3,"tol"))
      {
	buff4=strtok(NULL," \n\t");
	if(buff4==NULL)
	  error(63);
	
	simulCond->tolminim=atof(buff4);
      }
      else if(!strcmp(buff3,"maxcycle"))
      {
	buff4=strtok(NULL," \n\t");
	if(buff4==NULL)
	  error(63);
	
	simulCond->maxminst=atoi(buff4);
      }
      else if(!strcmp(buff3,"maxsize"))
      {
	buff4=strtok(NULL," \n\t");
	if(buff4==NULL)
	  error(63);
	
	simulCond->maxminsiz=atof(buff4);
      }
      else
	error(62);
    }
    else if(!strcmp(buff2,"timestep"))
    {
      buff3=strtok(NULL," \n\t");
      if(buff3==NULL)
	error(63);
      
      simulCond->timeStep=atof(buff3);
    }
    else if(!strcmp(buff2,"nsteps"))
    {
      buff3=strtok(NULL," \n\t");
      if(buff3==NULL)
	error(63);
      
      simulCond->nsteps=atoi(buff3);
    }
    else if(!strcmp(buff2,"cutoff"))
    {
      buff3=strtok(NULL," \n\t");
      if(buff3==NULL)
	error(63);
      
      simulCond->cutoff=atof(buff3);
    }
    else if(!strcmp(buff2,"cuton"))
    {
      buff3=strtok(NULL," \n\t");
      if(buff3==NULL)
	error(63);
      
      simulCond->cuton=atof(buff3);
    }
    else if(!strcmp(buff2,"delr"))
    {
      buff3=strtok(NULL," \n\t");
      if(buff3==NULL)
	error(63);
      
      simulCond->delr=atof(buff3);
    }
    else if(!strcmp(buff2,"elec"))
    {
      buff3=strtok(NULL," \n\t");
      if(buff3==NULL)
	error(63);
      
      nocase(buff3);
      
      if(!strcmp(buff3,"noelec"))
	simulCond->elecType=NOELEC;
      else if(!strcmp(buff3,"full"))
	simulCond->elecType=FULL;
      else if(!strcmp(buff3,"shift1"))
	simulCond->elecType=SHIFT1;
      else if(!strcmp(buff3,"shift2"))
	simulCond->elecType=SHIFT2;
      else if(!strcmp(buff3,"switch"))
	simulCond->elecType=SWITCH;
      else
	error(62);
    }
    else if(!strcmp(buff2,"vdw"))
    {
      buff3=strtok(NULL," \n\t");
      if(buff3==NULL)
	error(63);
      
      nocase(buff3);
      
      if(!strcmp(buff3,"novdw"))
	simulCond->vdwType=NOVDW;
      else if(!strcmp(buff3,"full"))
	simulCond->vdwType=VFULL;
      else if(!strcmp(buff3,"switch"))
	simulCond->vdwType=VSWITCH;
      else
	error(62);
    }
    else if(!strcmp(buff2,"nb14"))
    {
      simulCond->nb14=1;
      buff3=strtok(NULL," \n\t");
      if(buff3==NULL)
	error(63);
      
      ff->scal14=atof(buff3);
    }
    else if(!strcmp(buff2,"numforce"))
    {
      simulCond->numDeriv=1;
    }
    else if(!strcmp(buff2,"list"))
    {
      buff3=strtok(NULL," \n\t");
      if(buff3==NULL)
	error(63);
      
      simulCond->listupdate=atoi(buff3);
    }
    else if(!strcmp(buff2,"link"))
    {
      buff3=strtok(NULL," \n\t");
      if(buff3==NULL)
	error(63);
      
      simulCond->linkRatio=atoi(buff3);
    }
    else if(!strcmp(buff2,"nolink"))
    {
      simulCond->nolink=1;
    }
    else if(!strcmp(buff2,"integrator"))
    {
      buff3=strtok(NULL," \n\t");
      if(buff3==NULL)
	error(63);
      
      nocase(buff3);
      
      if(!strcmp(buff3,"leapfrog"))
	simulCond->integrator=0;
      else if(!strcmp(buff3,"velocity"))
	simulCond->integrator=1;
      else
	error(62);
    }
    else if(!strcmp(buff2,"ensemble"))
    {
      buff3=strtok(NULL," \n\t");
      if(buff3==NULL)
	error(63);
      
      nocase(buff3);
      
      if(!strcmp(buff3,"nve"))
	simulCond->ens=0;
      else if(!strcmp(buff3,"nvtb"))
	simulCond->ens=1;
      else if(!strcmp(buff3,"nptb"))
	simulCond->ens=2;
      else if(!strcmp(buff3,"nvth"))
	simulCond->ens=3;
      else if(!strcmp(buff3,"npth"))
	simulCond->ens=4;
      else
	error(62);
      
      if(simulCond->ens>0)
      {
	buff3=strtok(NULL," \n\t");
	if(buff3==NULL)
	error(63);
	
	simulCond->taut=atof(buff3);
      }
      
      if(simulCond->ens==2||simulCond->ens==4)
      {
	buff3=strtok(NULL," \n\t");
	if(buff3==NULL)
	error(63);
	
	simulCond->taup=atof(buff3);
      }
    }
    else if(!strcmp(buff2,"temperature"))
    {
      buff3=strtok(NULL," \n\t");
      if(buff3==NULL)
	error(63);
      
      simulCond->temp=atof(buff3);
    }
    else if(!strcmp(buff2,"pressure"))
    {
      buff3=strtok(NULL," \n\t");
      if(buff3==NULL)
	error(63);
      
      simulCond->press=atof(buff3)*bartoiu;
    }
    else if(!strcmp(buff2,"compressibility"))
    {
      buff3=strtok(NULL," \n\t");
      if(buff3==NULL)
	error(63);
      
      simulCond->compres=atof(buff3)/bartoiu;
    }
    else if(!strcmp(buff2,"consth"))
    {
      simulCond->keyconsth=1;
    }
    else if(!strcmp(buff2,"shake"))
    {
      buff3=strtok(NULL," \n\t");
      if(buff3==NULL)
	error(63);
      
      nocase(buff3);
      
      if(!strcmp(buff3,"tol"))
      {
	buff4=strtok(NULL," \n\t");
	if(buff4==NULL)
	  error(63);
	
	simulCond->tolshake=atof(buff4);
      }
      else if(!strcmp(buff3,"maxcycle"))
      {
	buff4=strtok(NULL," \n\t");
	if(buff4==NULL)
	  error(63);
	
	simulCond->maxcycle=atoi(buff4);
      }
      else
	error(62);
    }
    else if(!strcmp(buff2,"seed"))
    {
      buff3=strtok(NULL," \n\t");
	if(buff3==NULL)
	error(63);
	
	simulCond->keyrand=1;
	simulCond->seed=atoi(buff3);
    }
    else if(!strcmp(buff2,"print"))
    {
      buff3=strtok(NULL," \n\t");
      if(buff3==NULL)
	error(63);
      
      simulCond->printo=atoi(buff3);
    }
    else if(!strcmp(buff2,"prop"))
    {
      buff3=strtok(NULL," \n\t");
      if(buff3==NULL)
	error(63);
      
      simulCond->keyprop=1;
      simulCond->printpr=atoi(buff3);
    }
    else if(!strcmp(buff2,"traj"))
    {
      buff3=strtok(NULL," \n\t");
      if(buff3==NULL)
	error(63);
      
      simulCond->keytraj=1;
      simulCond->printtr=atoi(buff3);
    }
    else if(!strcmp(buff2,"resconf"))
    {
      buff3=strtok(NULL," \n\t");
      if(buff3==NULL)
	error(63);
      
      simulCond->fresconf=atoi(buff3);
    }
    else if(!strcmp(buff2,"write"))
    {
      buff3=strtok(NULL," \n\t");
      if(buff3==NULL)
	error(63);
      
      nocase(buff3);
      
      if(!strcmp(buff3,"field"))
	simulCond->keyforf=1;
      else
	error(62);
    }
    else if(!strcmp(buff2,"pbc"))
    {
      buff3=strtok(NULL," \n\t");
      if(buff3==NULL)
	error(63);
      
      box->type=atoi(buff3);
      
      if(fgets(buff1,1024,simuFile)!=NULL)
      {
	
	buff2=strtok(buff1," \n\t");
	if(buff2==NULL)
	  error(63);
	else if(isdigit(buff2[0])==0)
	  error(62);
	
	box->a1=atof(buff2);
	
	buff2=strtok(NULL," \n\t");
	if(buff2==NULL)
	  error(63);
	else if(isdigit(buff2[0])==0)
	  error(62);
	
	box->a2=atof(buff2);
	
	buff2=strtok(NULL," \n\t");
	if(buff2==NULL)
	  error(63);
	else if(isdigit(buff2[0])==0)
	  error(62);
	
	box->a3=atof(buff2);
	
      }
      else
	error(63);
	
      if(fgets(buff1,1024,simuFile)!=NULL)
      {
	
	buff2=strtok(buff1," \n\t");
	if(buff2==NULL)
	  error(63);
	else if(isdigit(buff2[0])==0)
	  error(62);
	
	box->b1=atof(buff2);
	
	buff2=strtok(NULL," \n\t");
	if(buff2==NULL)
	  error(63);
	else if(isdigit(buff2[0])==0)
	  error(62);
	
	box->b2=atof(buff2);
	
	buff2=strtok(NULL," \n\t");
	if(buff2==NULL)
	  error(63);
	else if(isdigit(buff2[0])==0)
	  error(62);
	
	box->b3=atof(buff2);
	
      }
      else
	error(63);
      
      if(fgets(buff1,1024,simuFile)!=NULL)
      {
	
	buff2=strtok(buff1," \n\t");
	if(buff2==NULL)
	  error(63);
	else if(isdigit(buff2[0])==0)
	  error(62);
	
	box->c1=atof(buff2);
	
	buff2=strtok(NULL," \n\t");
	if(buff2==NULL)
	  error(63);
	else if(isdigit(buff2[0])==0)
	  error(62);
	
	box->c2=atof(buff2);
	
	buff2=strtok(NULL," \n\t");
	if(buff2==NULL)
	  error(63);
	else if(isdigit(buff2[0])==0)
	  error(62);
	
	box->c3=atof(buff2);
	
      }
      else
	error(63);
    }
    else if(!strcmp(buff2,"end"))
      break;
    else
      error(61);
  }
  
}

void read_PSF(INPUTS *inp,ATOM **atom,FORCEFIELD *ff,SIMULPARAMS *simulCond,CONSTRAINT **constList)
{
  FILE *psfFile=NULL;
  char buff1[1024]="", *buff2=NULL, *buff3=NULL, *buff4=NULL, *buff5=NULL;
  int i,j,k,kk,kt,nalloc=100,nincr=100;
  int ia,ib;
  
  psfFile=fopen("PSF","r");
  
  if(psfFile==NULL)
  {
    error(20);
  }
  
  while(fgets(buff1,1024,psfFile)!=NULL)
  {
    if(strstr(buff1,"NATOM")!=NULL)
    {
      buff2=strtok(buff1," \n\t");
      simulCond->natom=atoi(buff2);
      break;
    }
  }
  
  inp->typesNum=(int*)malloc(nalloc*sizeof(*(inp->typesNum)));

  *atom=(ATOM*)malloc(simulCond->natom*sizeof(**atom));
  
  k=-1;
  for(i=0;i<simulCond->natom;i++)
  {
    if(fgets(buff1,1024,psfFile)!=NULL)
    {
      buff2=strtok(buff1," \n\t");
      buff3=strtok(NULL," \n\t");
      buff4=strtok(NULL," \n\t");
      buff5=strtok(NULL," \n\t");

      strcpy((*atom)[i].segn,buff3);
      (*atom)[i].resi=atoi(buff4);
      strcpy((*atom)[i].resn,buff5);
      
      buff2=strtok(NULL," \n\t");
      buff3=strtok(NULL," \n\t");
      buff4=strtok(NULL," \n\t");
      buff5=strtok(NULL," \n\t");
      
      strcpy((*atom)[i].label,buff2);
      
      kt=atoi(buff3)-1;
      k++;
      
      if(k>=nalloc)
      {
	nalloc+=nincr;
	inp->typesNum=(int*)realloc(inp->typesNum,nalloc*sizeof(*(inp->typesNum)));
      }
      
      kk=k;
      inp->typesNum[k]=kt;
      for(j=0;j<k;j++)
      {
	if(kt==inp->typesNum[j])
	{
	  k--;
	  kk=j;
	  break;
	}
      }
//       k=inp->typesNum[atoi(buff3)-1];
//       if(k==-1)
//  	error(21);

      (*atom)[i].type=kk;
      (*atom)[i].q=atof(buff4);
      (*atom)[i].m=atof(buff5);
    }
    else
    {
      error(22);
    }
  }
  
  inp->nTypes=k+1;
  inp->typesNum=(int*)realloc(inp->typesNum,(inp->nTypes*sizeof*(inp->typesNum)));
  
  while(fgets(buff1,1024,psfFile)!=NULL)
  {
    if(strstr(buff1,"NBOND")!=NULL)
    {
      buff2=strtok(buff1," \n\t");
      ff->nBond=atoi(buff2);
      break;
    }
  }
  
  simulCond->iBond=(int**)malloc(ff->nBond*sizeof(*(simulCond->iBond)));
  for(i=0;i<ff->nBond;i++)
    simulCond->iBond[i]=(int*)malloc(2*sizeof(**(simulCond->iBond)));
  
  if(simulCond->keyconsth)
  {
    for(i=0;i<simulCond->natom;i++)
      (*atom)[i].inconst=0;
    
    *constList=(CONSTRAINT*)malloc(ff->nBond*sizeof(**constList));
    simulCond->nconst=0;
  }
  
  i=0;
  while(fgets(buff1,1024,psfFile)!=NULL)
  {
    buff2=strtok(buff1," \n\t");
    buff3=strtok(NULL," \n\t");
    
    while(buff2!=NULL && buff3!=NULL)
    {
      
      if(simulCond->keyconsth)
      {
	ia=atoi(buff2)-1;
	ib=atoi(buff3)-1;
	if((*atom)[ia].label[0]=='H'||(*atom)[ib].label[0]=='H')
	{
	  (*constList)[simulCond->nconst].a=ia;
	  (*constList)[simulCond->nconst].b=ib;
	  
	  (*atom)[ia].inconst++;
	  (*atom)[ib].inconst++;
	  
	  buff2=strtok(NULL," \n\t");
	  buff3=strtok(NULL," \n\t");
	  
	  simulCond->nconst++;
	  ff->nBond--;
	}
	else
	{
	  simulCond->iBond[i][0]=ia;
	  simulCond->iBond[i][1]=ib;
	
	  buff2=strtok(NULL," \n\t");
	  buff3=strtok(NULL," \n\t");
	  i++;
	}
      }
      else
      {
	simulCond->iBond[i][0]=atoi(buff2)-1;
	simulCond->iBond[i][1]=atoi(buff3)-1;
	
	buff2=strtok(NULL," \n\t");
	buff3=strtok(NULL," \n\t");
	i++;
      }
      
    }
    if(buff2!=NULL && buff3==NULL)
      error(23);
    if(i==ff->nBond)
      break;
  }
  
  if(simulCond->keyconsth)
  {
    *constList=(CONSTRAINT*)realloc(*constList,simulCond->nconst*sizeof(**constList));
  }
  
  while(fgets(buff1,1024,psfFile)!=NULL)
  {
    if(strstr(buff1,"NTHETA")!=NULL)
    {
      buff2=strtok(buff1," \n\t");
      ff->nAngle=atoi(buff2);
      break;
    }
  }
  
  simulCond->iAngle=(int**)malloc(ff->nAngle*sizeof(*(simulCond->iAngle)));
  for(i=0;i<ff->nAngle;i++)
    simulCond->iAngle[i]=(int*)malloc(3*sizeof(**(simulCond->iAngle)));
  
  i=0;
  while(fgets(buff1,1024,psfFile)!=NULL)
  {
    buff2=strtok(buff1," \n\t");
    buff3=strtok(NULL," \n\t");
    buff4=strtok(NULL," \n\t");
    
    while(buff2!=NULL && buff3!=NULL && buff4!=NULL)
    {
      
      simulCond->iAngle[i][0]=atoi(buff2)-1;
      simulCond->iAngle[i][1]=atoi(buff3)-1;
      simulCond->iAngle[i][2]=atoi(buff4)-1;
      
      buff2=strtok(NULL," \n\t");
      buff3=strtok(NULL," \n\t");
      buff4=strtok(NULL," \n\t");
      i++;
    }
    if((buff2!=NULL && buff3==NULL)||(buff2!=NULL && buff4==NULL))
      error(24);
    if(i==ff->nAngle)
      break;
  }
  
  while(fgets(buff1,1024,psfFile)!=NULL)
  {
    if(strstr(buff1,"NPHI")!=NULL)
    {
      buff2=strtok(buff1," \n\t");
      ff->nDihedral=atoi(buff2);
      break;
    }
  }
  
  simulCond->iDihedral=(int**)malloc(ff->nAngle*sizeof(*(simulCond->iDihedral)));
  for(i=0;i<ff->nDihedral;i++)
    simulCond->iDihedral[i]=(int*)malloc(4*sizeof(**(simulCond->iDihedral)));
  
  i=0;
  while(fgets(buff1,1024,psfFile)!=NULL)
  {
    buff2=strtok(buff1," \n\t");
    buff3=strtok(NULL," \n\t");
    buff4=strtok(NULL," \n\t");
    buff5=strtok(NULL," \n\t");
    
    while(buff2!=NULL && buff3!=NULL && buff4!=NULL)
    {
      
      simulCond->iDihedral[i][0]=atoi(buff2)-1;
      simulCond->iDihedral[i][1]=atoi(buff3)-1;
      simulCond->iDihedral[i][2]=atoi(buff4)-1;
      simulCond->iDihedral[i][3]=atoi(buff5)-1;
      
      buff2=strtok(NULL," \n\t");
      buff3=strtok(NULL," \n\t");
      buff4=strtok(NULL," \n\t");
      buff5=strtok(NULL," \n\t");
      i++;
    }
    if((buff2!=NULL && buff3==NULL)||(buff2!=NULL && buff4==NULL)||(buff2!=NULL && buff5==NULL))
      error(25);
    if(i==ff->nDihedral)
      break;
  }
  
  while(fgets(buff1,1024,psfFile)!=NULL)
  {
    if(strstr(buff1,"NIMPHI")!=NULL)
    {
      buff2=strtok(buff1," \n\t");
      ff->nImproper=atoi(buff2);
      break;
    }
  }
  
  simulCond->iImproper=(int**)malloc(ff->nAngle*sizeof(*(simulCond->iImproper)));
  for(i=0;i<ff->nImproper;i++)
    simulCond->iImproper[i]=(int*)malloc(4*sizeof(**(simulCond->iImproper)));
  
  i=0;
  while(fgets(buff1,1024,psfFile)!=NULL)
  {
    buff2=strtok(buff1," \n\t");
    buff3=strtok(NULL," \n\t");
    buff4=strtok(NULL," \n\t");
    buff5=strtok(NULL," \n\t");
    
    while(buff2!=NULL && buff3!=NULL && buff4!=NULL)
    {
      
      simulCond->iImproper[i][0]=atoi(buff2)-1;
      simulCond->iImproper[i][1]=atoi(buff3)-1;
      simulCond->iImproper[i][2]=atoi(buff4)-1;
      simulCond->iImproper[i][3]=atoi(buff5)-1;
      
      buff2=strtok(NULL," \n\t");
      buff3=strtok(NULL," \n\t");
      buff4=strtok(NULL," \n\t");
      buff5=strtok(NULL," \n\t");
      i++;
    }
    if((buff2!=NULL && buff3==NULL)||(buff2!=NULL && buff4==NULL)||(buff2!=NULL && buff5==NULL))
      error(26);
    if(i==ff->nImproper)
      break;
  }
  
  fclose(psfFile);
}

void read_TOP(INPUTS *inp)
{
  
  FILE *topFile=NULL;
  char buff1[1024]="", *buff2=NULL, *buff3=NULL, *buff4=NULL ;
  int i,k;
  
  topFile=fopen("TOP","r");

  if (topFile==NULL)
  {
    error(10);
  }
  
  inp->types=(char**)malloc(inp->nTypes*sizeof(*(inp->types)));
  for(i=0;i<inp->nTypes;i++)
      inp->types[i]=(char*)malloc(5*sizeof(**(inp->types)));
  
  while(fgets(buff1,1024,topFile)!=NULL)
  {
    
    buff2=strtok(buff1," \n\t");
    
    if (buff2 != NULL)
    {	
      if (!strcmp(buff2,"MASS"))
      {
	buff3=strtok(NULL," \n\t");
	buff4=strtok(NULL," \n\t");
	k=atoi(buff3)-1;
	for(i=0;i<inp->nTypes;i++)
	{
	  if(k==inp->typesNum[i])
	  {
	    strcpy(inp->types[i],buff4);
	    break;
	  }
	}
      }
    }
  }
  
  fclose(topFile);
    
}

void read_PAR(INPUTS *inp)
{
  
  char buff1[1024]="", *buff2=NULL, *buff3=NULL, *buff4=NULL, *buff5=NULL, *buff6=NULL;

    FILE *parFile=NULL;
    int i,j,l,k=0;
    
    parFile=fopen("PAR","r");

    if(parFile==NULL)
    {
      error(30);
    }
    
    while(fgets(buff1,1024,parFile)!=NULL)
    {
      buff2=strtok(buff1," \n\t");
      
      if(buff2==NULL)
	continue;
      
      if(!strcmp(buff2,"BONDS"))
      {
	k=1;
	inp->nBondTypes=0;
	continue;
      }
      else if(!strcmp(buff2,"ANGLES"))
      {
	k=2;
	inp->nAngTypes=0;
	inp->nUbTypes=0;
	continue;
      }
      else if(!strcmp(buff2,"DIHEDRALS"))
      {
	k=3;
	inp->nDiheTypes=0;
	continue;
      }
      else if(!strcmp(buff2,"IMPROPER"))
      {
	k=4;
	inp->nImprTypes=0;
	continue;
      }
      else if(!strcmp(buff2,"NONBONDED"))
      {
	k=5;
	inp->nNonBonded=0;
	continue;
      }
      else if(!strcmp(buff2,"NBFIX"))
      {
	k=6;
	continue;
      }
      else if(!strcmp(buff2,"CMAP"))
      {
	k=7;
	continue;
      }
      else if(!strcmp(buff2,"HBOND"))
      {
	k=8;
	continue;
      }
      else if(!strcmp(buff2,"END"))
      {
	break;
      }
      
      
      if(k==1)      
      {
	if((buff2!=NULL)&&(buff2[0]!='!'))
	    inp->nBondTypes++;
      }
      else if(k==2)      
      {
	if((buff2!=NULL)&&(buff2[0]!='!'))
	{
	    inp->nAngTypes++;
	    
	    buff2=strtok(NULL," \n\t");
	    buff2=strtok(NULL," \n\t");
	    buff2=strtok(NULL," \n\t");
	    buff2=strtok(NULL," \n\t");
	    buff2=strtok(NULL," \n\t");
	    if((buff2!=NULL)&&(buff2[0]!='!'))
	      inp->nUbTypes++;
	}
      }
      else if(k==3)      
      {
	if((buff2!=NULL)&&(buff2[0]!='!'))
	    inp->nDiheTypes++;
      }
      else if(k==4)      
      {
	if((buff2!=NULL)&&(buff2[0]!='!'))
	    inp->nImprTypes++;
      }
      else if(k==5)      
      {
	if((buff2!=NULL)&&(buff2[0]!='!'))
	    inp->nNonBonded++;
      }
    }
    
//       Bond types arrays allocation
    
    inp->bondTypes=(int**)malloc(inp->nBondTypes*sizeof(*(inp->bondTypes)));
    inp->bondTypesParm=(double**)malloc(inp->nBondTypes*sizeof(*(inp->bondTypesParm)));
    
    for(i=0;i<inp->nBondTypes;i++)
    {
      inp->bondTypes[i]=(int*)malloc(2*sizeof(**(inp->bondTypes)));
      inp->bondTypesParm[i]=(double*)malloc(2*sizeof(**(inp->bondTypesParm)));
    }
    
//       Angle types arrays allocation
    
    inp->angTypes=(int**)malloc(inp->nAngTypes*sizeof(*(inp->angTypes)));
    inp->angTypesParm=(double**)malloc(inp->nAngTypes*sizeof(*(inp->angTypesParm)));
    
    for(i=0;i<inp->nAngTypes;i++)
    {
      inp->angTypes[i]=(int*)malloc(3*sizeof(**(inp->angTypes)));
      inp->angTypesParm[i]=(double*)malloc(2*sizeof(**(inp->angTypesParm)));
    }
    
//       Uray-Bradley types arrays allocation
    
    inp->ubTypes=(int**)malloc(inp->nUbTypes*sizeof(*(inp->ubTypes)));
    inp->ubTypesParm=(double**)malloc(inp->nUbTypes*sizeof(*(inp->ubTypesParm)));
        
    for(i=0;i<inp->nUbTypes;i++)
    {
      inp->ubTypes[i]=(int*)malloc(3*sizeof(**(inp->ubTypes)));
      inp->ubTypesParm[i]=(double*)malloc(2*sizeof(**(inp->ubTypesParm)));
    }
    
//       Dihedral types arrays allocation
    
    inp->diheTypes=(int**)malloc(inp->nDiheTypes*sizeof(*(inp->diheTypes)));
    inp->diheTypesParm=(double**)malloc(inp->nDiheTypes*sizeof(*(inp->diheTypesParm)));
    inp->nDiheTypesParm=(int*)malloc(inp->nDiheTypes*sizeof(*(inp->nDiheTypesParm)));
        
    for(i=0;i<inp->nDiheTypes;i++)
    {
      inp->diheTypes[i]=(int*)malloc(4*sizeof(**(inp->diheTypes)));
      inp->diheTypesParm[i]=(double*)malloc(3*sizeof(**(inp->diheTypesParm)));
    }
    
    for(i=0;i<inp->nDiheTypes;i++)
      inp->nDiheTypesParm[i]=1;
    
//       Improper types arrays allocation
    
    inp->imprTypes=(int**)malloc(inp->nImprTypes*sizeof(*(inp->imprTypes)));
    inp->imprTypesParm=(double**)malloc(inp->nImprTypes*sizeof(*(inp->imprTypesParm)));
    
    for(i=0;i<inp->nImprTypes;i++)
    {
      inp->imprTypes[i]=(int*)malloc(4*sizeof(**(inp->imprTypes)));
      inp->imprTypesParm[i]=(double*)malloc(3*sizeof(**(inp->imprTypesParm)));
    }
    
//       Non-bonding types arrays allocation
    
    inp->nonBondedTypesParm=(double**)malloc(inp->nTypes*sizeof(*(inp->nonBondedTypesParm)));
    
    for(i=0;i<inp->nTypes;i++)
      inp->nonBondedTypesParm[i]=(double*)malloc(6*sizeof(**(inp->nonBondedTypesParm)));
    
    for(i=0;i<inp->nTypes;i++)
    {
      for(j=0;j<6;j++)
	inp->nonBondedTypesParm[i][j]=-100.;
    }
      
    rewind(parFile);
    
    int ia,ib,ic,id,ii;
    k=0;
    
    while(fgets(buff1,1024,parFile)!=NULL)
    {
      
      buff2=strtok(buff1," \n\t");
      
      if(buff2==NULL)
	continue;
      
      if(!strcmp(buff2,"BONDS"))
      {
	k=1;
	i=0;
	inp->nBondTypes=0;
	continue;
      }
      else if(!strcmp(buff2,"ANGLES"))
      {
	k=2;
	i=0;
	j=0;
	inp->nAngTypes=0;
	inp->nUbTypes=0;
	continue;
      }
      else if(!strcmp(buff2,"DIHEDRALS"))
      {
	k=3;
	i=0;
	inp->nDiheTypes=0;
	continue;
      }
      else if(!strcmp(buff2,"IMPROPER"))
      {
	k=4;
	i=0;
	inp->nImprTypes=0;
	continue;
      }
      else if(!strcmp(buff2,"NONBONDED"))
      {
	k=5;
	i=0;
	inp->nNonBonded=0;
	continue;
      }
      else if(!strcmp(buff2,"NBFIX"))
      {
	k=6;
	i=0;
	continue;
      }
      else if(!strcmp(buff2,"CMAP"))
      {
	k=7;
	i=0;
	continue;
      }
      else if(!strcmp(buff2,"HBOND"))
      {
	k=8;
	i=0;
	continue;
      }
      else if(!strcmp(buff2,"END"))
      {
	break;
      }
      
      
      if(k==1)      
      {
	if((buff2!=NULL)&&(buff2[0]!='!'))
	{
	  buff3=strtok(NULL," \n\t");
	  buff4=strtok(NULL," \n\t");
	  buff5=strtok(NULL," \n\t");
	  
	  ia=-1;
	  ib=-1;
	  for(l=0;l<inp->nTypes;l++)
	  {
	    if(!strcmp(buff2,inp->types[l]))
	    {
	      ia=l;
	      break;
	    }
	  }
	  
	  if(ia==-1)
	    continue;
	  
	  for(l=0;l<inp->nTypes;l++)
	  {
	    if(!strcmp(buff3,inp->types[l]))
	    {
	      ib=l;
	      break;
	    }
	  }
	  
	  if(ib==-1)
	    continue;
	  
	  if(ib<ia)
	  {
	    inp->bondTypes[i][0]=ib;
	    inp->bondTypes[i][1]=ia;
	  }
	  else
	  {
	    inp->bondTypes[i][0]=ia;
	    inp->bondTypes[i][1]=ib;
	  }
	  
	  inp->bondTypesParm[i][0]=atof(buff4);
	  inp->bondTypesParm[i][1]=atof(buff5);
	  
	  i++;
	  inp->nBondTypes=i;
	}
      }
      else if(k==2)      
      {
	if((buff2!=NULL)&&(buff2[0]!='!'))
	{
	  buff3=strtok(NULL," \n\t");
	  buff4=strtok(NULL," \n\t");
	  buff5=strtok(NULL," \n\t");
	  buff6=strtok(NULL," \n\t");
	  
	  ia=-1;
	  ib=-1;
	  ic=-1;
	  
	  for(l=0;l<inp->nTypes;l++)
	  {
	    if(!strcmp(buff2,inp->types[l]))
	    {
	      ia=l;
	      break;
	    }
	  }
	  
	  if(ia==-1)
	    continue;
	  
	  for(l=0;l<inp->nTypes;l++)
	  {
	    if(!strcmp(buff3,inp->types[l]))
	    {
	      ib=l;
	      break;
	    }
	  }
	  
	  if(ib==-1)
	    continue;
	  
	  for(l=0;l<inp->nTypes;l++)
	  {
	    if(!strcmp(buff4,inp->types[l]))
	    {
	      ic=l;
	      break;
	    }
	  }
	  
	  if(ic==-1)
	    continue;
	  
	  if(ic<ia)
	  {
	    inp->angTypes[i][0]=ic;
	    inp->angTypes[i][1]=ib;
	    inp->angTypes[i][2]=ia;
	  }
	  else
	  {
	    inp->angTypes[i][0]=ia;
	    inp->angTypes[i][1]=ib;
	    inp->angTypes[i][2]=ic;
	  }
	  
	  inp->angTypesParm[i][0]=atof(buff5);
	  inp->angTypesParm[i][1]=atof(buff6);
	  
	  buff2=strtok(NULL," \n\t");
	  if((buff2!=NULL)&&(buff2[0]!='!'))
	  {
	    buff3=strtok(NULL," \n\t");
	    
	    inp->ubTypesParm[j][0]=atof(buff2);
	    inp->ubTypesParm[j][1]=atof(buff3);
	    
	    if(ic<ia)
	    {
	      inp->ubTypes[j][0]=ic;
	      inp->ubTypes[j][1]=ib;
	      inp->ubTypes[j][2]=ia;
	    }
	    else
	    {
	      inp->ubTypes[j][0]=ia;
	      inp->ubTypes[j][1]=ib;
	      inp->ubTypes[j][2]=ic;
	    }
	    j++;
	    inp->nUbTypes=j;
	  }
	  i++;
	  inp->nAngTypes=i;
	}
      }
      else if(k==3)      
      {
	int itype,index,i0,i1,i2,i3;
	
	if((buff2!=NULL)&&(buff2[0]!='!'))
	{
	  buff3=strtok(NULL," \n\t");
	  buff4=strtok(NULL," \n\t");
	  buff5=strtok(NULL," \n\t");
	  
	  ia=-1;
	  ib=-1;
	  ic=-1;
	  id=-1;
	  
	  if(!strcmp(buff2,"X"))
	  {
	    ia=inp->nTypes;
	  }
	  else
	  {
	    for(l=0;l<inp->nTypes;l++)
	    {
	      if(!strcmp(buff2,inp->types[l]))
	      {
		ia=l;
		break;
	      }
	    }
	  }
	  
	  if(ia==-1)
	    continue;
	  
	  for(l=0;l<inp->nTypes;l++)
	  {
	    if(!strcmp(buff3,inp->types[l]))
	    {
	      ib=l;
	      break;
	    }
	  }
	  
	  if(ib==-1)
	    continue;
	  
	  for(l=0;l<inp->nTypes;l++)
	  {
	    if(!strcmp(buff4,inp->types[l]))
	    {
	      ic=l;
	      break;
	    }
	  }
	  
	  if(ic==-1)
	    continue;
	  
	  if(!strcmp(buff5,"X"))
	  {
	    id=inp->nTypes;
	  }
	  else 
	  {
	    for(l=0;l<inp->nTypes;l++)
	    {
	      if(!strcmp(buff5,inp->types[l]))
	      {
		id=l;
		break;
	      }
	    }
	  }
	  
	  if(id==-1)
	    continue;
	  
	  if(id<ia)
	  {
	    i0=id;
	    i1=ic;
	    i2=ib;
	    i3=ia;
	  }
	  else if(ia<id)
	  {
	    i0=ia;
	    i1=ib;
	    i2=ic;
	    i3=id;
	  }
	  else if(ic<ib)
	  {
	    i0=id;
	    i1=ic;
	    i2=ib;
	    i3=ia;
	  }
	  else
	  {
	    i0=ia;
	    i1=ib;
	    i2=ic;
	    i3=id;
	  }
	  
	  itype=i;
	  for(l=0;l<inp->nDiheTypes;l++)
	  {
	    if( (i0==inp->diheTypes[l][0]) && (i1==inp->diheTypes[l][1]) && (i2==inp->diheTypes[l][2]) && (i3==inp->diheTypes[l][3]) )
	    {
	      itype=l;
	      break;
	    }
	  }
	  
	  buff3=strtok(NULL," \n\t");
	  buff4=strtok(NULL," \n\t");
	  buff5=strtok(NULL," \n\t");
	  
	  if(itype==i)
	  {
	    inp->diheTypes[l][0]=i0;
	    inp->diheTypes[l][1]=i1;
	    inp->diheTypes[l][2]=i2;
	    inp->diheTypes[l][3]=i3;
	    
	    inp->diheTypesParm[i][0]=atof(buff3);
	    inp->diheTypesParm[i][1]=atof(buff4);
	    inp->diheTypesParm[i][2]=atof(buff5);
	    
	    i++;
	    inp->nDiheTypes=i;
	    
	  }
	  else
	  {
	    inp->nDiheTypesParm[itype]++;
	    inp->diheTypesParm[itype]=(double*)realloc(inp->diheTypesParm[itype],
	      inp->nDiheTypesParm[itype]*3*sizeof(*(inp->diheTypesParm)));
	    
	    index=0+3*(inp->nDiheTypesParm[itype]-1);
	    
	    inp->diheTypesParm[itype][index]=atof(buff3);
	    inp->diheTypesParm[itype][index+1]=atof(buff4);
	    inp->diheTypesParm[itype][index+2]=atof(buff5);
	    
	  }
	}
      }
      else if(k==4)      
      {
	if((buff2!=NULL)&&(buff2[0]!='!'))
	{
	  buff3=strtok(NULL," \n\t");
	  buff4=strtok(NULL," \n\t");
	  buff5=strtok(NULL," \n\t");
	  
	  ia=-1;
	  ib=-1;
	  ic=-1;
	  id=-1;
	  
	  if(!strcmp(buff2,"X"))
	  {
	    ia=inp->nTypes;
	  }
	  else
	  {
	    for(l=0;l<inp->nTypes;l++)
	    {
	      if(!strcmp(buff2,inp->types[l]))
	      {
		ia=l;
		break;
	      }
	    }
	  }
	  
	  if(ia==-1)
	    continue;
	  
	  if(!strcmp(buff3,"X"))
	  {
	    ib=inp->nTypes;
	  }
	  else
	  {
	    for(l=0;l<inp->nTypes;l++)
	    {
	      if(!strcmp(buff3,inp->types[l]))
	      {
		ib=l;
		break;
	      }
	    }
	  }
	  
	  if(ib==-1)
	    continue;
	  
	  if(!strcmp(buff4,"X"))
	  {
	    ic=inp->nTypes;
	  }
	  else
	  {
	    for(l=0;l<inp->nTypes;l++)
	    {
	      if(!strcmp(buff4,inp->types[l]))
	      {
		ic=l;
		break;
	      }
	    }
	  }
	  
	  if(ic==-1)
	    continue;
	  
	  if(!strcmp(buff5,"X"))
	  {
	    id=inp->nTypes;
	  }
	  else 
	  {
	    for(l=0;l<inp->nTypes;l++)
	    {
	      if(!strcmp(buff5,inp->types[l]))
	      {
		id=l;
		break;
	      }
	    }
	  }
	  
	  if(id==-1)
	    continue;
	  
	  if(id<ia)
	  {
	    inp->imprTypes[i][0]=id;
	    inp->imprTypes[i][1]=ic;
	    inp->imprTypes[i][2]=ib;
	    inp->imprTypes[i][3]=ia;
	  }
	  else if(ia<id)
	  {
	    inp->imprTypes[i][0]=ia;
	    inp->imprTypes[i][1]=ib;
	    inp->imprTypes[i][2]=ic;
	    inp->imprTypes[i][3]=id;
	  }
	  else if(ic<ib)
	  {
	    inp->imprTypes[i][0]=id;
	    inp->imprTypes[i][1]=ic;
	    inp->imprTypes[i][2]=ib;
	    inp->imprTypes[i][3]=ia;
	  }
	  else
	  {
	    inp->imprTypes[i][0]=ia;
	    inp->imprTypes[i][1]=ib;
	    inp->imprTypes[i][2]=ic;
	    inp->imprTypes[i][3]=id;
	  }
	  
	  buff3=strtok(NULL," \n\t");
	  buff4=strtok(NULL," \n\t");
	  buff5=strtok(NULL," \n\t");
	
	  inp->imprTypesParm[i][0]=atof(buff3);
	  inp->imprTypesParm[i][1]=atof(buff4);
	  inp->imprTypesParm[i][2]=atof(buff5);
	  
	  i++;
	  inp->nImprTypes=i;
	}
      }
      else if(k==5)      
      {
	if((buff2!=NULL)&&(buff2[0]!='!'))
	{
	  buff3=strtok(NULL," \n\t");
	  buff4=strtok(NULL," \n\t");
	  buff5=strtok(NULL," \n\t");
	  
	  ii=-1;
	  for(l=0;l<inp->nTypes;l++)
	  {
	    if(!strcmp(buff2,inp->types[l]))
	    {
	      ii=l;
	      break;
	    }
	  }
	  
	  if(ii==-1)
	    continue;
	
	  inp->nonBondedTypesParm[ii][0]=atof(buff3);
	  inp->nonBondedTypesParm[ii][1]=atof(buff4);
	  inp->nonBondedTypesParm[ii][2]=atof(buff5);
	  
	  buff3=strtok(NULL," \n\t");
	  if((buff3!=NULL)&&(buff3[0]!='!'))
	  {
	    buff4=strtok(NULL," \n\t");
	    buff5=strtok(NULL," \n\t");
	    
	    inp->nonBondedTypesParm[ii][3]=atof(buff3);
	    inp->nonBondedTypesParm[ii][4]=atof(buff4);
	    inp->nonBondedTypesParm[ii][5]=atof(buff5);
	  }
	  i++;
	  inp->nNonBonded=i;
	}
      }
    }
    fclose(parFile);
}

void read_CONF(ATOM atom[],SIMULPARAMS *simulCond)
{
  
  char buff1[1024]="", *buff2=NULL;

    FILE *confFile=NULL;
    char ren[5],atl[5],sen[5];
    int i,atn,res,ire,natomCheck;
    double wei,xx,yy,zz;
    
    confFile=fopen("CONF","r");

    if (confFile==NULL)
    {
      error(40);
    }
    
    while(fgets(buff1,1024,confFile)!=NULL)
    {

      if(buff1[0]!='*')
	break;    
  
    }
    
    buff2=strtok(buff1," \n\t");
    natomCheck=atof(buff2);
    
    if(natomCheck!=simulCond->natom)
      error(41);
    
    for(i=0;i<natomCheck;i++)
    {
      fscanf(confFile,"%d %d %s %s %lf %lf %lf %s %d %lf",&atn,&ire,ren,atl,&xx,&yy,&zz,sen,&res,&wei);
      atom[i].x=xx;
      atom[i].y=yy;
      atom[i].z=zz;
      
      atom[i].ires=ire;

    }
    
    fclose(confFile);
}

void setup(INPUTS *inp,ATOM atom[],FORCEFIELD *ff,SIMULPARAMS *simulCond,CONSTRAINT *constList)
{
  int i,j,k,ia,ib,ic,id,i0,i1,i2,i3,itype;
  
  if(ff->nBond>0)
  {
    ff->parmBond=(double**)malloc(ff->nBond*sizeof(*(ff->parmBond)));
    for(i=0;i<ff->nBond;i++)
      ff->parmBond[i]=(double*)malloc(2*sizeof(**(ff->parmBond)));
  }
  
  if(ff->nAngle>0)
  {
    ff->parmAngle=(double**)malloc(ff->nAngle*sizeof(*(ff->parmAngle)));
    for(i=0;i<ff->nAngle;i++)
      ff->parmAngle[i]=(double*)malloc(3*sizeof(**(ff->parmAngle)));
  }
  
  if(ff->nDihedral>0)
  {
    simulCond->diheType=(int*)malloc(ff->nDihedral*sizeof(*(simulCond->diheType)));
    for(i=0;i<ff->nDihedral;i++)
      simulCond->diheType[i]=0;
    
    ff->parmDihe=(double**)malloc(ff->nDihedral*sizeof(*(ff->parmDihe)));
    
    ff->nParmDihe=(int*)malloc(ff->nDihedral*sizeof(*(ff->nParmDihe)));
  }
  
  if(ff->nImproper>0)
  {
    simulCond->imprType=(int*)malloc(ff->nImproper*sizeof(*(simulCond->imprType)));
    for(i=0;i<ff->nImproper;i++)
      simulCond->imprType[i]=0;
    
    ff->parmImpr=(double**)malloc(ff->nImproper*sizeof(*(ff->parmImpr)));
    for(i=0;i<ff->nImproper;i++)
      ff->parmImpr[i]=(double*)malloc(3*sizeof(**(ff->parmImpr)));
  }
  
  ff->parmVdw=(double**)malloc(simulCond->natom*sizeof(*(ff->parmVdw)));
  for(i=0;i<simulCond->natom;i++)
    ff->parmVdw[i]=(double*)malloc(6*sizeof(**(ff->parmVdw)));
  
  /*fprintf(outFile,"Let's check what there is in atom->atomType\n");
  for(i=0;i<simulCond->natom;i++)
  {
    fprintf(outFile,"%d %d %d\n",simulCond->natom,i,atom[i].type);
  }

  fprintf(outFile,"Let's check what there is in inp->types\n");
  for(i=0;i<inp->nTypes;i++)
  {
     fprintf(outFile,"%d %d %d %s\n",inp->nTypes,i,inp->typesNum[i],inp->types[i]);
  }
  
  fprintf(outFile,"Let's check what there is in iBond\n");
  for(i=0;i<ff->nBond;i++)
  {
     fprintf(outFile,"%d %d %d %d %d\n",ff->nBond,simulCond->iBond[i][0],simulCond->iBond[i][1],atom[simulCond->iBond[i][0]].type,atom[simulCond->iBond[i][1]].type);
  }
  
  fprintf(outFile,"Let's check what there is in bondTypes\n");
  for(j=0;j<inp->nBondTypes;j++)
  {
     fprintf(outFile,"%d %d %d\n",inp->nBondTypes,inp->bondTypes[j][0],inp->bondTypes[j][1]);
  }

  fprintf(outFile,"Let's check what there is in iAngle\n");
  for(i=0;i<ff->nAngle;i++)
  {
     fprintf(outFile,"%d %d %d %d\n",ff->nAngle,simulCond->iAngle[i][0],simulCond->iAngle[i][1],simulCond->iAngle[i][2]);
  }
  
  fprintf(outFile,"Let's check what there is in angTypes\n");
  for(j=0;j<inp->nAngTypes;j++)
  {
     fprintf(outFile,"%d %d %d %d\n",inp->nAngTypes,inp->angTypes[j][0],inp->angTypes[j][1],inp->angTypes[j][2]);
  }

  fprintf(outFile,"Let's check what there is in iDihedral\n");
  for(i=0;i<ff->nDihedral;i++)
  {
     fprintf(outFile,"%d %d %d %d %d\n",ff->nDihedral,simulCond->iDihedral[i][0],simulCond->iDihedral[i][1],simulCond->iDihedral[i][2],simulCond->iDihedral[i][3]);
  }
  
  fprintf(outFile,"Let's check what there is in diheTypes\n");
  for(j=0;j<inp->nDiheTypes;j++)
  {
     fprintf(outFile,"%d %d %d %d %d\n",inp->nDiheTypes,inp->diheTypes[j][0],inp->diheTypes[j][1],inp->diheTypes[j][2],inp->diheTypes[j][3]);
  }
  
  if(simulCond->keyconsth)
  {
    fprintf(outFile,"Let's check what there is in constList\n");
    fprintf(outFile,"pointer adress=%p\n",constList);
    for(i=0;i<simulCond->nconst;i++)
    {
      fprintf(outFile,"%d %d %d\n",simulCond->nconst,constList[i].a,constList[i].b);
    }
  }

  fprintf(outFile,"Now if the loop.\n");*/

  for(i=0;i<ff->nBond;i++)
  {
    ia=atom[simulCond->iBond[i][0]].type;
    ib=atom[simulCond->iBond[i][1]].type;
    
    if(ib<ia)
    {
      i0=ib;
      i1=ia;
    }
    else
    {
      i0=ia;
      i1=ib;
    }
    
    itype=-1;
    for(j=0;j<inp->nBondTypes;j++)
    {
      if( (i0==inp->bondTypes[j][0]) && (i1==inp->bondTypes[j][1]) )
      {
	itype=j;
	break;
      }
    }
    
    if(itype==-1)
      error(71);
    
    ff->parmBond[i][0]=inp->bondTypesParm[itype][0]*2.*kcaltoiu;
    ff->parmBond[i][1]=inp->bondTypesParm[itype][1];
    
  }
  
  if(simulCond->keyconsth)
  {
    for(i=0;i<simulCond->nconst;i++)
    {
      ia=atom[constList[i].a].type;
      ib=atom[constList[i].b].type;
      
      if(ib<ia)
      {
	i0=ib;
	i1=ia;
      }
      else
      {
	i0=ia;
	i1=ib;
      }
      
      itype=-1;
      for(j=0;j<inp->nBondTypes;j++)
      {
	if( (i0==inp->bondTypes[j][0]) && (i1==inp->bondTypes[j][1]) )
	{
	  itype=j;
	  break;
	}
      }
      
      if(itype==-1)
	error(71);
      
      constList[i].rc2=X2(inp->bondTypesParm[itype][1]);
    }
  }
  
  ff->nUb=0;
  for(i=0;i<ff->nAngle;i++)
  {
    ia=atom[simulCond->iAngle[i][0]].type;
    ib=atom[simulCond->iAngle[i][1]].type;
    ic=atom[simulCond->iAngle[i][2]].type;
    
    if(ic<ia)
    {
      i0=ic;
      i1=ib;
      i2=ia;
    }
    else
    {
      i0=ia;
      i1=ib;
      i2=ic;
    }
    
    itype=-1;
    for(j=0;j<inp->nAngTypes;j++)
    {
      if( (i0==inp->angTypes[j][0]) && (i1==inp->angTypes[j][1]) && (i2==inp->angTypes[j][2]) )
      {
	itype=j;
	break;
      }
    }
    
    if(itype==-1)
      error(72);
    
    ff->parmAngle[i][0]=inp->angTypesParm[itype][0]*2.*kcaltoiu;
    ff->parmAngle[i][1]=inp->angTypesParm[itype][1]*PI/180.;
    
    for(j=0;j<inp->nUbTypes;j++)
    {
      if( (i0==inp->ubTypes[j][0]) && (i1==inp->ubTypes[j][1]) && (i2==inp->ubTypes[j][2]) )
      {
	ff->nUb++;
	break;
      }
    }
  }
  
  if(ff->nUb>0)
  {
    ff->parmUb=(double**)malloc(ff->nUb*sizeof(*(ff->parmUb)));
    for(i=0;i<ff->nUb;i++)
      ff->parmUb[i]=(double*)malloc(2*sizeof(**(ff->parmUb)));
    
    simulCond->iUb=(int**)malloc(ff->nUb*sizeof(*(simulCond->iUb)));
    for(i=0;i<ff->nUb;i++)
      simulCond->iUb[i]=(int*)malloc(2*sizeof(**(simulCond->iUb)));
  }
  
  k=0;
  for(i=0;i<ff->nAngle;i++)
  {
    ia=atom[simulCond->iAngle[i][0]].type;
    ib=atom[simulCond->iAngle[i][1]].type;
    ic=atom[simulCond->iAngle[i][2]].type;
    
    if(ic<ia)
    {
      i0=ic;
      i1=ib;
      i2=ia;
    }
    else
    {
      i0=ia;
      i1=ib;
      i2=ic;
    }
    
    itype=-1;
    for(j=0;j<inp->nUbTypes;j++)
    {
      if( (i0==inp->ubTypes[j][0]) && (i1==inp->ubTypes[j][1]) && (i2==inp->ubTypes[j][2]) )
      {
	itype=j;
	break;
      }
    }
    
    if(itype==-1)
      continue;
    
    simulCond->iUb[k][0]=simulCond->iAngle[i][0];
    simulCond->iUb[k][1]=simulCond->iAngle[i][2];
    ff->parmUb[k][0]=inp->ubTypesParm[itype][0]*2.*kcaltoiu;
    ff->parmUb[k][1]=inp->ubTypesParm[itype][1];
    k++;
  }
  
  for(i=0;i<ff->nDihedral;i++)
  {
    ia=atom[simulCond->iDihedral[i][0]].type;
    ib=atom[simulCond->iDihedral[i][1]].type;
    ic=atom[simulCond->iDihedral[i][2]].type;
    id=atom[simulCond->iDihedral[i][3]].type;
    
    if(id<ia)
    {
      i0=id;
      i1=ic;
      i2=ib;
      i3=ia;
    }
    else if(ia<id)
    {
      i0=ia;
      i1=ib;
      i2=ic;
      i3=id;
    }
    else if(ic<ib)
    {
      i0=id;
      i1=ic;
      i2=ib;
      i3=ia;
    }
    else
    {
      i0=ia;
      i1=ib;
      i2=ic;
      i3=id;
    }
    
    itype=-1;
    for(j=0;j<inp->nDiheTypes;j++)
    {
      if( (i0==inp->diheTypes[j][0]) && (i1==inp->diheTypes[j][1]) && (i2==inp->diheTypes[j][2]) && (i3==inp->diheTypes[j][3]) )
      {
	itype=j;
	break;
      }
    }
    
    if(itype==-1)
    {
      if(ic<ib)
      {
	i0=inp->nTypes;
	i1=ic;
	i2=ib;
	i3=inp->nTypes;
      }
      else
      {
	i0=inp->nTypes;
	i1=ib;
	i2=ic;
	i3=inp->nTypes;
      }
    
      for(j=0;j<inp->nDiheTypes;j++)
      {
	if( (i0==inp->diheTypes[j][0]) && (i1==inp->diheTypes[j][1]) && (i2==inp->diheTypes[j][2]) && (i3==inp->diheTypes[j][3]) )
	{
	  itype=j;
	  break;
	}
      }
    }
    
    if(itype==-1)
      error(73);
    
    ff->nParmDihe[i]=inp->nDiheTypesParm[itype];
    
    ff->parmDihe[i]=(double*)malloc(3*inp->nDiheTypesParm[itype]*sizeof(**(ff->parmDihe)));
    
    for(j=0;j<ff->nParmDihe[i];j++)
    {
      k=3*j;
      ff->parmDihe[i][k]=inp->diheTypesParm[itype][k]*kcaltoiu;
      ff->parmDihe[i][k+1]=inp->diheTypesParm[itype][k+1];
      ff->parmDihe[i][k+2]=inp->diheTypesParm[itype][k+2]*PI/180.;
    }
    
    if(ff->nParmDihe[i]==1&&ff->parmDihe[i][1]<0.5)
    {
      ff->parmDihe[i][0]=ff->parmDihe[i][0]*2.;
      ff->parmDihe[i][1]=ff->parmDihe[i][2];
      simulCond->diheType[i]=2;
    }
    else if(ff->nParmDihe[i]>1&&ff->parmDihe[i][1]<0.5)
      error(50);
    else
      simulCond->diheType[i]=1;
  }
  
  for(i=0;i<ff->nImproper;i++)
  {
    ia=atom[simulCond->iImproper[i][0]].type;
    ib=atom[simulCond->iImproper[i][1]].type;
    ic=atom[simulCond->iImproper[i][2]].type;
    id=atom[simulCond->iImproper[i][3]].type;
    
    if(id<ia)
    {
      i0=id;
      i1=ic;
      i2=ib;
      i3=ia;
    }
    else if(ia<id)
    {
      i0=ia;
      i1=ib;
      i2=ic;
      i3=id;
    }
    else if(ic<ib)
    {
      i0=id;
      i1=ic;
      i2=ib;
      i3=ia;
    }
    else
    {
      i0=ia;
      i1=ib;
      i2=ic;
      i3=id;
    }
    
    itype=-1;
    for(j=0;j<inp->nImprTypes;j++)
    {
      if( (i0==inp->imprTypes[j][0]) && (i1==inp->imprTypes[j][1]) && (i2==inp->imprTypes[j][2]) && (i3==inp->imprTypes[j][3]) )
      {
	itype=j;
	break;
      }
    }
    
    if(itype==-1)
    {
      
      i0=id;
      i1=ic;
      i2=ib;
      i3=inp->nTypes;
	
      for(j=0;j<inp->nImprTypes;j++)
      {
	if( (i0==inp->imprTypes[j][0]) && (i1==inp->imprTypes[j][1]) && (i2==inp->imprTypes[j][2]) && (i3==inp->imprTypes[j][3]) )
	{
	  itype=j;
	  break;
	}
      }
    }
    
   if(itype==-1)
    {
      
      i0=ia;
      i1=ib;
      i2=ic;
      i3=inp->nTypes;
	
      for(j=0;j<inp->nImprTypes;j++)
      {
	if( (i0==inp->imprTypes[j][0]) && (i1==inp->imprTypes[j][1]) && (i2==inp->imprTypes[j][2]) && (i3==inp->imprTypes[j][3]) )
	{
	  itype=j;
	  break;
	}
      }
    }
    
    if(itype==-1)
    {
      if(ic<ib)
      {
	i0=inp->nTypes;
	i1=ic;
	i2=ib;
	i3=inp->nTypes;
      }
      else
      {
	i0=inp->nTypes;
	i1=ib;
	i2=ic;
	i3=inp->nTypes;
      }
    
      for(j=0;j<inp->nImprTypes;j++)
      {
	if( (i0==inp->imprTypes[j][0]) && (i1==inp->imprTypes[j][1]) && (i2==inp->imprTypes[j][2]) && (i3==inp->imprTypes[j][3]) )
	{
	  itype=j;
	  break;
	}
      }
    }
    
    if(itype==-1)
    {
      if(id<ia)
      {
	i0=id;
	i1=inp->nTypes;
	i2=inp->nTypes;
	i3=ia;
      }
      else
      {
	i0=ia;
	i1=inp->nTypes;
	i2=inp->nTypes;
	i3=id;
      }
    
      for(j=0;j<inp->nImprTypes;j++)
      {
	if( (i0==inp->imprTypes[j][0]) && (i1==inp->imprTypes[j][1]) && (i2==inp->imprTypes[j][2]) && (i3==inp->imprTypes[j][3]) )
	{
	  itype=j;
	  break;
	}
      }
    }
    
    if(itype==-1)
    {
      
      i0=id;
      i1=ic;
      i2=inp->nTypes;
      i3=inp->nTypes;
	
      for(j=0;j<inp->nImprTypes;j++)
      {
	if( (i0==inp->imprTypes[j][0]) && (i1==inp->imprTypes[j][1]) && (i2==inp->imprTypes[j][2]) && (i3==inp->imprTypes[j][3]) )
	{
	  itype=j;
	  break;
	}
      }
    }
    
   if(itype==-1)
    {
      
      i0=ia;
      i1=ib;
      i2=inp->nTypes;
      i3=inp->nTypes;
	
      for(j=0;j<inp->nImprTypes;j++)
      {
	if( (i0==inp->imprTypes[j][0]) && (i1==inp->imprTypes[j][1]) && (i2==inp->imprTypes[j][2]) && (i3==inp->imprTypes[j][3]) )
	{
	  itype=j;
	  break;
	}
      }
    }
    
    if(itype==-1)
      error(74);
    
    ff->parmImpr[i][0]=inp->imprTypesParm[itype][0]*kcaltoiu;
    ff->parmImpr[i][1]=inp->imprTypesParm[itype][1];
    ff->parmImpr[i][2]=inp->imprTypesParm[itype][2]*PI/180.;
    
    if(ff->parmImpr[i][1]<0.5)
    {
      ff->parmImpr[i][0]=ff->parmImpr[i][0]*2.;
      ff->parmImpr[i][1]=ff->parmImpr[i][2];
      simulCond->imprType[i]=2;
    }
    else
      simulCond->imprType[i]=1;
  }
  
  for(i=0;i<simulCond->natom;i++)
  {
    
    ff->parmVdw[i][0]=sqrt(-kcaltoiu*inp->nonBondedTypesParm[atom[i].type][1]);
    ff->parmVdw[i][1]=inp->nonBondedTypesParm[atom[i].type][2]/sq6rt2;
    ff->parmVdw[i][2]=inp->nonBondedTypesParm[atom[i].type][0];
    
    if(inp->nonBondedTypesParm[atom[i].type][5]>=0.)
    {
      ff->parmVdw[i][3]=sqrt(-kcaltoiu*inp->nonBondedTypesParm[atom[i].type][4]);
      ff->parmVdw[i][4]=inp->nonBondedTypesParm[atom[i].type][5]/sq6rt2;
      ff->parmVdw[i][5]=inp->nonBondedTypesParm[atom[i].type][3];
    }
    else
    {
      ff->parmVdw[i][3]=ff->parmVdw[i][0];
      ff->parmVdw[i][4]=ff->parmVdw[i][1];
      ff->parmVdw[i][5]=ff->parmVdw[i][2];
    }
  }
  
}

void write_CONF(ATOM atom[],SIMULPARAMS *simulCond)
{
  
  FILE *confFile=NULL;
  int i;
  double wei=0.;

  confFile=fopen("RESCONF","w");

  if (confFile==NULL)
  {
    error(47);
  }

  fprintf(confFile,"*Restart Configuration\n");
  fprintf(confFile,"* Written by MDBas\n");
  fprintf(confFile,"*\n");
  fprintf(confFile,"%5d\n",simulCond->natom);
  
  for(i=0;i<simulCond->natom;i++)
  {
    fprintf(confFile,"%5d %4d %-4s %-4s%10.5lf%10.5lf%10.5lf %-4s %-4d%10.5lf\n",
	    i+1,atom[i].ires,atom[i].resn,atom[i].label,atom[i].x,atom[i].y,
	    atom[i].z,atom[i].segn,atom[i].resi,wei);
    
  }

  fclose(confFile);
}

void write_prop(SIMULPARAMS *simulCond,ENERGY *ener,PBC *box)
{
  FILE *propFile=fopen("PROP","ab");
  
  double buffer[29];
  double temp,press;
  
  temp=2.*ener->kin/((double)simulCond->degfree*rboltzui);
  if(box->type>0)
    press=(2.*ener->kin-ener->virtot)/(3.*box->vol*bartoiu);
  else
    press=0.;

  buffer[0]=(double)simulCond->step;
  buffer[1]=(double)simulCond->step*simulCond->timeStep;
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
  
  fwrite(buffer,sizeof(double),29,propFile);
  
  fclose(propFile);
  
}

void write_rest(SIMULPARAMS *simulCond,ENERGY *ener,ATOM *atom)
{
  FILE *restFile=fopen("RESTART","wb");
  
  /*/double buffer[20];
  
  fwrite(&(simulCond->step),sizeof(int),1,restFile);

  buffer[0]=ener->tot;
  buffer[1]=ener->kin;
  buffer[2]=ener->pot;
  buffer[3]=ener->elec;
  buffer[4]=ener->vdw;
  buffer[5]=ener->bond;
  buffer[6]=ener->ang;
  buffer[7]=ener->ub;
  buffer[8]=ener->dihe;
  buffer[9]=ener->impr;
  buffer[10]=ener->virtot;
  buffer[11]=ener->virpot;
  buffer[12]=ener->virelec;
  buffer[13]=ener->virvdw;
  buffer[14]=ener->virbond;
  buffer[15]=ener->virub;
  buffer[16]=ener->virshake;
  buffer[17]=ener->conint;
  buffer[18]=simulCond->lambdat;
  buffer[19]=simulCond->lambdat;
  
  fwrite(buffer,sizeof(double),20,propFile);
  
  int i;
  double *x,*y,*z;
  double *vx,*vy,*vz;
  double *fx,*fy,*fz;
  
  
  
  for(i=0;i<simulCond;i++)
  {
    x[i]=
  }*/
  
  size_t ret;
  
  ret=fwrite(&(simulCond->step),sizeof(int),1,restFile);
  if(ret!=1)
    error(501);
  
  ret=fwrite(&(simulCond->natom),sizeof(double),1,restFile);
  if(ret!=1)
    error(501);
  
  ret=fwrite(atom,sizeof(ATOM),simulCond->natom,restFile);
  if(ret!=simulCond->natom)
    error(501);
  
  ret=fwrite(ener,sizeof(ENERGY),1,restFile);
  if(ret!=1)
    error(501);
  
  ret=fwrite(&(simulCond->lambdat),sizeof(double),1,restFile);
  if(ret!=1)
    error(501);
  
  ret=fwrite(&(simulCond->gammap),sizeof(double),1,restFile);
  if(ret!=1)
    error(501);
  
  fclose(restFile);
  
}

void read_rest(SIMULPARAMS *simulCond,ENERGY *ener,ATOM *atom)
{
  FILE *restFile=fopen("RESTART","rb");
  
  size_t ret;
  
  ret=fread(&(simulCond->step),sizeof(int),1,restFile);
  if(ret!=1)
    error(502);
  
  ret=fread(&(simulCond->natom),sizeof(double),1,restFile);
  if(ret!=1)
    error(502);
  
  ret=fread(atom,sizeof(ATOM),simulCond->natom,restFile);
  if(ret!=simulCond->natom)
    error(502);
  
  ret=fread(ener,sizeof(ENERGY),1,restFile);
  if(ret!=1)
    error(502);
  
  ret=fread(&(simulCond->lambdat),sizeof(double),1,restFile);
  if(ret!=1)
    error(502);
  
  ret=fread(&(simulCond->gammap),sizeof(double),1,restFile);
  if(ret!=1)
    error(502);
  
  fclose(restFile);
  
}

void write_FORF(INPUTS *inp,ATOM atom[],FORCEFIELD *ff,SIMULPARAMS *simulCond)
{
  FILE *forfFile;
  int i,j,k,l,ia,ib,ic,id,nd;
  double kf,r0,ncos;
  
  forfFile=fopen("FORF","w");
  
  fprintf(forfFile,"ATOMS %d\n",simulCond->natom);
  
  for(i=0;i<simulCond->natom;i++)
  {
    k=atom[i].type;
    fprintf(forfFile,"%s %lf %lf 1\n",inp->types[k],atom[i].m,atom[i].q);
  }
  
  fprintf(forfFile,"BONDS %d\n",ff->nBond+ff->nUb);
  
  for(i=0;i<ff->nBond;i++)
  {
    ia=simulCond->iBond[i][0]+1;
    ib=simulCond->iBond[i][1]+1;
    kf=ff->parmBond[i][0]/kcaltoiu;
    r0=ff->parmBond[i][1];
    
    fprintf(forfFile,"harm %d %d %lf %lf\n",ia,ib,kf,r0);
       
  }
  
  for(i=0;i<ff->nUb;i++)
  {
    ia=simulCond->iUb[i][0]+1;
    ib=simulCond->iUb[i][1]+1;
    kf=ff->parmUb[i][0]/kcaltoiu;
    r0=ff->parmUb[i][1];
    
    fprintf(forfFile,"harm %d %d %lf %lf\n",ia,ib,kf,r0);
       
  }
  
  fprintf(forfFile,"ANGLES %d\n",ff->nAngle);
  
  for(i=0;i<ff->nAngle;i++)
  {
    ia=simulCond->iAngle[i][0]+1;
    ib=simulCond->iAngle[i][1]+1;
    ic=simulCond->iAngle[i][2]+1;
    kf=ff->parmAngle[i][0]/kcaltoiu;
    r0=ff->parmAngle[i][1]*180./PI;
    
    fprintf(forfFile,"harm %d %d %d %lf %lf\n",ia,ib,ic,kf,r0);
  }
  
  nd=0;
  for(i=0;i<ff->nDihedral;i++)
  {
    for(j=0;j<ff->nParmDihe[i];j++)
    {
      nd++;
    }
  }
  
  fprintf(forfFile,"DIHEDRALS %d\n",nd+ff->nImproper);
  
  for(i=0;i<ff->nDihedral;i++)
  {
    ia=simulCond->iDihedral[i][0]+1;
    ib=simulCond->iDihedral[i][1]+1;
    ic=simulCond->iDihedral[i][2]+1;
    id=simulCond->iDihedral[i][3]+1;
    
    if(simulCond->diheType[i]==1)
    {
      for(j=0;j<ff->nParmDihe[i];j++)
      {
	k=3*j;
	kf=ff->parmDihe[i][k]/kcaltoiu;
	ncos=ff->parmDihe[i][k+1];
	r0=ff->parmDihe[i][k+2]*180./PI;
	
	fprintf(forfFile,"cos %d %d %d %d %lf %lf %lf\n",ia,ib,ic,id,kf,r0,ncos);
      }
    }
    else if(simulCond->diheType[i]==2)
    {
      kf=ff->parmDihe[i][0]/kcaltoiu;
      r0=ff->parmDihe[i][1]*180./PI;
      
      fprintf(forfFile,"harm %d %d %d %d %lf %lf\n",ia,ib,ic,id,kf,r0);
    }
    
  }
  
  for(i=0;i<ff->nImproper;i++)
  {
    ia=simulCond->iImproper[i][0]+1;
    ib=simulCond->iImproper[i][1]+1;
    ic=simulCond->iImproper[i][2]+1;
    id=simulCond->iImproper[i][3]+1;
    
    kf=ff->parmDihe[i][0]/kcaltoiu;
    r0=ff->parmDihe[i][1]*180./PI;
    
    fprintf(forfFile,"harm %d %d %d %d %lf %lf\n",ia,ib,ic,id,kf,r0);
  }
  
  nd=simulCond->natom*(simulCond->natom+1)/2;
  fprintf(forfFile,"VDW %d\n",nd);
  
  for(i=0;i<simulCond->natom;i++)
  {
    for(j=i;j<simulCond->natom;j++)
    {
      k=atom[i].type;
      l=atom[j].type;
      kf=ff->parmVdw[i][0]*ff->parmVdw[j][0]/kcaltoiu;
      r0=ff->parmVdw[i][1]+ff->parmVdw[j][1];
      
      fprintf(forfFile,"%s %s %lf %lf 1\n",inp->types[k],inp->types[l],kf,r0);
    }
  }
  
  fclose(forfFile);
}

void write_DCD_header(SIMULPARAMS *simulCond, PBC *box)
{
    FILE *dcdf=fopen("TRAJ","wb");
    
    char HDR[4]={'C','O','R','D'};
    
    int ICNTRL[20]={0};
    ICNTRL[0]= simulCond->nsteps/simulCond->printtr;
    ICNTRL[1]= simulCond->printtr;
    ICNTRL[2]= simulCond->printtr;
    ICNTRL[3]= simulCond->nsteps;
    ICNTRL[7]= simulCond->degfree;
    ICNTRL[8]= 0; //no frozen atom
    //ICNTRL[9]= timestep in akma but in 32 bits mode
    ICNTRL[10]= (box->type == NOBOX)?0:1;
    ICNTRL[19]= 37;
    
    int NATOM=simulCond->natom;
    
    int NTITLE=1;
    char TITLE[80]="";
    
    time_t rawtime;
    struct tm *timeinfo;

    time(&rawtime);
    timeinfo = localtime(&rawtime);
  
    sprintf(TITLE,"* DCD writing time : %s",asctime(timeinfo));
    
    unsigned int size;
    
    // first : HDR + ICNTRL
    size = 4*sizeof(char)+20*sizeof(int);
    fwrite(&size,sizeof(unsigned int),1,dcdf);
    fwrite(HDR,sizeof(char),4,dcdf);
    fwrite(ICNTRL,sizeof(int),20,dcdf);
    fwrite(&size,sizeof(unsigned int),1,dcdf);
    
    // second : NTITLE + TITLE
    size = sizeof(int) + 80*sizeof(char);
    fwrite(&size,sizeof(unsigned int),1,dcdf);
    fwrite(&NTITLE,sizeof(int),1,dcdf);
    fwrite(TITLE,sizeof(char),80,dcdf);
    fwrite(&size,sizeof(unsigned int),1,dcdf);
    
    // third : natom
    size = sizeof(int);
    fwrite(&size,sizeof(unsigned int),1,dcdf);
    fwrite(&NATOM,sizeof(int),1,dcdf);
    fwrite(&size,sizeof(unsigned int),1,dcdf);
    
    fclose(dcdf);
}

void write_DCD_traj(ATOM atom[], SIMULPARAMS *simulCond, PBC *box)
{
    FILE *dcdf=fopen("TRAJ","ab");
    
    unsigned int size;
    
    // if some PBC write a part of the matrix
    if (box->type != NOBOX)
    {
        double box_matrix[6]; 
        box_to_crystal(box,box_matrix);
        
        size = 6*sizeof(double);
        fwrite(&size,sizeof(unsigned int),1,dcdf);
        fwrite(box_matrix,sizeof(double),6,dcdf);
        fwrite(&size,sizeof(unsigned int),1,dcdf);
    }
    
    // alloc of crdinates array + loading of coordinates ; but in float mode !
    float *X = (float*)malloc(simulCond->natom*sizeof(float));
    float *Y = (float*)malloc(simulCond->natom*sizeof(float));
    float *Z = (float*)malloc(simulCond->natom*sizeof(float));
    
    int i;
    for(i=0;i<simulCond->natom;i++)
    {
        X[i] = (float) atom[i].x;
        Y[i] = (float) atom[i].y;
        Z[i] = (float) atom[i].z;
    }
    
    size = simulCond->natom*sizeof(float);
    
    // writing X coordinates
    fwrite(&size,sizeof(unsigned int),1,dcdf);
    fwrite(X,sizeof(float),simulCond->natom,dcdf);
    fwrite(&size,sizeof(unsigned int),1,dcdf);
    
    // writing Y coordinates
    fwrite(&size,sizeof(unsigned int),1,dcdf);
    fwrite(Y,sizeof(float),simulCond->natom,dcdf);
    fwrite(&size,sizeof(unsigned int),1,dcdf);
    
    // writing Z coordinates
    fwrite(&size,sizeof(unsigned int),1,dcdf);
    fwrite(Z,sizeof(float),simulCond->natom,dcdf);
    fwrite(&size,sizeof(unsigned int),1,dcdf);
    
    fclose(dcdf);
}

void free_temp_array(INPUTS *inp)
{
  free(inp->typesNum);
  free(inp->nDiheTypesParm);
  free_2D(inp->nBondTypes,inp->bondTypes,inp->bondTypesParm,NULL);
  free_2D(inp->nAngTypes,inp->angTypes,inp->angTypesParm,NULL);
  free_2D(inp->nUbTypes,inp->ubTypes,inp->ubTypesParm,NULL);
  free_2D(inp->nDiheTypes,inp->diheTypes,inp->diheTypesParm,NULL);
  free_2D(inp->nImprTypes,inp->imprTypes,inp->imprTypesParm,NULL);
  free_2D(inp->nTypes,inp->nonBondedTypesParm,NULL);
}

void error(int errorNumber)
{
  fprintf(outFile,"MDBas failed due to error number: %d\n",errorNumber);
  switch (errorNumber)
  {
  case 10:
    fprintf(outFile,"MDBas cannot find or open topology file TOP.\n");
    fprintf(outFile,"Most likely, it is not properly named. Please check.\n");
    break;
  case 20:
    fprintf(outFile,"MDBas cannot find or open structure file PSF.\n");
    fprintf(outFile,"Most likely, it is not properly named. Please check.\n");
    break;
  case 22:
    fprintf(outFile,"MDBas encountered a problem while reading atomic properties\n");
    fprintf(outFile,"in the PSF file. There is an unexpected line there. Please\n");
    fprintf(outFile,"consult the manual for further details about PSF file\n");
    break;
  case 23:
    fprintf(outFile,"There is problem in bonds sequence in the PSF file. Please\n");
    fprintf(outFile,"consult the manual for further details about PSF file\n");
    break;
  case 24:
    fprintf(outFile,"There is problem in angles sequence in the PSF file. Please\n");
    fprintf(outFile,"consult the manual for further details about PSF file\n");
    break;
  case 25:
    fprintf(outFile,"There is problem in dihedrals sequence in the PSF file. Please\n");
    fprintf(outFile,"consult the manual for further details about PSF file\n");
    break;
  case 26:
    fprintf(outFile,"There is problem in improper angles sequence in the PSF file.\n");
    fprintf(outFile,"Please consult the manual for further details about PSF file\n");
    break;
  case 30:
    fprintf(outFile,"MDBas cannot find or open parameter file PAR.\n");
    fprintf(outFile,"Most likely, it is not properly named. Please check.\n");
    break;
  case 40:
    fprintf(outFile,"MDBas cannot find or open configuration file CONF.\n");
    fprintf(outFile,"Most likely, it is not properly named. Please check.\n");
    break;
  case 41:
    fprintf(outFile,"MDBas found a different number of atoms in CONF file and\n");
    fprintf(outFile,"in PSF file. Structure does not match configuration.\n");
    fprintf(outFile,"Check carefully these files.\n");
    break;
  case 47:
    fprintf(outFile,"MDBas cannot open configuration file RESCONF.\n");
    break;
  case 50:
    fprintf(outFile,"A dihedral angle is specified as a Fourier series but\n");
    fprintf(outFile,"with one of the component being an harmonic potential.\n");
    fprintf(outFile,"Check in PAR file.\n");
    break;
  case 60:
    fprintf(outFile,"MDBas cannot find or open simulation file SIMU.\n");
    fprintf(outFile,"Most likely, it is not properly named. Please check.\n");
    break;
  case 61:
    fprintf(outFile,"MDBas does not recognise a keyword specified in SIMU.\n");
    fprintf(outFile,"Please check SIMU file and the manual for the list of\n");
    fprintf(outFile,"allowed keywords.\n");
    break;
  case 62:
    fprintf(outFile,"MDBas does not recognise a parameter specified in SIMU.\n");
    fprintf(outFile,"Please check SIMU file and the manual for the list of\n");
    fprintf(outFile,"allowed keywords and their associated parameters.\n");
    break;
  case 63:
    fprintf(outFile,"MDBas does not find a required parameter in SIMU.\n");
    fprintf(outFile,"Please check SIMU file and the manual for the list of\n");
    fprintf(outFile,"allowed keywords and their associated parameters.\n");
    break;
  case 71:
    fprintf(outFile,"There is an undefined bond in the PSF. Most likely,\n");
    fprintf(outFile,"there are missing parameters in the PAR file. Please check\n");
    break;
  case 72:
    fprintf(outFile,"There is an undefined angle in the PSF. Most likely,\n");
    fprintf(outFile,"there are missing parameters in the PAR file. Please check\n");
    break;
  case 73:
    fprintf(outFile,"There is an undefined dihedral angle in the PSF. Most likely,\n");
    fprintf(outFile,"there are missing parameters in the PAR file. Please check\n");
    break;
  case 74:
    fprintf(outFile,"There is an undefined improper angle in the PSF. Most likely,\n");
    fprintf(outFile,"there are missing parameters in the PAR file. Please check\n");
    break;
  case 110:
    fprintf(outFile,"MDBas found a too many non-parameterised dihedral angles:\n");
    fprintf(outFile,"4*nDihedrals. nDihedrals comes from the value specified\n");
    fprintf(outFile,"in PSF file. Please check in PAR file. If such a number is\n");
    fprintf(outFile,"normal for your simulation, you have to enter list.c to\n");
    fprintf(outFile,"increase the size of the 1-4 pairs array from 5*nDihedrals to\n");
    fprintf(outFile,"the size you really need. Then recompile MDBas.\n");
    break;
  case 111:
    fprintf(outFile,"MDBas encountered a problem while setting the excluded atoms\n");
    fprintf(outFile,"list. The last atom has exclusion which should not happen. This\n");
    fprintf(outFile,"a bit annoying for there is no simple explanation for this.\n");
    fprintf(outFile,"Maybe an error in one of the input files which is not detected\n");
    fprintf(outFile,"by MDBas. Sorry for the trouble.\n");
    break;
  case 112:
    fprintf(outFile,"MDBas encountered a problem while setting the excluded atoms\n");
    fprintf(outFile,"list. The total excluded atoms does not match of the sum of\n");
    fprintf(outFile,"excluded atoms for each atom. This a bit annoying for there is\n");
    fprintf(outFile,"no simple explanation for this. Maybe an error in one of the\n");
    fprintf(outFile,"input files which is not detected by MDBas. Sorry for the trouble.\n");
    break;
  case 120:
    fprintf(outFile,"The number of neighbours around an atom is larger than the\n");
    fprintf(outFile,"maximum allocated memory in the verList array (default 2048).\n");
    fprintf(outFile,"This can occur fo very large cutoffs or in very heterogeneous\n");
    fprintf(outFile,"systems. This can be fixed by changing the variable MAXLIST in\n");
    fprintf(outFile,"global.h and recompilation of MDBas. An easy way to achieve this\n");
    fprintf(outFile,"without editing the header file, is to add the compilation option\n");
    fprintf(outFile,"-DMAXLIST=X, X being the new size of the array, for example 3072.\n");
    break;
  case 201:
    fprintf(outFile,"Unknown electrostatic potential. This is most likely due to an\n");
    fprintf(outFile,"error in the SIMU file. Please check this file and the manual\n");
    fprintf(outFile,"for the list of keywords and available potentials.\n");
    break;
  case 202:
    fprintf(outFile,"Unknown van der Waals potential. This is most likely due to an\n");
    fprintf(outFile,"error in the SIMU file. Please check this file and the manual\n");
    fprintf(outFile,"for the list of keywords and available potentials.\n");
    break;
  case 310:
    fprintf(outFile,"Velocities quenching convergence failure, most likely due a non suitable\n");
    fprintf(outFile,"initial configuration. If not, you can try increasing the number of cycles or\n");
    fprintf(outFile,"make Shake convergence criterion more tolerant. Please check the manual.\n");
    break;
  case 311:
    fprintf(outFile,"Shake convergence failure, most likely due a non suitable initial\n");
    fprintf(outFile,"configuration. If not, you can try increasing the number of cycles or\n");
    fprintf(outFile,"make Shake convergence criterion more tolerant. Please check the manual.\n");
    break;
  case 411:
    fprintf(outFile,"Too large ratio of the cutoff for the link cell is asked. Maximum is 5.\n");
    break;
  case 412:
    fprintf(outFile,"Too few link cells are created compared to the required ratio of the\n");
    fprintf(outFile,"cutoff. Please check first that you cutoff is smaller than half of the\n");
    fprintf(outFile,"smallest lattice parameter of your simulation. Alternatively, you can\n");
    fprintf(outFile,"change the ratio or force the use of the standard neighbour list algorithm.\n");
    break;
  default:
    fprintf(outFile,"MDBas failed due to unknown error number: %d\n",errorNumber);
    fprintf(outFile,"Reading the manual will not help you. You are by yourself.\n");
    fprintf(outFile,"Errare humanum est.\n");
    break;
  }
  
  exit(errorNumber);
    
}
