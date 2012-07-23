#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "global.h"
#include "io.h"
#include "memory.h"
#include "utils.h"

void read_SIMU(SIMULPARAMS *simulCond,FORCEFIELD *ff)
{
  char buff1[1024]="", *buff2=NULL, *buff3=NULL, *buff4=NULL;
  FILE *simuFile;
  
  simuFile=fopen("SIMU","r");
  
  if (simuFile==NULL)
  {
    error(60);
  }
  
  simulCond->keymd=1;
  
  simulCond->step=0;
  simulCond->firstener=1;
  
  simulCond->mdNature=1;
  simulCond->timeStep=0.001;
  simulCond->nsteps=0;
 
  simulCond->cutoff=12.0;
  simulCond->cuton=10.0;
  simulCond->delr=2.0;

  simulCond->temp=300.0;
  
  simulCond->keyrand=0;
  simulCond->seed=12345;
  
  simulCond->elecType=1;
  simulCond->vdwType=1;
  simulCond->nb14=0;
  ff->scal14=1.0;
  simulCond->numDeriv=0;

  simulCond->integrator=1;
  simulCond->ens=0;
  simulCond->enstime=1.0;
  simulCond->tolshake=1e-8;
  simulCond->maxcycle=100;
  simulCond->keyconsth=0;
  simulCond->nconst=0;
  
  simulCond->keyener=0;
  simulCond->keytraj=0;
  simulCond->keyforf=0;
  simulCond->printo=1000;
  simulCond->printtr=1000;
  
  simulCond->periodicType=0;
  simulCond->periodicBox[0][0]=0.;
  simulCond->periodicBox[1][1]=0.;
  simulCond->periodicBox[2][2]=0.;
  
  while(fgets(buff1,1024,simuFile)!=NULL)
  {
    buff2=strtok(buff1," \n\t");
    
    if(buff2==NULL)
      continue;
    
    nocase(buff2);
    puts(buff2);
    
    if(!strcmp(buff2,"mdbas"))
      simulCond->mdNature=0;
    
    else if(!strcmp(buff2,"charmm"))
      simulCond->mdNature=1;
    
    else if(!strcmp(buff2,"namd"))
      simulCond->mdNature=2;
    
    else if(!strcmp(buff2,"nomd"))
      simulCond->keymd=0;
    
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
	simulCond->elecType=0;
      else if(!strcmp(buff3,"full"))
	simulCond->elecType=1;
      else if(!strcmp(buff3,"shift1"))
	simulCond->elecType=2;
      else if(!strcmp(buff3,"shift2"))
	simulCond->elecType=3;
      else if(!strcmp(buff3,"switch"))
	simulCond->elecType=4;
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
	simulCond->vdwType=0;
      else if(!strcmp(buff3,"full"))
	simulCond->vdwType=1;
      else if(!strcmp(buff3,"switch"))
	simulCond->vdwType=2;
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
      else if(!strcmp(buff3,"nvt"))
	simulCond->ens=1;
      else if(!strcmp(buff3,"npt"))
	simulCond->ens=3;
      else
	error(62);
      
      if(simulCond->ens>0)
      {
	buff3=strtok(NULL," \n\t");
	if(buff3==NULL)
	error(63);
	
	simulCond->enstime=atof(buff3);
      }
    }
    else if(!strcmp(buff2,"temperature"))
    {
      buff3=strtok(NULL," \n\t");
      if(buff3==NULL)
	error(63);
      
      simulCond->temp=atof(buff3);
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
    else if(!strcmp(buff2,"ener"))
    {
      buff3=strtok(NULL," \n\t");
      if(buff3==NULL)
	error(63);
      
      simulCond->keyener=1;
      simulCond->printo=atoi(buff3);
    }
    else if(!strcmp(buff2,"traj"))
    {
      buff3=strtok(NULL," \n\t");
      if(buff3==NULL)
	error(63);
      
      simulCond->keytraj=1;
      simulCond->printtr=atoi(buff3);
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
      
      simulCond->periodicType=atoi(buff3);
      
      if(fgets(buff1,1024,simuFile)!=NULL)
      {
	
	buff2=strtok(buff1," \n\t");
	if(buff2==NULL)
	  error(63);
	else if(isdigit(buff2[0])==0)
	  error(62);
	
	simulCond->periodicBox[0][0]=atof(buff2);
	
	buff2=strtok(buff1," \n\t");
	if(buff2==NULL)
	  error(63);
	else if(isdigit(buff2[0])==0)
	  error(62);
	
	simulCond->periodicBox[0][1]=atof(buff2);
	
	buff2=strtok(buff1," \n\t");
	if(buff2==NULL)
	  error(63);
	else if(isdigit(buff2[0])==0)
	  error(62);
	
	simulCond->periodicBox[0][2]=atof(buff2);
	
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
	
	simulCond->periodicBox[1][0]=atof(buff2);
	
	buff2=strtok(buff1," \n\t");
	if(buff2==NULL)
	  error(63);
	else if(isdigit(buff2[0])==0)
	  error(62);
	
	simulCond->periodicBox[1][1]=atof(buff2);
	
	buff2=strtok(buff1," \n\t");
	if(buff2==NULL)
	  error(63);
	else if(isdigit(buff2[0])==0)
	  error(62);
	
	simulCond->periodicBox[1][2]=atof(buff2);
	
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
	
	simulCond->periodicBox[2][0]=atof(buff2);
	
	buff2=strtok(buff1," \n\t");
	if(buff2==NULL)
	  error(63);
	else if(isdigit(buff2[0])==0)
	  error(62);
	
	simulCond->periodicBox[2][1]=atof(buff2);
	
	buff2=strtok(buff1," \n\t");
	if(buff2==NULL)
	  error(63);
	else if(isdigit(buff2[0])==0)
	  error(62);
	
	simulCond->periodicBox[2][2]=atof(buff2);
	
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

void read_PSF(INPUTS *inp,ATOM *atom,FORCEFIELD *ff,SIMULPARAMS *simulCond,CONSTRAINT *constList)
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
      atom->natom=atoi(buff2);
      break;
    }
  }
  
  inp->typesNum=(int*)malloc(nalloc*sizeof(*(inp->typesNum)));

  atom->x=(double*)malloc(atom->natom*sizeof(*(atom->x)));
  atom->y=(double*)malloc(atom->natom*sizeof(*(atom->y)));
  atom->z=(double*)malloc(atom->natom*sizeof(*(atom->z)));
    
  atom->vx=(double*)malloc(atom->natom*sizeof(*(atom->vx)));
  atom->vy=(double*)malloc(atom->natom*sizeof(*(atom->vy)));
  atom->vz=(double*)malloc(atom->natom*sizeof(*(atom->vz)));
  
  atom->fx=(double*)malloc(atom->natom*sizeof(*(atom->fx)));
  atom->fy=(double*)malloc(atom->natom*sizeof(*(atom->fy)));
  atom->fz=(double*)malloc(atom->natom*sizeof(*(atom->fz)));
  
  atom->atomType=(int*)malloc(atom->natom*sizeof(*(atom->atomType)));
  atom->resn=(int*)malloc(atom->natom*sizeof(*(atom->resn)));
  atom->m=(double*)malloc(atom->natom*sizeof(*(atom->m)));
  
  atom->atomLabel=(char**)malloc(atom->natom*sizeof(*(atom->atomLabel)));
  atom->segi=(char**)malloc(atom->natom*sizeof(*(atom->segi)));
  atom->resi=(char**)malloc(atom->natom*sizeof(*(atom->resi)));
    
  for(i=0;i<atom->natom;i++)
  {
    atom->atomLabel[i]=(char*)malloc(5*sizeof(**(atom->atomLabel)));
    atom->segi[i]=(char*)malloc(5*sizeof(**(atom->segi)));
    atom->resi[i]=(char*)malloc(5*sizeof(**(atom->resi)));
  }
  
  ff->q=(double*)malloc(atom->natom*sizeof(*(ff->q)));
  
  k=-1;
  for(i=0;i<atom->natom;i++)
  {
    if(fgets(buff1,1024,psfFile)!=NULL)
    {
      buff2=strtok(buff1," \n\t");
      buff3=strtok(NULL," \n\t");
      buff4=strtok(NULL," \n\t");
      buff5=strtok(NULL," \n\t");

      strcpy(atom->segi[i],buff3);
      atom->resn[i]=atoi(buff4);
      strcpy(atom->resi[i],buff5);
      
      buff2=strtok(NULL," \n\t");
      buff3=strtok(NULL," \n\t");
      buff4=strtok(NULL," \n\t");
      buff5=strtok(NULL," \n\t");
      
      strcpy(atom->atomLabel[i],buff2);
      
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

      atom->atomType[i]=kk;
      ff->q[i]=atof(buff4);
      atom->m[i]=atof(buff5);
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
    atom->inconst=(int*)malloc(atom->natom*sizeof(*(atom->inconst)));
    for(i=0;i<atom->natom;i++)
      atom->inconst[i]=0;
    
    printf("pointer adress=%p\n",constList);
  
    constList=(CONSTRAINT*)malloc(ff->nBond*sizeof(*constList));
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
	if(atom->atomLabel[ia][0]=='H'||atom->atomLabel[ib][0]=='H')
	{
	  constList[simulCond->nconst].a=ia;
	  constList[simulCond->nconst].b=ib;
	  
	  atom->inconst[ia]++;
	  atom->inconst[ib]++;
	  
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
    constList=(CONSTRAINT*)realloc(constList,simulCond->nconst*sizeof(*constList));
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

void read_CONF(ATOM *atom)
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
    
    if(natomCheck!=atom->natom)
      error(41);
    
    for(i=0;i<natomCheck;i++)
    {
      fscanf(confFile,"%d %d %s %s %lf %lf %lf %s %d %lf",&atn,&ire,ren,atl,&xx,&yy,&zz,sen,&res,&wei);
      atom->x[i]=xx;
      atom->y[i]=yy;
      atom->z[i]=zz;

    }
    
    fclose(confFile);
}

void setup(INPUTS *inp,ATOM *atom,FORCEFIELD *ff,SIMULPARAMS *simulCond,CONSTRAINT *constList)
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
  
  ff->parmVdw=(double**)malloc(atom->natom*sizeof(*(ff->parmVdw)));
  for(i=0;i<atom->natom;i++)
    ff->parmVdw[i]=(double*)malloc(6*sizeof(**(ff->parmVdw)));
  
  printf("Let's check what there is in atom->atomType\n");
  for(i=0;i<atom->natom;i++)
  {
    printf("%d %d %d\n",atom->natom,i,atom->atomType[i]);
  }

  printf("Let's check what there is in inp->types\n");
  for(i=0;i<inp->nTypes;i++)
  {
     printf("%d %d %d %s\n",inp->nTypes,i,inp->typesNum[i],inp->types[i]);
  }
  
  printf("Let's check what there is in iBond\n");
  for(i=0;i<ff->nBond;i++)
  {
     printf("%d %d %d %d %d\n",ff->nBond,simulCond->iBond[i][0],simulCond->iBond[i][1],atom->atomType[simulCond->iBond[i][0]],atom->atomType[simulCond->iBond[i][1]]);
  }
  
  printf("Let's check what there is in bondTypes\n");
  for(j=0;j<inp->nBondTypes;j++)
  {
     printf("%d %d %d\n",inp->nBondTypes,inp->bondTypes[j][0],inp->bondTypes[j][1]);
  }

  printf("Let's check what there is in iAngle\n");
  for(i=0;i<ff->nAngle;i++)
  {
     printf("%d %d %d %d\n",ff->nAngle,simulCond->iAngle[i][0],simulCond->iAngle[i][1],simulCond->iAngle[i][2]);
  }
  
  printf("Let's check what there is in angTypes\n");
  for(j=0;j<inp->nAngTypes;j++)
  {
     printf("%d %d %d %d\n",inp->nAngTypes,inp->angTypes[j][0],inp->angTypes[j][1],inp->angTypes[j][2]);
  }

  printf("Let's check what there is in iDihedral\n");
  for(i=0;i<ff->nDihedral;i++)
  {
     printf("%d %d %d %d %d\n",ff->nDihedral,simulCond->iDihedral[i][0],simulCond->iDihedral[i][1],simulCond->iDihedral[i][2],simulCond->iDihedral[i][3]);
  }
  
  printf("Let's check what there is in diheTypes\n");
  for(j=0;j<inp->nDiheTypes;j++)
  {
     printf("%d %d %d %d %d\n",inp->nDiheTypes,inp->diheTypes[j][0],inp->diheTypes[j][1],inp->diheTypes[j][2],inp->diheTypes[j][3]);
  }
  
  if(simulCond->keyconsth)
  {
    printf("Let's check what there is in constList\n");
    printf("pointer adress=%p\n",constList);
    for(i=1;i<simulCond->nconst;i++)
    {
      printf("%d %d %d\n",simulCond->nconst,constList[i].a,constList[i].b);
    }
  }

  printf("Now if the loop.\n");

  for(i=0;i<ff->nBond;i++)
  {
    ia=atom->atomType[simulCond->iBond[i][0]];
    ib=atom->atomType[simulCond->iBond[i][1]];
    
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
      ia=atom->atomType[constList[i].a];
      ib=atom->atomType[constList[i].b];
      
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
    ia=atom->atomType[simulCond->iAngle[i][0]];
    ib=atom->atomType[simulCond->iAngle[i][1]];
    ic=atom->atomType[simulCond->iAngle[i][2]];
    
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
    ia=atom->atomType[simulCond->iAngle[i][0]];
    ib=atom->atomType[simulCond->iAngle[i][1]];
    ic=atom->atomType[simulCond->iAngle[i][2]];
    
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
    ia=atom->atomType[simulCond->iDihedral[i][0]];
    ib=atom->atomType[simulCond->iDihedral[i][1]];
    ic=atom->atomType[simulCond->iDihedral[i][2]];
    id=atom->atomType[simulCond->iDihedral[i][3]];
    
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
    ia=atom->atomType[simulCond->iImproper[i][0]];
    ib=atom->atomType[simulCond->iImproper[i][1]];
    ic=atom->atomType[simulCond->iImproper[i][2]];
    id=atom->atomType[simulCond->iImproper[i][3]];
    
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
  
  for(i=0;i<atom->natom;i++)
  {
    
    ff->parmVdw[i][0]=sqrt(-kcaltoiu*inp->nonBondedTypesParm[atom->atomType[i]][1]);
    ff->parmVdw[i][1]=inp->nonBondedTypesParm[atom->atomType[i]][2]/sq6rt2;
    ff->parmVdw[i][2]=inp->nonBondedTypesParm[atom->atomType[i]][0];
    
    if(inp->nonBondedTypesParm[atom->atomType[i]][5]>=0.)
    {
      ff->parmVdw[i][3]=sqrt(-kcaltoiu*inp->nonBondedTypesParm[atom->atomType[i]][4]);
      ff->parmVdw[i][4]=inp->nonBondedTypesParm[atom->atomType[i]][5]/sq6rt2;
      ff->parmVdw[i][5]=inp->nonBondedTypesParm[atom->atomType[i]][3];
    }
    else
    {
      ff->parmVdw[i][3]=ff->parmVdw[i][0];
      ff->parmVdw[i][4]=ff->parmVdw[i][1];
      ff->parmVdw[i][5]=ff->parmVdw[i][2];
    }
  }
  
}

void write_FORF(INPUTS *inp,ATOM *atom,FORCEFIELD *ff,SIMULPARAMS *simulCond)
{
  FILE *forfFile;
  int i,j,k,l,ia,ib,ic,id,nd;
  double kf,r0,ncos;
  
  forfFile=fopen("FORF","w");
  
  fprintf(forfFile,"ATOMS %d\n",atom->natom);
  
  for(i=0;i<atom->natom;i++)
  {
    k=atom->atomType[i];
    fprintf(forfFile,"%s %lf %lf 1\n",inp->types[k],atom->m[i],ff->q[i]);
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
  
  nd=atom->natom*(atom->natom+1)/2;
  fprintf(forfFile,"VDW %d\n",nd);
  
  for(i=0;i<atom->natom;i++)
  {
    for(j=i;j<atom->natom;j++)
    {
      k=atom->atomType[i];
      l=atom->atomType[j];
      kf=ff->parmVdw[i][0]*ff->parmVdw[j][0]/kcaltoiu;
      r0=ff->parmVdw[i][1]+ff->parmVdw[j][1];
      
      fprintf(forfFile,"%s %s %lf %lf 1\n",inp->types[k],inp->types[l],kf,r0);
    }
  }
  
  fclose(forfFile);
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
  printf("MDBas failed due to error number: %d\n",errorNumber);
  switch (errorNumber)
  {
  case 10:
    printf("MDBas cannot find or open topology file TOP.\n");
    printf("Most likely, it is not properly named. Please check.\n");
    break;
  case 20:
    printf("MDBas cannot find or open structure file PSF.\n");
    printf("Most likely, it is not properly named. Please check.\n");
    break;
  case 22:
    printf("MDBas encountered a problem while reading atomic properties\n");
    printf("in the PSF file. There is an unexpected line there. Please\n");
    printf("consult the manual for further details about PSF file\n");
    break;
  case 23:
    printf("There is problem in bonds sequence in the PSF file. Please\n");
    printf("consult the manual for further details about PSF file\n");
    break;
  case 24:
    printf("There is problem in angles sequence in the PSF file. Please\n");
    printf("consult the manual for further details about PSF file\n");
    break;
  case 25:
    printf("There is problem in dihedrals sequence in the PSF file. Please\n");
    printf("consult the manual for further details about PSF file\n");
    break;
  case 26:
    printf("There is problem in improper angles sequence in the PSF file.\n");
    printf("Please consult the manual for further details about PSF file\n");
    break;
  case 30:
    printf("MDBas cannot find or open parameter file PAR.\n");
    printf("Most likely, it is not properly named. Please check.\n");
    break;
  case 40:
    printf("MDBas cannot find or open configuration file CONF.\n");
    printf("Most likely, it is not properly named. Please check.\n");
    break;
  case 41:
    printf("MDBas found a different number of atoms in CONF file and\n");
    printf("in PSF file. Structure does not match configuration.\n");
    printf("Check carefully these files.\n");
    break;
  case 50:
    printf("A dihedral angle is specified as a Fourier series but\n");
    printf("with one of the component being an harmonic potential.\n");
    printf("Check in PAR file.\n");
    break;
  case 60:
    printf("MDBas cannot find or open simulation file SIMU.\n");
    printf("Most likely, it is not properly named. Please check.\n");
    break;
  case 61:
    printf("MDBas does not recognise a keyword specified in SIMU.\n");
    printf("Please check SIMU file and the manual for the list of\n");
    printf("allowed keywords.\n");
    break;
  case 62:
    printf("MDBas does not recognise a parameter specified in SIMU.\n");
    printf("Please check SIMU file and the manual for the list of\n");
    printf("allowed keywords and their associated parameters.\n");
    break;
  case 63:
    printf("MDBas does not find a required parameter in SIMU.\n");
    printf("Please check SIMU file and the manual for the list of\n");
    printf("allowed keywords and their associated parameters.\n");
    break;
  case 71:
    printf("There is an undefined bond in the PSF. Most likely,\n");
    printf("there are missing parameters in the PAR file. Please check\n");
    break;
  case 72:
    printf("There is an undefined angle in the PSF. Most likely,\n");
    printf("there are missing parameters in the PAR file. Please check\n");
    break;
  case 73:
    printf("There is an undefined dihedral angle in the PSF. Most likely,\n");
    printf("there are missing parameters in the PAR file. Please check\n");
    break;
  case 74:
    printf("There is an undefined improper angle in the PSF. Most likely,\n");
    printf("there are missing parameters in the PAR file. Please check\n");
    break;
  case 110:
    printf("MDBas found a too many non-parameterised dihedral angles:\n");
    printf("4*nDihedrals. nDihedrals comes from the value specified\n");
    printf("in PSF file. Please check in PAR file. If such a number is\n");
    printf("normal for your simulation, you have to enter list.c to\n");
    printf("increase the size of the 1-4 pairs array from 5*nDihedrals to\n");
    printf("the size you really need. Then recompile MDBas.\n");
    break;
  case 111:
    printf("MDBas encountered a problem while setting the excluded atoms\n");
    printf("list. The last atom has exclusion which should not happen. This\n");
    printf("a bit annoying for there is no simple explanation for this.\n");
    printf("Maybe an error in one of the input files which is not detected\n");
    printf("by MDBas. Sorry for the trouble.\n");
    break;
  case 112:
    printf("MDBas encountered a problem while setting the excluded atoms\n");
    printf("list. The total excluded atoms does not match of the sum of\n");
    printf("excluded atoms for each atom. This a bit annoying for there is\n");
    printf("no simple explanation for this. Maybe an error in one of the\n");
    printf("input files which is not detected by MDBas. Sorry for the trouble.\n");
    break;
  case 201:
    printf("Unknown electrostatic potential. This is most likely due to an\n");
    printf("error in the SIMU file. Please check this file and the manual\n");
    printf("for the list of keywords and available potentials.\n");
    break;
  case 202:
    printf("Unknown van der Waals potential. This is most likely due to an\n");
    printf("error in the SIMU file. Please check this file and the manual\n");
    printf("for the list of keywords and available potentials.\n");
    break;
  case 311:
    printf("Shake convergence failure, most likely due a non suitable initial\n");
    printf("configuration. If not, you can try increasing the number of cycles or\n");
    printf("make Shake convergence criterion more tolerant. Please check the manual.\n");
    break;
  default:
    printf("MDBas failed due to unknown error number: %d\n",errorNumber);
    printf("Reading the manual will not help you. You are by yourself.\n");
    printf("Errare humanum est.\n");
    break;
  }
  
  exit(errorNumber);
    
}
