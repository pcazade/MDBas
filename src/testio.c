#include <stdio.h>
#include <stdlib.h>
#include "global.h"
#include "io.h"

void read_PSF2(INPUTS *inp,ATOM *atom,FORCEFIELD *ff,ENERGYFORCE *enerFor,SIMULPARAMS *simulCond)
{
  FILE *psfFile=NULL;
  char buff1[1024]="", *buff2=NULL, *buff3=NULL, *buff4=NULL, *buff5=NULL;
  int i,j,k,kk,kt,nalloc=100;nincr=100;
  
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
  inp->typesNum=(int*)realloc(inp->typesNum,*sizeof(inp->nTypes*(inp->typesNum)));
  
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
  
  i=0;
  while(fgets(buff1,1024,psfFile)!=NULL)
  {
    buff2=strtok(buff1," \n\t");
    buff3=strtok(NULL," \n\t");
    
    while(buff2!=NULL && buff3!=NULL)
    {
      
      simulCond->iBond[i][0]=atoi(buff2)-1;
      simulCond->iBond[i][1]=atoi(buff3)-1;
      
      buff2=strtok(NULL," \n\t");
      buff3=strtok(NULL," \n\t");
      i++;
    }
    if(buff2!=NULL && buff3==NULL)
      error(23);
    if(i==ff->nBond)
      break;
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

void read_TOP2(INPUTS *inp)
{
  
  FILE *topFile=NULL;
  char buff1[1024]="", *buff2=NULL, *buff3=NULL, *buff4=NULL ;
  int i,j,k;
  
  topFile=fopen("TOP","r");

  if (topFile==NULL)
  {
    error(10);
  }
  
  inp->types=(char**)malloc(inp->nTypes*sizeof(*(inp->types)));
  
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


void read_PAR2(INPUTS *inp)
{
  
  char buff1[1024]="", *buff2=NULL, *buff3=NULL, *buff4=NULL, *buff5=NULL, *buff6=NULL;

    FILE *parFile=NULL;
    int i,j,l,n,k=0;
    
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
    
    int ia,ib,ic,id,ii,jj,kk;
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
	    break;
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
	    break;
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
	  
	  if(!strcmp(buff2,"X"))
	  {
	    ia=inp->nTypes;
	    break;
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
	    break;
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
	    break;
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
	    break;
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
	  inp->nImproper=i;
	}
      }
      else if(k==5)      
      {
	if((buff2!=NULL)&&(buff2[0]!='!'))
	{
	  buff3=strtok(NULL," \n\t");
	  buff4=strtok(NULL," \n\t");
	  buff5=strtok(NULL," \n\t");
	  
	  ii=-1
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

