#include <stdlib.h>
#include "global.h"
#include "utils.h"
#include "memory.h"
#include "io.h"
#include "list.h"

void makelist(SIMULPARAMS *simulCond,ATOM *atom,FORCEFIELD *ff,CONSTRAINT *constList)
{
  if(simulCond->firstener==1)
  {
//     printf("Building exclusion and Verlet lists.\n");
    exclude_list(simulCond,atom,ff,constList);
    verlet_list(simulCond,atom,ff);
    
    simulCond->firstener=0;
    
  }
  else if( (simulCond->step>0) && (simulCond->step%simulCond->listupdate==0 ) )
  {
//     printf("Updating Verlet list.\n");
    verlet_list_update(simulCond,atom,ff);
  }
}

void exclude_list(SIMULPARAMS *simulCond,ATOM *atom,FORCEFIELD *ff,CONSTRAINT *constList)
{
  int i,j,k,l,ii,jj,kk,ia,ib,ic,id,exclude;
  int **tempAtom,**tempVer14,**tempConnect,*tempConnectNum;
  int nAlloc=12,nIncr=10,nConnect=12,nExclude;
  
  simulCond->excludeNum=(int*)malloc(atom->natom*sizeof(*(simulCond->excludeNum)));
  for(i=0;i<atom->natom;i++)
    simulCond->excludeNum[i]=0;
    
  tempAtom=(int**)malloc(atom->natom*sizeof(*(tempAtom)));
  for(i=0;i<atom->natom;i++)
    tempAtom[i]=(int*)malloc(nAlloc*sizeof(**(tempAtom)));
  
/* The following temporary arrays are designed for chechking the
 * connectivity of system and make sure that all existing but not
 * parameterised angles and dihedrals are properly excluded. */
  
  tempConnectNum=(int*)malloc(atom->natom*sizeof(*tempConnectNum));
  for(i=0;i<atom->natom;i++)
    tempConnectNum[i]=0;
  
  tempConnect=(int**)malloc(atom->natom*sizeof(*tempConnect));
  for(i=0;i<atom->natom;i++)
    tempConnect[i]=(int*)malloc(nConnect*sizeof(**tempConnect));
  
    
  for(i=0;i<ff->nBond;i++)
  {
    ia=simulCond->iBond[i][0];
    ib=simulCond->iBond[i][1];
    
    if(ib<ia)
    {
      ii=ib;
      jj=ia;
    }
    else
    {
      ii=ia;
      jj=ib;
    }
  
    if(simulCond->excludeNum[ii]>=nAlloc)
    {
      nAlloc+=nIncr;
      for(j=0;j<atom->natom;j++)
	tempAtom[j]=(int*)realloc(tempAtom[j],nAlloc*sizeof(**(tempAtom)));
    }
    
    if(tempConnectNum[ii]>=nConnect||tempConnectNum[jj]>=nConnect)
    {
      nConnect+=nIncr;
      for(j=0;j<atom->natom;j++)
	tempConnect[j]=(int*)realloc(tempConnect[j],nConnect*sizeof(**tempConnect));
    }
    
    tempConnect[ii][tempConnectNum[ii]]=jj;
    tempConnect[jj][tempConnectNum[jj]]=ii;
    tempConnectNum[ii]++;
    tempConnectNum[jj]++;
    
    tempAtom[ii][simulCond->excludeNum[ii]]=jj;
    simulCond->excludeNum[ii]++;      
  }
  
  if(simulCond->keyconsth)
  {
    for(i=0;i<simulCond->nconst;i++)
    {
      ia=constList[i].a;
      ib=constList[i].b;
      
      if(ib<ia)
      {
	ii=ib;
	jj=ia;
      }
      else
      {
	ii=ia;
	jj=ib;
      }
    
      if(simulCond->excludeNum[ii]>=nAlloc)
      {
	nAlloc+=nIncr;
	for(j=0;j<atom->natom;j++)
	  tempAtom[j]=(int*)realloc(tempAtom[j],nAlloc*sizeof(**(tempAtom)));
      }
      
      if(tempConnectNum[ii]>=nConnect||tempConnectNum[jj]>=nConnect)
      {
	nConnect+=nIncr;
	for(j=0;j<atom->natom;j++)
	  tempConnect[j]=(int*)realloc(tempConnect[j],nConnect*sizeof(**tempConnect));
      }
      
      tempConnect[ii][tempConnectNum[ii]]=jj;
      tempConnect[jj][tempConnectNum[jj]]=ii;
      tempConnectNum[ii]++;
      tempConnectNum[jj]++;
      
      tempAtom[ii][simulCond->excludeNum[ii]]=jj;
      simulCond->excludeNum[ii]++;      
    }
  }
  
  for(i=0;i<ff->nAngle;i++)
  {
    ia=simulCond->iAngle[i][0];
    ib=simulCond->iAngle[i][1];
    ic=simulCond->iAngle[i][2];
    
    if(ib<ia)
    {
      ii=ib;
      jj=ia;
    }
    else
    {
      ii=ia;
      jj=ib;
    }
  
    if(simulCond->excludeNum[ii]>=nAlloc)
    {
      nAlloc+=nIncr;
      for(j=0;j<atom->natom;j++)
	tempAtom[j]=(int*)realloc(tempAtom[j],nAlloc*sizeof(**(tempAtom)));
    }
    
    exclude=1;
    for(j=0;j<simulCond->excludeNum[ii];j++)
    {
      if(tempAtom[ii][j]==jj)
      {
	exclude=0;
	break;
      }
    }
    
    if(exclude)
    {
      tempAtom[ii][simulCond->excludeNum[ii]]=jj;
      simulCond->excludeNum[ii]++;
    }
    
    if(ic<ia)
    {
      ii=ic;
      jj=ia;
    }
    else
    {
      ii=ia;
      jj=ic;
    }
  
    if(simulCond->excludeNum[ii]>=nAlloc)
    {
      nAlloc+=nIncr;
      for(j=0;j<atom->natom;j++)
	tempAtom[j]=(int*)realloc(tempAtom[j],nAlloc*sizeof(**(tempAtom)));
    }
    
    exclude=1;
    for(j=0;j<simulCond->excludeNum[ii];j++)
    {
      if(tempAtom[ii][j]==jj)
      {
	exclude=0;
	break;
      }
    }
    
    if(exclude)
    {
      tempAtom[ii][simulCond->excludeNum[ii]]=jj;
      simulCond->excludeNum[ii]++;
    }
    
    if(ic<ib)
    {
      ii=ic;
      jj=ib;
    }
    else
    {
      ii=ib;
      jj=ic;
    }
  
    if(simulCond->excludeNum[ii]>=nAlloc)
    {
      nAlloc+=nIncr;
      for(j=0;j<atom->natom;j++)
	tempAtom[j]=(int*)realloc(tempAtom[j],nAlloc*sizeof(**(tempAtom)));
    }
    
    exclude=1;
    for(j=0;j<simulCond->excludeNum[ii];j++)
    {
      if(tempAtom[ii][j]==jj)
      {
	exclude=0;
	break;
      }
    }
    
    if(exclude)
    {
      tempAtom[ii][simulCond->excludeNum[ii]]=jj;
      simulCond->excludeNum[ii]++;
    }
       
  }
  
  tempVer14=(int**)malloc(5*ff->nDihedral*sizeof(*tempVer14));
  for(i=0;i<5*ff->nDihedral;i++)
    tempVer14[i]=(int*)malloc(2*sizeof(**tempVer14));
  
  ff->npr14=0;
  
  for(i=0;i<ff->nDihedral;i++)
  {
    ia=simulCond->iDihedral[i][0];
    ib=simulCond->iDihedral[i][1];
    ic=simulCond->iDihedral[i][2];
    id=simulCond->iDihedral[i][3];
    
    if(ib<ia)
    {
      ii=ib;
      jj=ia;
    }
    else
    {
      ii=ia;
      jj=ib;
    }
  
    if(simulCond->excludeNum[ii]>=nAlloc)
    {
      nAlloc+=nIncr;
      for(j=0;j<atom->natom;j++)
	tempAtom[j]=(int*)realloc(tempAtom[j],nAlloc*sizeof(**(tempAtom)));
    }
    
    exclude=1;
    for(j=0;j<simulCond->excludeNum[ii];j++)
    {
      if(tempAtom[ii][j]==jj)
      {
	exclude=0;
	break;
      }
    }
    
    if(exclude)
    {
      tempAtom[ii][simulCond->excludeNum[ii]]=jj;
      simulCond->excludeNum[ii]++;
    }
    
    if(ic<ia)
    {
      ii=ic;
      jj=ia;
    }
    else
    {
      ii=ia;
      jj=ic;
    }
  
    if(simulCond->excludeNum[ii]>=nAlloc)
    {
      nAlloc+=nIncr;
      for(j=0;j<atom->natom;j++)
	tempAtom[j]=(int*)realloc(tempAtom[j],nAlloc*sizeof(**(tempAtom)));
    }
    
    exclude=1;
    for(j=0;j<simulCond->excludeNum[ii];j++)
    {
      if(tempAtom[ii][j]==jj)
      {
	exclude=0;
	break;
      }
    }
    
    if(exclude)
    {
      tempAtom[ii][simulCond->excludeNum[ii]]=jj;
      simulCond->excludeNum[ii]++;
    }
    
    if(id<ia)
    {
      ii=id;
      jj=ia;
    }
    else
    {
      ii=ia;
      jj=id;
    }
  
    if(simulCond->excludeNum[ii]>=nAlloc)
    {
      nAlloc+=nIncr;
      for(j=0;j<atom->natom;j++)
	tempAtom[j]=(int*)realloc(tempAtom[j],nAlloc*sizeof(**(tempAtom)));
    }
    
    exclude=1;
    for(j=0;j<simulCond->excludeNum[ii];j++)
    {
      if(tempAtom[ii][j]==jj)
      {
	exclude=0;
	break;
      }
    }
    
    if(exclude)
    {
      tempVer14[ff->npr14][0]=ia;
      tempVer14[ff->npr14][1]=id;
      ff->npr14++;
      
      tempAtom[ii][simulCond->excludeNum[ii]]=jj;
      simulCond->excludeNum[ii]++;
    }
    
    if(ic<ib)
    {
      ii=ic;
      jj=ib;
    }
    else
    {
      ii=ib;
      jj=ic;
    }
  
    if(simulCond->excludeNum[ii]>=nAlloc)
    {
      nAlloc+=nIncr;
      for(j=0;j<atom->natom;j++)
	tempAtom[j]=(int*)realloc(tempAtom[j],nAlloc*sizeof(**(tempAtom)));
    }
    
    exclude=1;
    for(j=0;j<simulCond->excludeNum[ii];j++)
    {
      if(tempAtom[ii][j]==jj)
      {
	exclude=0;
	break;
      }
    }
    
    if(exclude)
    {
      tempAtom[ii][simulCond->excludeNum[ii]]=jj;
      simulCond->excludeNum[ii]++;
    }
    
    if(id<ib)
    {
      ii=id;
      jj=ib;
    }
    else
    {
      ii=ib;
      jj=id;
    }
  
    if(simulCond->excludeNum[ii]>=nAlloc)
    {
      nAlloc+=nIncr;
      for(j=0;j<atom->natom;j++)
	tempAtom[j]=(int*)realloc(tempAtom[j],nAlloc*sizeof(**(tempAtom)));
    }
    
    exclude=1;
    for(j=0;j<simulCond->excludeNum[ii];j++)
    {
      if(tempAtom[ii][j]==jj)
      {
	exclude=0;
	break;
      }
    }
    
    if(exclude)
    {
      tempAtom[ii][simulCond->excludeNum[ii]]=jj;
      simulCond->excludeNum[ii]++;
    }
    
    if(id<ic)
    {
      ii=id;
      jj=ic;
    }
    else
    {
      ii=ic;
      jj=id;
    }
  
    if(simulCond->excludeNum[ii]>=nAlloc)
    {
      nAlloc+=nIncr;
      for(j=0;j<atom->natom;j++)
	tempAtom[j]=(int*)realloc(tempAtom[j],nAlloc*sizeof(**(tempAtom)));
    }
    
    exclude=1;
    for(j=0;j<simulCond->excludeNum[ii];j++)
    {
      if(tempAtom[ii][j]==jj)
      {
	exclude=0;
	break;
      }
    }
    
    if(exclude)
    {
      tempAtom[ii][simulCond->excludeNum[ii]]=jj;
      simulCond->excludeNum[ii]++;
    }
       
  }  
    
  for(i=0;i<ff->nImproper;i++)
  {
    ia=simulCond->iImproper[i][0];
    ib=simulCond->iImproper[i][1];
    ic=simulCond->iImproper[i][2];
    id=simulCond->iImproper[i][3];
    
    if(ib<ia)
    {
      ii=ib;
      jj=ia;
    }
    else
    {
      ii=ia;
      jj=ib;
    }
  
    if(simulCond->excludeNum[ii]>=nAlloc)
    {
      nAlloc+=nIncr;
      for(j=0;j<atom->natom;j++)
	tempAtom[j]=(int*)realloc(tempAtom[j],nAlloc*sizeof(**(tempAtom)));
    }
    
    exclude=1;
    for(j=0;j<simulCond->excludeNum[ii];j++)
    {
      if(tempAtom[ii][j]==jj)
      {
	exclude=0;
	break;
      }
    }
    
    if(exclude)
    {
      tempAtom[ii][simulCond->excludeNum[ii]]=jj;
      simulCond->excludeNum[ii]++;
    }
    
    if(ic<ia)
    {
      ii=ic;
      jj=ia;
    }
    else
    {
      ii=ia;
      jj=ic;
    }
  
    if(simulCond->excludeNum[ii]>=nAlloc)
    {
      nAlloc+=nIncr;
      for(j=0;j<atom->natom;j++)
	tempAtom[j]=(int*)realloc(tempAtom[j],nAlloc*sizeof(**(tempAtom)));
    }
    
    exclude=1;
    for(j=0;j<simulCond->excludeNum[ii];j++)
    {
      if(tempAtom[ii][j]==jj)
      {
	exclude=0;
	break;
      }
    }
    
    if(exclude)
    {
      tempAtom[ii][simulCond->excludeNum[ii]]=jj;
      simulCond->excludeNum[ii]++;
    }
    
    if(id<ia)
    {
      ii=id;
      jj=ia;
    }
    else
    {
      ii=ia;
      jj=id;
    }
  
    if(simulCond->excludeNum[ii]>=nAlloc)
    {
      nAlloc+=nIncr;
      for(j=0;j<atom->natom;j++)
	tempAtom[j]=(int*)realloc(tempAtom[j],nAlloc*sizeof(**(tempAtom)));
    }
    
    exclude=1;
    for(j=0;j<simulCond->excludeNum[ii];j++)
    {
      if(tempAtom[ii][j]==jj)
      {
	exclude=0;
	break;
      }
    }
    
    if(exclude)
    {
      tempAtom[ii][simulCond->excludeNum[ii]]=jj;
      simulCond->excludeNum[ii]++;
    }
    
    if(ic<ib)
    {
      ii=ic;
      jj=ib;
    }
    else
    {
      ii=ib;
      jj=ic;
    }
  
    if(simulCond->excludeNum[ii]>=nAlloc)
    {
      nAlloc+=nIncr;
      for(j=0;j<atom->natom;j++)
	tempAtom[j]=(int*)realloc(tempAtom[j],nAlloc*sizeof(**(tempAtom)));
    }
    
    exclude=1;
    for(j=0;j<simulCond->excludeNum[ii];j++)
    {
      if(tempAtom[ii][j]==jj)
      {
	exclude=0;
	break;
      }
    }
    
    if(exclude)
    {
      tempAtom[ii][simulCond->excludeNum[ii]]=jj;
      simulCond->excludeNum[ii]++;
    }
    
    if(id<ib)
    {
      ii=id;
      jj=ib;
    }
    else
    {
      ii=ib;
      jj=id;
    }
  
    if(simulCond->excludeNum[ii]>=nAlloc)
    {
      nAlloc+=nIncr;
      for(j=0;j<atom->natom;j++)
	tempAtom[j]=(int*)realloc(tempAtom[j],nAlloc*sizeof(**(tempAtom)));
    }
    
    exclude=1;
    for(j=0;j<simulCond->excludeNum[ii];j++)
    {
      if(tempAtom[ii][j]==jj)
      {
	exclude=0;
	break;
      }
    }
    
    if(exclude)
    {
      tempAtom[ii][simulCond->excludeNum[ii]]=jj;
      simulCond->excludeNum[ii]++;
    }
    
    if(id<ic)
    {
      ii=id;
      jj=ic;
    }
    else
    {
      ii=ic;
      jj=id;
    }
  
    if(simulCond->excludeNum[ii]>=nAlloc)
    {
      nAlloc+=nIncr;
      for(j=0;j<atom->natom;j++)
	tempAtom[j]=(int*)realloc(tempAtom[j],nAlloc*sizeof(**(tempAtom)));
    }
    
    exclude=1;
    for(j=0;j<simulCond->excludeNum[ii];j++)
    {
      if(tempAtom[ii][j]==jj)
      {
	exclude=0;
	break;
      }
    }
    
    if(exclude)
    {
      tempAtom[ii][simulCond->excludeNum[ii]]=jj;
      simulCond->excludeNum[ii]++;
    }
       
  }
  
  for(ia=0;ia<atom->natom;ia++)
  {
    for(j=0;j<tempConnectNum[ia];j++)
    {
      ib=tempConnect[ia][j];
      for(k=0;k<tempConnectNum[ib];k++)
      {
	ic=tempConnect[ib][k];
	
	if(ic==ia)
	  continue;
	
	if(ic<ia)
	{
	  ii=ic;
	  jj=ia;
	}
	else
	{
	  ii=ia;
	  jj=ic;
	}
      
	if(simulCond->excludeNum[ii]>=nAlloc)
	{
	  nAlloc+=nIncr;
	  for(kk=0;kk<atom->natom;kk++)
	    tempAtom[kk]=(int*)realloc(tempAtom[kk],nAlloc*sizeof(**(tempAtom)));
	}
	
	exclude=1;
	for(kk=0;kk<simulCond->excludeNum[ii];kk++)
	{
	  if(tempAtom[ii][kk]==jj)
	  {
	    exclude=0;
	    break;
	  }
	}
	
	if(exclude)
	{
	  tempAtom[ii][simulCond->excludeNum[ii]]=jj;
	  simulCond->excludeNum[ii]++;
	}
	
	for(l=0;l<tempConnectNum[ic];l++)
	{
	  id=tempConnect[ic][l];
	  
	  if(id==ia||id==ib)
	    continue;
	  
	  if(id<ia)
	  {
	    ii=id;
	    jj=ia;
	  }
	  else
	  {
	    ii=ia;
	    jj=id;
	  }
	
	  if(simulCond->excludeNum[ii]>=nAlloc)
	  {
	    nAlloc+=nIncr;
	    for(kk=0;kk<atom->natom;kk++)
	      tempAtom[kk]=(int*)realloc(tempAtom[kk],nAlloc*sizeof(**(tempAtom)));
	  }
	  
	  if(ff->npr14>=5*ff->nDihedral)
	    error(110);
	  
	  exclude=1;
	  for(kk=0;kk<simulCond->excludeNum[ii];kk++)
	  {
	    if(tempAtom[ii][kk]==jj)
	    {
	      exclude=0;
	      break;
	    }
	  }
	  
	  if(exclude)
	  {
	    
	    tempVer14[ff->npr14][0]=ia;
	    tempVer14[ff->npr14][1]=id;
	    ff->npr14++;
	    
	    tempAtom[ii][simulCond->excludeNum[ii]]=jj;
	    simulCond->excludeNum[ii]++;
	  }
	  
	}
      }
    }
  }
  
  ff->ver14=(int**)malloc(ff->npr14*sizeof(*(ff->ver14)));
  for(i=0;i<ff->npr14;i++)
    ff->ver14[i]=(int*)malloc(2*sizeof(**(ff->ver14)));
  
  for(i=0;i<ff->npr14;i++)
  {
    ff->ver14[i][0]=tempVer14[i][0];
    ff->ver14[i][1]=tempVer14[i][1];
  }
  
  free_2D(5*ff->nDihedral,tempVer14,NULL);
  
  free_2D(atom->natom,tempConnect,NULL);
  free(tempConnectNum);
  
    
  nExclude=0;
  for(i=0;i<atom->natom;i++)
    nExclude+=simulCond->excludeNum[i];
  
  if(simulCond->excludeNum[atom->natom-1]!=0)
    error(111);
  
  simulCond->excludeAtom=(int*)malloc(nExclude*sizeof(*(simulCond->excludeAtom)));
  
  k=0;
  for(i=0;i<atom->natom;i++)
  {
    for(j=0;j<simulCond->excludeNum[i];j++)
    {
      simulCond->excludeAtom[k]=tempAtom[i][j];
      k++;
    }
  }
  
  if(k!=nExclude)
    error(112);
  
  free_2D(atom->natom,tempAtom,NULL);

}

void verlet_list(SIMULPARAMS *simulCond,ATOM *atom,FORCEFIELD *ff)
{
  int i,j,k,kv,exclude,nalloc,incr=1000;
  double r,cutnb,delta[3];
  
  cutnb=simulCond->cutoff+simulCond->delr;
  
  ff->verPair=(int*)malloc((atom->natom-1)*sizeof(*(ff->verPair)));
  for(i=0;i<atom->natom-1;i++)
    ff->verPair[i]=0;
  
  nalloc=((atom->natom*atom->natom/64)-(atom->natom/8))/2;
  ff->verList=(int*)malloc(nalloc*sizeof(*(ff->verList)));
  for(i=0;i<nalloc;i++)
    ff->verList[i]=0;
  
  kv=0;
  ff->npr=0;
  for(i=0;i<atom->natom-1;i++)
  {
    for(j=i+1;j<atom->natom;j++)
    {
      r=distance(i,j,atom,delta,simulCond);
      if(r<=cutnb)
      {
	exclude=0;
	for (k=0;k<simulCond->excludeNum[i];k++)
	{
	  if(simulCond->excludeAtom[kv+k]==j)
	  {
	    exclude=1;
	    break;
	  }
	}
	if(!exclude)
	{
	  if(ff->npr>=nalloc)
	  {
	    nalloc+=incr;
	    ff->verList=(int*)realloc(ff->verList,nalloc*sizeof(*(ff->verList)));
	  }
	  ff->verList[ff->npr]=j;
	  ff->verPair[i]++;
	  ff->npr++;
	}
      }
    }
    kv+=simulCond->excludeNum[i];
  }
  ff->verList=(int*)realloc(ff->verList,ff->npr*sizeof(*(ff->verList)));
  
}

void verlet_list_update(SIMULPARAMS *simulCond,ATOM *atom,FORCEFIELD *ff)
{
  int i,j,k,kv,exclude,nalloc,incr=1000;
  double r,cutnb,delta[3];
  
  cutnb=simulCond->cutoff+simulCond->delr;
  nalloc=ff->npr;
    
  kv=0;
  ff->npr=0;
  for(i=0;i<atom->natom-1;i++)
  {
    ff->verPair[i]=0;
    for(j=i+1;j<atom->natom;j++)
    {
      r=distance(i,j,atom,delta,simulCond);
      if(r<=cutnb)
      {
	exclude=0;
	for (k=0;k<simulCond->excludeNum[i];k++)
	{
	  if(simulCond->excludeAtom[kv+k]==j)
	  {
	    exclude=1;
	    break;
	  }
	}
	if(!exclude)
	{
	  if(ff->npr>=nalloc)
	  {
	    nalloc+=incr;
	    ff->verList=(int*)realloc(ff->verList,nalloc*sizeof(*(ff->verList)));
	  }
	  ff->verList[ff->npr]=j;
	  ff->verPair[i]++;
	  ff->npr++;
	}
      }
    }
    kv+=simulCond->excludeNum[i];
  }
  ff->verList=(int*)realloc(ff->verList,ff->npr*sizeof(*(ff->verList)));
}
