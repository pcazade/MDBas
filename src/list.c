#include <stdlib.h>
#include "global.h"
#include "utils.h"
#include "memory.h"
#include "io.h"
#include "list.h"

void makelist(SIMULPARAMS *simulCond,ATOM *atom,FORCEFIELD *ff,CONSTRAINT *constList,PBC *box)
{
  if(simulCond->firstener==1)
  {
    //printf("Building exclusion and Verlet lists.\n");
    exclude_list(simulCond,atom,ff,constList);
    //printf("Exclude list done\n");
    verlet_list(simulCond,atom,ff,box);
    //printf("Verlet list done\n");
    
    simulCond->firstener=0;
    
  }
  else if( /*(simulCond->step>0) &&*/ (simulCond->step%simulCond->listupdate==0 ) )
  {
    //printf("Updating Verlet list at step %d.\n",simulCond->step);
    verlet_list_update(simulCond,atom,ff,box);
    //printf("Update done\n\n");
  }
}

void exclude_list(SIMULPARAMS *simulCond,ATOM *atom,FORCEFIELD *ff,CONSTRAINT *constList)
{
  int i,j,k,l,ii,jj,kk,ia,ib,ic,id,exclude;
  int **tempAtom=NULL,**tempVer14=NULL,**tempConnect=NULL,*tempConnectNum=NULL;
  int nAlloc=12,nIncr=10,nConnect=12,nExclude;
  
  simulCond->excludeNum=(int*)malloc(simulCond->natom*sizeof(*(simulCond->excludeNum)));
  for(i=0;i<simulCond->natom;i++)
    simulCond->excludeNum[i]=0;
    
  tempAtom=(int**)malloc(simulCond->natom*sizeof(*(tempAtom)));
  for(i=0;i<simulCond->natom;i++)
    tempAtom[i]=(int*)malloc(nAlloc*sizeof(**(tempAtom)));
  
/* The following temporary arrays are designed for chechking the
 * connectivity of system and make sure that all existing but not
 * parameterised angles and dihedrals are properly excluded. */
  
  tempConnectNum=(int*)malloc(simulCond->natom*sizeof(*tempConnectNum));
  for(i=0;i<simulCond->natom;i++)
    tempConnectNum[i]=0;
  
  tempConnect=(int**)malloc(simulCond->natom*sizeof(*tempConnect));
  for(i=0;i<simulCond->natom;i++)
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
      for(j=0;j<simulCond->natom;j++)
	tempAtom[j]=(int*)realloc(tempAtom[j],nAlloc*sizeof(**(tempAtom)));
    }
    
    if(tempConnectNum[ii]>=nConnect||tempConnectNum[jj]>=nConnect)
    {
      nConnect+=nIncr;
      for(j=0;j<simulCond->natom;j++)
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
	for(j=0;j<simulCond->natom;j++)
	  tempAtom[j]=(int*)realloc(tempAtom[j],nAlloc*sizeof(**(tempAtom)));
      }
      
      if(tempConnectNum[ii]>=nConnect||tempConnectNum[jj]>=nConnect)
      {
	nConnect+=nIncr;
	for(j=0;j<simulCond->natom;j++)
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
      for(j=0;j<simulCond->natom;j++)
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
      for(j=0;j<simulCond->natom;j++)
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
      for(j=0;j<simulCond->natom;j++)
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
      for(j=0;j<simulCond->natom;j++)
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
      for(j=0;j<simulCond->natom;j++)
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
      for(j=0;j<simulCond->natom;j++)
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
      for(j=0;j<simulCond->natom;j++)
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
      for(j=0;j<simulCond->natom;j++)
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
      for(j=0;j<simulCond->natom;j++)
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
      for(j=0;j<simulCond->natom;j++)
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
      for(j=0;j<simulCond->natom;j++)
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
      for(j=0;j<simulCond->natom;j++)
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
      for(j=0;j<simulCond->natom;j++)
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
      for(j=0;j<simulCond->natom;j++)
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
      for(j=0;j<simulCond->natom;j++)
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
  
  for(ia=0;ia<simulCond->natom;ia++)
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
	  for(kk=0;kk<simulCond->natom;kk++)
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
	    for(kk=0;kk<simulCond->natom;kk++)
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
  
  free_2D(simulCond->natom,tempConnect,NULL);
  free(tempConnectNum);
  
    
  nExclude=0;
  for(i=0;i<simulCond->natom;i++)
    nExclude+=simulCond->excludeNum[i];
  
  if(simulCond->excludeNum[simulCond->natom-1]!=0)
    error(111);
  
  simulCond->excludeAtom=(int*)malloc(nExclude*sizeof(*(simulCond->excludeAtom)));
  
  k=0;
  for(i=0;i<simulCond->natom;i++)
  {
    for(j=0;j<simulCond->excludeNum[i];j++)
    {
      simulCond->excludeAtom[k]=tempAtom[i][j];
      k++;
    }
  }
  
  if(k!=nExclude)
    error(112);
  
  free_2D(simulCond->natom,tempAtom,NULL);

}

void verlet_list(SIMULPARAMS *simulCond,ATOM *atom,FORCEFIELD *ff,PBC *box)
{
  int i,j,k,kv,exclude,nalloc,incr=262144;
  double r,cutnb,delta[3];
  
  cutnb=simulCond->cutoff+simulCond->delr;
  
  ff->verPair=(int*)malloc((simulCond->natom-1)*sizeof(*(ff->verPair)));
  for(i=0;i<simulCond->natom-1;i++)
    ff->verPair[i]=0;

//   nalloc=((simulCond->natom*simulCond->natom/64)-(simulCond->natom/8))/2;
//   nalloc = (int) simulCond->natom * ( simulCond->natom/8 - 1  ) / 2 ;
  nalloc = (int) simulCond->natom * 1024 ;
  ff->verList=(int*)malloc(nalloc*sizeof(*(ff->verList)));
//   for(i=0;i<nalloc;i++)
//     ff->verList[i]=0;
  
  kv=0;
  ff->npr=0;
  for(i=0;i<simulCond->natom-1;i++)
  {
    for(j=i+1;j<simulCond->natom;j++)
    {
      r=distance(i,j,atom,delta,simulCond,box);
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
    
//     printf("Verpair[%d] : %d\n",i,ff->verPair[i]);
  }
  ff->verList=(int*)realloc(ff->verList,ff->npr*sizeof(*(ff->verList)));

#ifdef _OPENMP
  /** Verlet cumulatedSum list (useful for parallelisation)**/
  ff->verCumSum=(int*)malloc((simulCond->natom-1)*sizeof(*(ff->verCumSum)));
  ff->verCumSum[0] = 0;
  for(i=1;i<simulCond->natom-1;i++)
    ff->verCumSum[i] = ff->verCumSum[i-1] + ff->verPair[i-1];
  /** **/
#endif

}

void verlet_list_update(SIMULPARAMS *simulCond,ATOM *atom,FORCEFIELD *ff,PBC *box)
{
  int i,j,k,kv,exclude,nalloc,incr=1000;
  double r,cutnb,delta[3];
  
  cutnb=simulCond->cutoff+simulCond->delr;
  nalloc=ff->npr;
    
  kv=0;
  ff->npr=0;
  for(i=0;i<simulCond->natom-1;i++)
  {
    ff->verPair[i]=0;
    for(j=i+1;j<simulCond->natom;j++)
    {
      r=distance(i,j,atom,delta,simulCond,box);
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

#ifdef _OPENMP
  /** Verlet cumulatedSum list (useful for parallelisation)**/
  ff->verCumSum[0] = 0;
  for(i=1;i<simulCond->natom-1;i++)
    ff->verCumSum[i] = ff->verCumSum[i-1] + ff->verPair[i-1];
  /** **/
#endif

}

/*void link_cell_verlet_list(SIMULPARAMS *simulCond,ATOM *atom,FORCEFIELD *ff)
{
  
  int transx[508]={0,1,0,0,-1,1,0,-1,1,0,-1,1,-1,1,2,0,0,-2,2,-1,1,0,-2,2,0,
    0,-1,1,0,-1,1,-2,2,-2,2,-1,1,-1,1,-1,1,-2,2,0,-2,2,0,-2,2,-2,2,
    -1,1,-2,2,-2,2,-1,1,-2,2,-2,2,3,0,0,-3,3,-1,1,0,-3,3,0,0,-1,1,0,
    -1,1,-3,3,-3,3,-1,1,-1,1,-1,1,-3,3,-2,2,0,-3,3,0,0,-2,2,0,-2,2,
    -3,3,-3,3,-2,2,-1,1,-3,3,-3,3,-1,1,-1,1,-2,2,-2,2,-1,1,-2,2,-3,3,
    -3,3,-2,2,-2,2,-2,2,-3,3,0,-3,3,0,-3,3,-3,3,-1,1,-3,3,-3,3,-1,1,
    -3,3,-3,3,-2,2,-3,3,-3,3,-2,2,-3,3,-3,3,4,0,0,-4,4,-1,1,0,-4,4,0,
    0,-1,1,0,-1,1,-4,4,-4,4,-1,1,-1,1,-1,1,-4,4,-2,2,0,-4,4,0,0,-2,2,
    0,-2,2,-4,4,-4,4,-2,2,-1,1,-4,4,-4,4,-1,1,-1,1,-2,2,-2,2,-1,1,-2,
    2,-4,4,-4,4,-2,2,-2,2,-2,2,-4,4,-3,3,0,-4,4,0,0,-3,3,0,-3,3,-4,4,
    -4,4,-3,3,-1,1,-4,4,-4,4,-1,1,-1,1,-3,3,-3,3,-1,1,-3,3,-4,4,-4,4,
    -3,3,-2,2,-4,4,-4,4,-2,2,-2,2,-3,3,-3,3,-2,2,-3,3,-4,4,-4,4,-3,3,
    -3,3,-3,3,-4,4,0,-4,4,0,-4,4,-4,4,-1,1,-4,4,-4,4,-1,1,-4,4,-4,4,
    -2,2,-4,4,-4,4,-2,2,-4,4,-4,4,-3,3,-4,4,-4,4,-3,3,5,0,0,-5,5,-1,
    1,0,-5,5,0,0,-1,1,0,-1,1,-5,5,-5,5,-1,1,-1,1,-1,1,-5,5,-2,2,0,-5,
    5,0,0,-2,2,0,-2,2,-5,5,-5,5,-2,2,-1,1,-5,5,-5,5,-1,1,-1,1,-2,2,
    -2,2,-1,1,-2,2,-5,5,-5,5,-2,2,-2,2,-2,2,-5,5,-3,3,0,-5,5,0,0,-3,
    3,0,-3,3,-5,5,-5,5,-3,3,-1,1,-5,5,-5,5,-1,1,-1,1,-3,3,-3,3,-1,1,
    -3,3,-5,5,-5,5,-3,3,-2,2,-5,5,-5,5,-2,2,-2,2,-3,3,-3,3,-2,2,-3,3,
    -5,5,-5,5,-3,3,-3,3,-3,3};
    
  int transy[508]={0,0,1,0,1,1,-1,0,0,1,-1,-1,1,1,0,2,0,1,1,2,2,-2,0,0,2,
    -1,0,0,1,-2,-2,-1,-1,1,1,2,2,-1,-1,1,1,2,2,-2,0,0,2,-2,-2,2,2,-2,
    -2,-1,-1,1,1,2,2,-2,-2,2,2,0,3,0,1,1,3,3,-3,0,0,3,-1,0,0,1,-3,-3,
    -1,-1,1,1,3,3,-1,-1,1,1,2,2,3,3,-3,0,0,3,-2,0,0,2,-3,-3,-2,-2,2,
    2,3,3,-3,-3,-1,-1,1,1,3,3,-2,-2,-1,-1,1,1,2,2,-3,-3,-2,-2,2,2,3,
    3,-2,-2,2,2,3,3,-3,0,0,3,-3,-3,3,3,-3,-3,-1,-1,1,1,3,3,-3,-3,3,3,
    -3,-3,-2,-2,2,2,3,3,-3,-3,3,3,0,4,0,1,1,4,4,-4,0,0,4,-1,0,0,1,-4,
    -4,-1,-1,1,1,4,4,-1,-1,1,1,2,2,4,4,-4,0,0,4,-2,0,0,2,-4,-4,-2,-2,
    2,2,4,4,-4,-4,-1,-1,1,1,4,4,-2,-2,-1,-1,1,1,2,2,-4,-4,-2,-2,2,2,
    4,4,-2,-2,2,2,3,3,4,4,-4,0,0,4,-3,0,0,3,-4,-4,-3,-3,3,3,4,4,-4,
    -4,-1,-1,1,1,4,4,-3,-3,-1,-1,1,1,3,3,-4,-4,-3,-3,3,3,4,4,-4,-4,
    -2,-2,2,2,4,4,-3,-3,-2,-2,2,2,3,3,-4,-4,-3,-3,3,3,4,4,-3,-3,3,3,
    4,4,-4,0,0,4,-4,-4,4,4,-4,-4,-1,-1,1,1,4,4,-4,-4,4,4,-4,-4,-2,-2,
    2,2,4,4,-4,-4,4,4,-4,-4,-3,-3,3,3,4,4,0,5,0,1,1,5,5,-5,0,0,5,-1,
    0,0,1,-5,-5,-1,-1,1,1,5,5,-1,-1,1,1,2,2,5,5,-5,0,0,5,-2,0,0,2,-5,
    -5,-2,-2,2,2,5,5,-5,-5,-1,-1,1,1,5,5,-2,-2,-1,-1,1,1,2,2,-5,-5,
    -2,-2,2,2,5,5,-2,-2,2,2,3,3,5,5,-5,0,0,5,-3,0,0,3,-5,-5,-3,-3,3,
    3,5,5,-5,-5,-1,-1,1,1,5,5,-3,-3,-1,-1,1,1,3,3,-5,-5,-3,-3,3,3,5,
    5,-5,-5,-2,-2,2,2,5,5,-3,-3,-2,-2,2,2,3,3,-5,-5,-3,-3,3,3,5,5,-3,
    -3,3,3};
    
  int transz[508]={0,0,0,1,0,0,1,1,1,1,1,1,1,1,0,0,2,0,0,0,0,1,1,1,1,2,2,2,
    2,1,1,1,1,1,1,1,1,2,2,2,2,0,0,2,2,2,2,1,1,1,1,2,2,2,2,2,2,2,2,2,
    2,2,2,0,0,3,0,0,0,0,1,1,1,1,3,3,3,3,1,1,1,1,1,1,1,1,3,3,3,3,0,0,
    0,0,2,2,2,2,3,3,3,3,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,
    3,3,2,2,2,2,2,2,2,2,3,3,3,3,0,0,3,3,3,3,1,1,1,1,3,3,3,3,3,3,3,3,
    2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,0,0,4,0,0,0,0,1,1,1,1,4,4,4,4,1,
    1,1,1,1,1,1,1,4,4,4,4,0,0,0,0,2,2,2,2,4,4,4,4,1,1,1,1,1,1,1,1,2,
    2,2,2,2,2,2,2,4,4,4,4,4,4,4,4,2,2,2,2,2,2,2,2,4,4,4,4,0,0,0,0,3,
    3,3,3,4,4,4,4,1,1,1,1,1,1,1,1,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,2,
    2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,4,
    4,4,4,0,0,4,4,4,4,1,1,1,1,4,4,4,4,4,4,4,4,2,2,2,2,4,4,4,4,4,4,4,
    4,3,3,3,3,4,4,4,4,4,4,4,4,0,0,5,0,0,0,0,1,1,1,1,5,5,5,5,1,1,1,1,
    1,1,1,1,5,5,5,5,0,0,0,0,2,2,2,2,5,5,5,5,1,1,1,1,1,1,1,1,2,2,2,2,
    2,2,2,2,5,5,5,5,5,5,5,5,2,2,2,2,2,2,2,2,5,5,5,5,0,0,0,0,3,3,3,3,
    5,5,5,5,1,1,1,1,1,1,1,1,3,3,3,3,3,3,3,3,5,5,5,5,5,5,5,5,2,2,2,2,
    2,2,2,2,3,3,3,3,3,3,3,3,5,5,5,5,5,5,5,5,3,3,3,3,3,3,3,3,5,5,5,5};
    
  int cellCheck,nlcx,nlcy,nlcz;
  int i,j,k,kv,exclude,nalloc,incr=262144;
  double r,cutnb,dnlcx,dnlcy,dnlcz,delta[3];
  
  cutnb=simulCond->cutoff+simulCond->delr;
  
  ff->verPair=(int*)malloc((simulCond->natom-1)*sizeof(*(ff->verPair)));
  for(i=0;i<simulCond->natom-1;i++)
    ff->verPair[i]=0;

//   nalloc=((simulCond->natom*simulCond->natom/64)-(simulCond->natom/8))/2;
//   nalloc = (int) simulCond->natom * ( simulCond->natom/8 - 1  ) / 2 ;
//   nalloc = simulCond->natom*4./3.*PI*X3(cutnb)/boxvolume;
  nalloc = (int) simulCond->natom * 1024 ;
  ff->verList=(int*)malloc(nalloc*sizeof(*(ff->verList)));
//   for(i=0;i<nalloc;i++)
//     ff->verList[i]=0;

  for(i=0;i<simulCond->natom;i++)
  {
    link[i]=-1;
  }
  
  if(simulCond->linkRatio==1)
    cellCheck=14;
  else if(simulCond->linkRatio==2)
    cellCheck=63;
  else if(simulCond->linkRatio==3)
    cellCheck=156;
  else if(simulCond->linkRatio==4)
    cellCheck=307;
  else if(simulCond->linkRatio==5)
    cellCheck=508;
  else
    error(411);

  nlcx=int(box->a1*(double)simulCond->linkRatio/cutnb);
  nlcy=int(box->b2*(double)simulCond->linkRatio/cutnb);
  nlcz=int(box->c3*(double)simulCond->linkRatio/cutnb);
  
  dnlcx=(double)nlcx;
  dnlcy=(double)nlcy;
  dnlcz=(double)nlcz;
  
  enoughLinkCells=1;
  if(nlcx<2*simulCond->ratio+1)
    enoughLinkCells=0;
  if(nlcy<2*simulCond->ratio+1)
    enoughLinkCells=0;
  if(nlcz<2*simulCond->ratio+1)
    enoughLinkCells=0;
  if(!enoughLinkCells)
    error(412);
  
  ff->ncells=nlcx*nlcy*nlcz;
  
  if(ff->ncells>simulCond->maxcells)
    error(413);
  
  for(i=0;i<simulCond->natom;i++)
  {
    head[i]=-1;
  }
  
  for(i=0;i<simulCond->natom;i++)
  {
    ix=(int)dnlcx*atom[i].x/box->a1;
    ix=MIN(ix,nlcx-1);
    iy=(int)dnlcy*atom[i].y/box->b2;
    iy=MIN(iy,nlcy-1);
    iz=(int)dnlcz*atom[i].z/box->c3;
    iz=MIN(iz,nlcz-1);
    
    icell=ix+nlcx*(iy+nlcy*iz);
    
    j=head[icell];
    head[icell]=i;
    link[i]=j;
  }
  
  ix=0;
  iy=0;
  iz=0;
  
  ff->npr=0;
  
  for(k=0;k<ff->ncells;k++)
  {
    i=head[k];
    
    if(i>-1)
    {
      for(l=0;l<cellCheck;l++)
      {
      
	cx=0.;
	cy=0.;
	cz=0.;
	
	jx=ix+transx(l);
	jy=iy+transy(l);
	jz=iz+transz(l);
	
	if(jx>nlcx-1)
	{
	  jx=jx-nlcx;
	  cx=1.;
	}
	else if(jx<0)
	{
	  jx=jx+nlcx
	  cx=-1.;
	}
	
	if(jy>nlcy-1)
	{
	  jy=jy-nlcy;
	  cy=1.;
	}
	else if(jy<0)
	{
	  jy=jy+nlcy
	  cy=-1.;
	}
	
	if(jz>nlcz-1)
	{
	  jz=jz-nlcz;
	  cz=1.;
	}
	else if(jz<0)
	{
	  jz=jz+nlcz
	  cz=-1.;
	}
	
	ll=jx+nlcx*(jy+nlcy*jz);
	j=head(ll);
	
	if(j>-1)
	{
	  while(i!=-1)
	  {
	    ff->npr=0;
	    if(k==ll)j=link[i];
	    
	    if(j>-1)
	    {
	      while(j!=-1)
	      {
		xt=atom[j].x-atom[i].x+(cx*box->a1);
		yt=atom[j].y-atom[i].y+(cy*box->b2);
		zt=atom[j].z-atom[i].z+(cz*box->c3);
		
		r=X2(xt)+X2(yt)+X2(zt);
		
		if(r<=cutnb)
		{
		  exclude=0;
		  for (k=0;k<simulCond->excludeNum[i];k++)
		  {
		    if(simulCond->excludeAtom[i][k]==j)
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
		      ff->verList[i]=(int*)realloc(ff->verList[i],nalloc*sizeof(**(ff->verList)));
		    }
		    ff->verList[i][ff->verPair[i]]=j;
		    ff->verPair[i]++;
		  } //exclude
		}//cutnb
		
		j=link[j];
	      
	      }//while j
	      
	    }//if j>-1
	    
	    j=head[ll];
	    i=link[i];
	    
	  }//while i
	}//if j>-1
      }//for checkcell
    }// if i >-
    
    ix=ix+1;
    if(ix>nlcx-1)
    {       
      ix=1;
      iy=iy+1;
            
      if(iy>nlcy-1)
      {        
	iy=1;
	iz=iz+1;
      }
      
    }
    
  }//for ncells
  
  
  ff->verList[i]=(int*)realloc(ff->verList,ff->npr*sizeof(*(ff->verList)));

}*/
