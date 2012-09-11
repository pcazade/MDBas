#include <stdio.h>
#include <float.h>
#include <math.h>
#include "global.h"
#include "utils.h"

void bond_energy(ATOM *atom,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond,PBC *box)
{
  int i,j,ll;
  double r,fx,fy,fz,dbond,virbond=0.;
  double delta[3],stress[6]={0.};
  
  for(ll=0;ll<ff->nBond;ll++)
  {
    
    i=simulCond->iBond[ll][0];
    j=simulCond->iBond[ll][1];
    
    r=distance(i,j,atom,delta,simulCond,box);
    
    ener->bond+=0.5*ff->parmBond[ll][0]*X2(r-ff->parmBond[ll][1]);
    dbond=ff->parmBond[ll][0]*(r-ff->parmBond[ll][1]);
    
    virbond+=dbond*r;
    
    fx=dbond*delta[0]/r;
    fy=dbond*delta[1]/r;
    fz=dbond*delta[2]/r;
    
    atom[i].fx+=fx;
    atom[i].fy+=fy;
    atom[i].fz+=fz;
    
    atom[j].fx+=-fx;
    atom[j].fy+=-fy;
    atom[j].fz+=-fz;
    
    /*stress[0]-=fx*delta[0];
    stress[1]-=fy*delta[0];
    stress[2]-=fz*delta[0];
    stress[3]-=fy*delta[1];
    stress[4]-=fz*delta[1];
    stress[5]-=fz*delta[2];*/
    
  }
  
  ener->virbond+=virbond;
  
  /*box->stress1+=stress[0];
  box->stress2+=stress[1];
  box->stress3+=stress[2];
  box->stress4+=stress[1];
  box->stress5+=stress[3];
  box->stress6+=stress[4];
  box->stress7+=stress[2];
  box->stress8+=stress[4];
  box->stress9+=stress[5];*/
  
}

void ub_energy(ATOM *atom,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond,PBC *box)
{
  int i,j,ll;
  double r,fx,fy,fz,dub,virub=0.;
  double delta[3],stress[6]={0.};
  
  for(ll=0;ll<ff->nUb;ll++)
  {
    i=simulCond->iUb[ll][0];
    j=simulCond->iUb[ll][1];
    
    r=distance(i,j,atom,delta,simulCond,box);
    
    ener->ub+=0.5*ff->parmUb[ll][0]*X2(r-ff->parmUb[ll][1]);
    dub=ff->parmUb[ll][0]*(r-ff->parmUb[ll][1]);
    
    virub+=dub*r;

    fx=dub*delta[0]/r;
    fy=dub*delta[1]/r;
    fz=dub*delta[2]/r;
    
    atom[i].fx+=fx;
    atom[i].fy+=fy;
    atom[i].fz+=fz;
    
    atom[j].fx+=-fx;
    atom[j].fy+=-fy;
    atom[j].fz+=-fz;
    
    /*stress[0]-=fx*delta[0];
    stress[1]-=fy*delta[0];
    stress[2]-=fz*delta[0];
    stress[3]-=fy*delta[1];
    stress[4]-=fz*delta[1];
    stress[5]-=fz*delta[2];*/
    
  }
  
  ener->virub+=virub;
  
  /*box->stress1+=stress[0];
  box->stress2+=stress[1];
  box->stress3+=stress[2];
  box->stress4+=stress[1];
  box->stress5+=stress[3];
  box->stress6+=stress[4];
  box->stress7+=stress[2];
  box->stress8+=stress[4];
  box->stress9+=stress[5];*/
  
}

void angle_energy(ATOM *atom,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond,PBC *box)
{
  int i,j,k,ll;
  double dangle,rab,rbc,cost,sint,theta;
  double dab[3],dbc[3],stress[6];
  double fxa,fya,fza,fxc,fyc,fzc;
  
  for(ll=0;ll<ff->nAngle;ll++)
  {
    
    i=simulCond->iAngle[ll][0];
    j=simulCond->iAngle[ll][1];
    k=simulCond->iAngle[ll][2];
    
    rab=distance(j,i,atom,dab,simulCond,box);
    rbc=distance(j,k,atom,dbc,simulCond,box);
    
    cost=(dab[0]*dbc[0]+dab[1]*dbc[1]+dab[2]*dbc[2])/(rab*rbc);
    sint=MAX(1.e-8,sqrt(1.-(cost*cost)));
    theta=acos(cost);
    
    ener->ang+=0.5*ff->parmAngle[ll][0]*X2(theta-ff->parmAngle[ll][1]);
    dangle=ff->parmAngle[ll][0]*(theta-ff->parmAngle[ll][1])/sint;
    
    fxa=dangle*(dbc[0]/rbc-dab[0]*cost/rab)/rab;
    fya=dangle*(dbc[1]/rbc-dab[1]*cost/rab)/rab;
    fza=dangle*(dbc[2]/rbc-dab[2]*cost/rab)/rab;
    
    fxc=dangle*(dab[0]/rab-dbc[0]*cost/rbc)/rbc;
    fyc=dangle*(dab[1]/rab-dbc[1]*cost/rbc)/rbc;
    fzc=dangle*(dab[2]/rab-dbc[2]*cost/rbc)/rbc;
    
    atom[i].fx+=fxa;
    atom[i].fy+=fya;
    atom[i].fz+=fza;
    
    atom[j].fx+=-fxa-fxc;
    atom[j].fy+=-fya-fyc;
    atom[j].fz+=-fza-fzc;
    
    atom[k].fx+=fxc;
    atom[k].fy+=fyc;
    atom[k].fz+=fzc;
    
    /*stress[0]=rab*fxa*dab[0]+rbc*fxc*dbc[0];
    stress[1]=rab*fya*dab[0]+rbc*fyc*dbc[0];
    stress[2]=rab*fza*dab[0]+rbc*fzc*dbc[0];
    stress[3]=rab*fya*dab[1]+rbc*fyc*dbc[1];
    stress[4]=rab*fza*dab[1]+rbc*fzc*dbc[1];
    stress[5]=rab*fza*dab[2]+rbc*fzc*dbc[2];*/
    
  }
  
  /*box->stress1+=stress[0];
  box->stress2+=stress[1];
  box->stress3+=stress[2];
  box->stress4+=stress[1];
  box->stress5+=stress[3];
  box->stress6+=stress[4];
  box->stress7+=stress[2];
  box->stress8+=stress[4];
  box->stress9+=stress[5];*/
  
}

void dihedral_energy(ATOM *atom,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond,PBC *box)
{
  int i,j,k,l,ll,nd,ind;
  double pi,twopi;
  double edihe=0.,ddihe=0.;
  double cosp,sinp,phi;
  double /*rab,*/rbc,/*rcd,*/rpb,rpc,r2pb,r2pc,pbpc;
  double dab[3],dbc[3],dcd[3],/*dac[3],*/pb[3],pc[3],stress[6]={0.};
  double fax,fay,faz,fbx,fby,fbz,fcx,fcy,fcz,fdx,fdy,fdz;
  
  pi=atan(-1.);
  twopi=2.*pi;
  
  for(ll=0;ll<ff->nDihedral;ll++)
  {
    i=simulCond->iDihedral[ll][0];
    j=simulCond->iDihedral[ll][1];
    k=simulCond->iDihedral[ll][2];
    l=simulCond->iDihedral[ll][3];
    
    /*rab=*/distance(j,i,atom,dab,simulCond,box);
    rbc=distance(k,j,atom,dbc,simulCond,box);
    /*rcd=*/distance(l,k,atom,dcd,simulCond,box);
    
    /*dac[0]=dab[0]+dbc[0];
    dac[1]=dab[1]+dbc[1];
    dac[2]=dab[2]+dbc[2];*/
    
// construct first dihedral vector
    
    pb[0]=dab[1]*dbc[2]-dab[2]*dbc[1];
    pb[1]=dab[2]*dbc[0]-dab[0]*dbc[2];
    pb[2]=dab[0]*dbc[1]-dab[1]*dbc[0];
    
    r2pb=X2(pb[0])+X2(pb[1])+X2(pb[2]);
    rpb=sqrt(r2pb);
    
// construct second dihedral vector
    
    pc[0]=dbc[1]*dcd[2]-dbc[2]*dcd[1];
    pc[1]=dbc[2]*dcd[0]-dbc[0]*dcd[2];
    pc[2]=dbc[0]*dcd[1]-dbc[1]*dcd[0];
    
    r2pc=X2(pc[0])+X2(pc[1])+X2(pc[2]);
    rpc=sqrt(r2pc);
    
// determine dihedral angle 
    
    pbpc=pb[0]*pc[0]+pb[1]*pc[1]+pb[2]*pc[2];
    cosp=pbpc/(rpb*rpc);
    
    sinp=(dbc[0]*(pc[1]*pb[2]-pc[2]*pb[1])+dbc[1]*(pb[0]*pc[2]-pb[2]*pc[0])+
          dbc[2]*(pc[0]*pb[1]-pc[1]*pb[0]))/(rpb*rpc*rbc);
        
    phi=atan2(sinp,cosp);
    
// avoid singularity in sinp
    
    if(sinp>=0.)
    {
      sinp=MAX(DBL_EPSILON,fabs(sinp));
    }
    else
    {
      sinp=-(MAX(DBL_EPSILON,fabs(sinp)));
    }
    
// calculate potential energy and scalar force term
    
    if(simulCond->diheType[ll]==COSDIH)
    {

// key=1 for cosine dihedral
      
      edihe=0.;
      ddihe=0.;
      
      for(nd=0;nd<ff->nParmDihe[ll];nd++)
      {
	ind=3*nd;
	edihe+=ff->parmDihe[ll][ind]*(1.+cos(ff->parmDihe[ll][ind+1]*phi-ff->parmDihe[ll][ind+2]));
	
	ddihe+=-ff->parmDihe[ll][ind]*ff->parmDihe[ll][ind+1]*
	      sin(ff->parmDihe[ll][ind+1]*phi-ff->parmDihe[ll][ind+2])/(rpb*rpc*sinp);
      }
    }
    
    else if(simulCond->diheType[ll]==HARMDIH)
    {
      
// key=2 for harmonic improper dihedral
      
      phi=phi-ff->parmDihe[ll][1];
      phi=phi-nint(phi/twopi)*twopi;
      
      edihe=0.5*ff->parmDihe[ll][0]*(phi*phi);
      
      ddihe=ff->parmDihe[ll][0]*phi/(rpb*rpc*sinp);
    }
    
// calculate potential energy
      
    ener->dihe+=edihe;
    
    fax=ddihe*((-pc[1]*dbc[2]+pc[2]*dbc[1])-pbpc*(-pb[1]*dbc[2]+pb[2]*dbc[1])/r2pb);
    fay=ddihe*(( pc[0]*dbc[2]-pc[2]*dbc[0])-pbpc*( pb[0]*dbc[2]-pb[2]*dbc[0])/r2pb);
    faz=ddihe*((-pc[0]*dbc[1]+pc[1]*dbc[0])-pbpc*(-pb[0]*dbc[1]+pb[1]*dbc[0])/r2pb);
	
    fcx=ddihe*((-pc[1]*dab[2]+pc[2]*dab[1])-pbpc*(-pb[1]*dab[2]+pb[2]*dab[1])/r2pb);
    fcy=ddihe*(( pc[0]*dab[2]-pc[2]*dab[0])-pbpc*( pb[0]*dab[2]-pb[2]*dab[0])/r2pb);
    fcz=ddihe*((-pc[0]*dab[1]+pc[1]*dab[0])-pbpc*(-pb[0]*dab[1]+pb[1]*dab[0])/r2pb);
	
    fbx=ddihe*((-pb[1]*dcd[2]+pb[2]*dcd[1])-pbpc*(-pc[1]*dcd[2]+pc[2]*dcd[1])/r2pc);
    fby=ddihe*(( pb[0]*dcd[2]-pb[2]*dcd[0])-pbpc*( pc[0]*dcd[2]-pc[2]*dcd[0])/r2pc);
    fbz=ddihe*((-pb[0]*dcd[1]+pb[1]*dcd[0])-pbpc*(-pc[0]*dcd[1]+pc[1]*dcd[0])/r2pc);
	
    fdx=ddihe*((-pb[1]*dbc[2]+pb[2]*dbc[1])-pbpc*(-pc[1]*dbc[2]+pc[2]*dbc[1])/r2pc);
    fdy=ddihe*(( pb[0]*dbc[2]-pb[2]*dbc[0])-pbpc*( pc[0]*dbc[2]-pc[2]*dbc[0])/r2pc);
    fdz=ddihe*((-pb[0]*dbc[1]+pb[1]*dbc[0])-pbpc*(-pc[0]*dbc[1]+pc[1]*dbc[0])/r2pc);

    atom[i].fx+=fax;
    atom[i].fy+=fay;
    atom[i].fz+=faz;
	
    atom[j].fx+=-fax-fcx+fbx;
    atom[j].fy+=-fay-fcy+fby;
    atom[j].fz+=-faz-fcz+fbz;
	
    atom[k].fx+=fcx-fbx-fdx;
    atom[k].fy+=fcy-fby-fdy;
    atom[k].fz+=fcz-fbz-fdz;
	
    atom[l].fx+=fdx;
    atom[l].fy+=fdy;
    atom[l].fz+=fdz;
    
    /*stress[0]+=dab[0]*fax+dbc[0]*(fbx-fcx)-dcd[0]*fdx;
    stress[1]+=dab[1]*fax+dbc[1]*(fbx-fcx)-dcd[1]*fdx; 
    stress[2]+=dab[2]*fax+dbc[2]*(fbx-fcx)-dcd[2]*fdx; 
    stress[3]+=dab[1]*fay+dbc[1]*(fby-fcy)-dcd[1]*fdy; 
    stress[4]+=dab[1]*faz+dbc[1]*(fbz-fcz)-dcd[1]*fdz; 
    stress[5]+=dab[2]*faz+dbc[2]*(fbz-fcz)-dcd[2]*fdz; */
      
  }
  
  /*box->stress1+=stress[0];
  box->stress2+=stress[1];
  box->stress3+=stress[2];
  box->stress4+=stress[1];
  box->stress5+=stress[3];
  box->stress6+=stress[4];
  box->stress7+=stress[2];
  box->stress8+=stress[4];
  box->stress9+=stress[5];*/
  
}

void improper_energy(ATOM *atom,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond,PBC *box)
{
  int i,j,k,l,ll;
  double pi,twopi;
  double edihe=0.,ddihe=0.;
  double cosp,sinp,phi;
  double /*rab,*/rbc,/*rcd,*/rpb,rpc,r2pb,r2pc,pbpc;
  double dab[3],dbc[3],dcd[3],/*dac[3],*/pb[3],pc[3],stress[6]={0.};
  double fax,fay,faz,fbx,fby,fbz,fcx,fcy,fcz,fdx,fdy,fdz;
  
  pi=atan(-1.);
  twopi=2.*pi;
  
  for(ll=0;ll<ff->nImproper;ll++)
  {
    i=simulCond->iImproper[ll][0];
    j=simulCond->iImproper[ll][1];
    k=simulCond->iImproper[ll][2];
    l=simulCond->iImproper[ll][3];
    
    /*rab=*/distance(j,i,atom,dab,simulCond,box);
    rbc=distance(k,j,atom,dbc,simulCond,box);
    /*rcd=*/distance(l,k,atom,dcd,simulCond,box);
    
    /*dac[0]=dab[0]+dbc[0];
    dac[1]=dab[1]+dbc[1];
    dac[2]=dab[2]+dbc[2];*/
    
// construct first dihedral vector
    
    pb[0]=dab[1]*dbc[2]-dab[2]*dbc[1];
    pb[1]=dab[2]*dbc[0]-dab[0]*dbc[2];
    pb[2]=dab[0]*dbc[1]-dab[1]*dbc[0];
    
    r2pb=X2(pb[0])+X2(pb[1])+X2(pb[2]);
    rpb=sqrt(r2pb);
    
// construct second dihedral vector
    
    pc[0]=dbc[1]*dcd[2]-dbc[2]*dcd[1];
    pc[1]=dbc[2]*dcd[0]-dbc[0]*dcd[2];
    pc[2]=dbc[0]*dcd[1]-dbc[1]*dcd[0];
    
    r2pc=X2(pc[0])+X2(pc[1])+X2(pc[2]);
    rpc=sqrt(r2pc);
    
// determine dihedral angle 
    
    pbpc=pb[0]*pc[0]+pb[1]*pc[1]+pb[2]*pc[2];
    cosp=pbpc/(rpb*rpc);
    
    sinp=(dbc[0]*(pc[1]*pb[2]-pc[2]*pb[1])+dbc[1]*(pb[0]*pc[2]-pb[2]*pc[0])+
          dbc[2]*(pc[0]*pb[1]-pc[1]*pb[0]))/(rpb*rpc*rbc);
        
    phi=atan2(sinp,cosp);
    
// avoid singularity in sinp
    
    if(sinp>=0)
    {
      sinp=MAX(DBL_EPSILON,fabs(sinp));
    }
    else
    {
      sinp=-(MAX(DBL_EPSILON,fabs(sinp)));
    }
    
// calculate potential energy and scalar force term
    
    if(simulCond->imprType[ll]==COSDIH)
    {

// key=1 for cosine dihedral
      
      edihe=ff->parmImpr[ll][0]*(1.+cos(ff->parmImpr[ll][1]*phi-ff->parmImpr[ll][2]));
      ddihe=-ff->parmImpr[ll][0]*ff->parmImpr[ll][1]*
            sin(ff->parmImpr[ll][1]*phi-ff->parmImpr[ll][2])/(rpb*rpc*sinp);  
    }
    
    else if(simulCond->imprType[ll]==HARMDIH)
    {
      
// key=2 for harmonic improper dihedral
      
      phi=phi-ff->parmImpr[ll][1];
      phi=phi-nint(phi/twopi)*twopi;
      edihe=0.5*ff->parmImpr[ll][0]*(phi*phi);
      ddihe=ff->parmImpr[ll][0]*phi/(rpb*rpc*sinp);
    }
    
// calculate potential energy
      
    ener->impr+=edihe;
    
    fax=ddihe*((-pc[1]*dbc[2]+pc[2]*dbc[1])-pbpc*(-pb[1]*dbc[2]+pb[2]*dbc[1])/r2pb);
    fay=ddihe*(( pc[0]*dbc[2]-pc[2]*dbc[0])-pbpc*( pb[0]*dbc[2]-pb[2]*dbc[0])/r2pb);
    faz=ddihe*((-pc[0]*dbc[1]+pc[1]*dbc[0])-pbpc*(-pb[0]*dbc[1]+pb[1]*dbc[0])/r2pb);
	
    fcx=ddihe*((-pc[1]*dab[2]+pc[2]*dab[1])-pbpc*(-pb[1]*dab[2]+pb[2]*dab[1])/r2pb);
    fcy=ddihe*(( pc[0]*dab[2]-pc[2]*dab[0])-pbpc*( pb[0]*dab[2]-pb[2]*dab[0])/r2pb);
    fcz=ddihe*((-pc[0]*dab[1]+pc[1]*dab[0])-pbpc*(-pb[0]*dab[1]+pb[1]*dab[0])/r2pb);
	
    fbx=ddihe*((-pb[1]*dcd[2]+pb[2]*dcd[1])-pbpc*(-pc[1]*dcd[2]+pc[2]*dcd[1])/r2pc);
    fby=ddihe*(( pb[0]*dcd[2]-pb[2]*dcd[0])-pbpc*( pc[0]*dcd[2]-pc[2]*dcd[0])/r2pc);
    fbz=ddihe*((-pb[0]*dcd[1]+pb[1]*dcd[0])-pbpc*(-pc[0]*dcd[1]+pc[1]*dcd[0])/r2pc);
	
    fdx=ddihe*((-pb[1]*dbc[2]+pb[2]*dbc[1])-pbpc*(-pc[1]*dbc[2]+pc[2]*dbc[1])/r2pc);
    fdy=ddihe*(( pb[0]*dbc[2]-pb[2]*dbc[0])-pbpc*( pc[0]*dbc[2]-pc[2]*dbc[0])/r2pc);
    fdz=ddihe*((-pb[0]*dbc[1]+pb[1]*dbc[0])-pbpc*(-pc[0]*dbc[1]+pc[1]*dbc[0])/r2pc);
	
    atom[i].fx+=fax;
    atom[i].fy+=fay;
    atom[i].fz+=faz;
	
    atom[j].fx+=-fax-fcx+fbx;
    atom[j].fy+=-fay-fcy+fby;
    atom[j].fz+=-faz-fcz+fbz;
	
    atom[k].fx+=fcx-fbx-fdx;
    atom[k].fy+=fcy-fby-fdy;
    atom[k].fz+=fcz-fbz-fdz;
	
    atom[l].fx+=fdx;
    atom[l].fy+=fdy;
    atom[l].fz+=fdz;
    
    /*stress[0]+=dab[0]*fax+dbc[0]*(fbx-fcx)-dcd[0]*fdx;
    stress[1]+=dab[1]*fax+dbc[1]*(fbx-fcx)-dcd[1]*fdx; 
    stress[2]+=dab[2]*fax+dbc[2]*(fbx-fcx)-dcd[2]*fdx; 
    stress[3]+=dab[1]*fay+dbc[1]*(fby-fcy)-dcd[1]*fdy; 
    stress[4]+=dab[1]*faz+dbc[1]*(fbz-fcz)-dcd[1]*fdz; 
    stress[5]+=dab[2]*faz+dbc[2]*(fbz-fcz)-dcd[2]*fdz;*/
        
  }
  
  /*box->stress1+=stress[0];
  box->stress2+=stress[1];
  box->stress3+=stress[2];
  box->stress4+=stress[1];
  box->stress5+=stress[3];
  box->stress6+=stress[4];
  box->stress7+=stress[2];
  box->stress8+=stress[4];
  box->stress9+=stress[5];*/
  
}
