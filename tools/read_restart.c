#include <stdio.h>
#include <stdlib.h>

typedef struct
{
    int step,nSteps;
    real timeStep;

    int nAtom,nDegFree;

    int nBond,nAngle,nDihedral,nImproper,nUb;

    int nConst,maxCycle;
    real tolShake;

    real chargeConst,cutOff,cutOn,delr,scal14;

    real temp0,press0,kinTemp0;

    real tolMinim,maxminsiz,maxminst;
} PARAM;

typedef struct
{
    real tot,pot,kin;
    real elec,vdw;
    real bond,ang,ub,dihe,impr;
    real conint,consv;
    real virelec,virvdw,virbond,virub;
    real virshake,virpot,virtot;
} ENERGY;

typedef struct
{
    real tauT,chiT;

    real tauP,chiP,compress;
} BATH;

typedef struct
{
    char label[5],segn[5],resn[5];
    int type,resi,ires,inconst;

    real x,y,z;
    real vx,vy,vz;
    real fx,fy,fz;

    real m,q;
} ATOM;

#define _FORMATN1_ "%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\n\n"
#define _FORMATL1_ "%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\n"

#define _FORMATN2_ "%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\n\n"
#define _FORMATL2_ "%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\n"

int main(int argc, char* argv[])
{
    FILE *restFile=fopen(argv[1],"rb");

    FILE *corFile=fopen("COR","w");
    FILE *velFile=fopen("VEL","w");
    FILE *forFile=fopen("FOR","w");

    PARAM param;
    ENERGY ener;
    BATH bath;
    ATOM *atom=NULL;

    size_t ret;
    real wei;
    int i;

    char label1[14][7]= { {"Etot"} , {"Ekin"}  , {"Epot"} , {"Ecoul"} , {"Evdw"} ,
        {"Ebond"} , {"Eangle"} , {"Eub"} , {"Edihe"} , {"Eimpr"} ,
        {"Conint"} , {"Consv"} , {"chiT"} , {"chiP"}
    };

    char label2[7][7]= { {"Virtot"} , {"Virpot"} , {"Virele"} , {"Virvdw"} ,
        {"Virbnd"} , {"Virub"} , {"Virshk"}
    };

    ret=fread(&(param.step),sizeof(int),1,restFile);
    if(ret!=1)
        printf("Problem while reading the file.\n");

    ret=fread(&(param.nAtom),sizeof(real),1,restFile);
    if(ret!=1)
        printf("Problem while reading the file.\n");

    atom=(ATOM*)malloc(param.nAtom*sizeof(ATOM));

    ret=fread(atom,sizeof(ATOM),param.nAtom,restFile);
    if(ret!=param.nAtom)
        printf("Problem while reading the file.\n");

    ret=fread(&ener,sizeof(ENERGY),1,restFile);
    if(ret!=1)
        printf("Problem while reading the file.\n");

    ret=fread(&(bath.chiT),sizeof(real),1,restFile);
    if(ret!=1)
        printf("Problem while reading the file.\n");

    ret=fread(&(bath.chiP),sizeof(real),1,restFile);
    if(ret!=1)
        printf("Problem while reading the file.\n");

    for(i=0; i<param.nAtom; i++)
    {
        fprintf(corFile,"%5d%5d %-4s %-4s%10.5lf%10.5lf%10.5lf %-4s %-4d%10.5lf\n",
                i+1,atom[i].ires,atom[i].resn,atom[i].label,atom[i].x,atom[i].y,
                atom[i].z,atom[i].segn,atom[i].resi,wei);

    }

    for(i=0; i<param.nAtom; i++)
    {
        fprintf(velFile,"%5d%5d %-4s %-4s%10.5lf%10.5lf%10.5lf %-4s %-4d%10.5lf\n",
                i+1,atom[i].ires,atom[i].resn,atom[i].label,atom[i].vx,atom[i].vy,
                atom[i].vz,atom[i].segn,atom[i].resi,wei);

    }

    for(i=0; i<param.nAtom; i++)
    {
        fprintf(forFile,"%5d%5d %-4s %-4s%10.5lf%10.5lf%10.5lf %-4s %-4d%10.5lf\n",
                i+1,atom[i].ires,atom[i].resn,atom[i].label,atom[i].fx,atom[i].fy,
                atom[i].fz,atom[i].segn,atom[i].resi,wei);

    }

    printf(_FORMATL1_,label1[0],label1[1],label1[2],label1[3],label1[4],
           label1[5],label1[6],label1[7],label1[8],label1[9],
           label1[10],label1[11],label1[12],label1[13]);

    printf(_FORMATN1_,ener.tot,ener.pot,ener.kin,ener.elec,ener.vdw,
           ener.bond,ener.ang,ener.ub,ener.dihe,ener.impr,
           ener.conint,ener.consv,bath.chiT,bath.chiP);

    printf(_FORMATL2_,label2[0],label2[1],label2[2],label2[3],
           label2[4],label2[5],label2[6]);

    printf(_FORMATN2_,ener.virtot,ener.virpot,ener.virelec,ener.virvdw,
           ener.virbond,ener.virub,ener.virshake);

    fclose(restFile);

}
