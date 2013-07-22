#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <float.h>
#include <stddef.h>
#include <stdarg.h>
#include <assert.h>

typedef struct
{
    char label[5],segn[5],resn[5];
    int type,resi,ires,inconst,frozen;

    double x,y,z;
    double vx,vy,vz;
    double fx,fy,fz;

    double m,rm,q;
} ATOM;

typedef struct
{
    char **types;
    int *typesNum;
    int nTypes,nBondTypes,nAngTypes,nUbTypes,nDiheTypes,nImprTypes,nNonBonded;
    int **bondTypes,**angTypes,**ubTypes,*nDiheTypesParm,**diheTypes,**imprTypes;
    double **bondTypesParm,**angTypesParm,**ubTypesParm;
    double **diheTypesParm,**imprTypesParm,**nonBondedTypesParm;
} INPUTS;

typedef struct
{
    int nAtom,nDegFree,nFrozen;

    int nBond,nAngle,nDihedral,nImproper,nUb;

    int nConst;
} PARAM;

typedef struct
{
    int type;
    int a,b;
    double k,r0,beta;
} BOND;

typedef struct
{
    int type;
    int a,b,c;
    double k,theta0;
} ANGLE;

typedef struct
{
    int type,order;
    int a,b,c,d;
    double *k,*phi0,*mult;
} DIHE;

typedef struct
{
    int type;
    double eps,sig,bet;
    double eps14,sig14,bet14;
} VDW;

typedef struct
{
    int a,b;
    double r0;
} CONSTRAINT;

/** Formats defintition */

#define _FORMAT1_ "%d\t%s\t%d\t% #lf\t% #lf\t%d\n"
#define _FORMAT2_ "%d\t%d\t%d\t% #lf\t% #lf\t% #lf\n"
#define _FORMAT3_ "%d\t%d\t% #lf\n"
#define _FORMAT4_ "%d\t%d\t% #lf\t% #lf\t% #lf\n"
#define _FORMAT5_ "%d\t%d\t%d\t%d\t% #lf\t% #lf\n"
#define _FORMAT6_ "%d\t%d\t%d\t%d\t%d\t%d\t% #lf\t% #lf\t% #lf\n"
#define _FORMAT7_ "%d\t%d\t% #lf\t% #lf\t% #lf\t% #lf\t% #lf\t% #lf\n"

/** End of formats definition */

#define sq6rt2  (1.12246204830937)

/** Beginning functions prototypes */

void const_SELE(PARAM *param,BOND **bond,CONSTRAINT **constList);

void write_FORF(PARAM *param,ATOM atom[],BOND bond[],ANGLE angle[],
                DIHE dihe[],DIHE impr[],CONSTRAINT constList[],
                BOND ub[],VDW vdw[]);

void read_PSF(INPUTS *inp,PARAM *param,ATOM **atom,BOND **bond,
              ANGLE **angle,DIHE **dihe,DIHE **impr,CONSTRAINT **constList);

void read_TOP(INPUTS *inp);

void read_PAR(INPUTS *inp);

void read_CONF(PARAM *param,ATOM atom[]);

void setup(INPUTS *inp,PARAM *param,ATOM atom[],BOND bond[],
           ANGLE angle[],DIHE dihe[],DIHE impr[],CONSTRAINT constList[],
           BOND **ub,VDW **vdw);

void free_temp_array(INPUTS *inp);

void error(int errorNumber);

void** calloc_2D(int dim1, int dim2, size_t si);

void free_2D(int dim1, ...);

void*** calloc_3D(int dim1, int dim2, int dim3, size_t si);

void free_3D(int dim1, int dim2, ...);

/** End functions prototypes */

int main(int argc, char* argv[])
{

    /** Beginning of structures declaration. */

    INPUTS inp;

    PARAM param;

    ATOM *atom=NULL;

    BOND *bond=NULL,*ub=NULL;
    ANGLE *angle=NULL;
    DIHE *dihe=NULL,*impr=NULL;

    VDW *vdw=NULL;

    CONSTRAINT *constList=NULL;

    /** End of structures declarations. */

    read_PSF(&inp,&param,&atom,&bond,&angle,&dihe,&impr,&constList);

    read_TOP(&inp);

    read_PAR(&inp);

    read_CONF(&param,atom);

    const_SELE(&param,&bond,&constList);

    setup(&inp,&param,atom,bond,angle,dihe,impr,constList,&ub,&vdw);

    write_FORF(&param,atom,bond,angle,dihe,impr,constList,ub,vdw);

    free_temp_array(&inp);

    return 0;

}

void const_SELE(PARAM *param,BOND **bond,CONSTRAINT **constList)
{
    int i,j,k;
    int ia,ib;
    int keep=1,provide=0;
    double rc;

    printf("How many constraint bonds?\n");
    scanf("%d",&(param->nConst));

    if(param->nConst>0)
    {

        printf("Do you want to use CHARMM distance (0) or provide it (1)?\n");
        scanf("%d",&provide);

        *constList=(CONSTRAINT*)malloc(param->nConst*sizeof(CONSTRAINT));

        if(provide)
        {
            printf("Give the list of constrained pairs and the ref. distance:\n");

            for(i=0; i<param->nConst; i++)
            {
                scanf("%d %d %lf",&ia,&ib,&rc);
                (*constList)[i].a=ia;
                (*constList)[i].b=ib;
                (*constList)[i].r0=rc;
            }
        }
        else
        {
            printf("Give the list of constrained pairs:\n");

            for(i=0; i<param->nConst; i++)
            {
                (*constList)[i].r0=-1.0;
                scanf("%d %d",&ia,&ib);
                (*constList)[i].a=ia;
                (*constList)[i].b=ib;
            }
        }

        printf("ConstList done\n");

        BOND *temp;

        temp=(BOND*)malloc(param->nBond*sizeof(BOND));

        k=0;

        for(i=0; i<param->nBond; i++)
        {

            ia=(*bond)[i].a;
            ib=(*bond)[i].b;

            keep=1;

            for(j=0; j<param->nConst; j++)
            {

                if( (ia==(*constList)[j].a) && (ib==(*constList)[j].b)  )
                {
                    keep=0;
                    param->nBond--;
                    break;
                }
                else if( (ia==(*constList)[j].b) && (ib==(*constList)[j].a)  )
                {
                    keep=0;
                    param->nBond--;
                    break;
                }

            }

            if(keep)
            {
                temp[k].a=(*bond)[i].a;
                temp[k].b=(*bond)[i].b;

                k++;
            }

        }

        free(*bond);

        param->nBond=k;

        *bond=(BOND*)malloc(param->nBond*sizeof(BOND));

        for(i=0; i<param->nBond; i++)
        {
            (*bond)[i].a=temp[i].a;
            (*bond)[i].b=temp[i].b;
        }

        free(temp);

    }
}

void write_FORF(PARAM *param,ATOM atom[],BOND bond[],ANGLE angle[],
                DIHE dihe[],DIHE impr[],CONSTRAINT constList[],
                BOND ub[],VDW vdw[])
{
    FILE *forfFile=fopen("FORF","w");
    int i,j,nDihe,nImpr;

    fprintf(forfFile,"# Force Filed file for MDBas\n");

    fprintf(forfFile,"Atoms\t%d\n",param->nAtom);
    fprintf(forfFile,"# index\tlabel\ttype\tcharge\tmass\tfrozen\n");

    for(i=0; i<param->nAtom; i++)
        fprintf(forfFile,_FORMAT1_,i+1,atom[i].label,atom[i].type+1,atom[i].q,atom[i].m,atom[i].frozen);

    fprintf(forfFile,"Bonds\t%d\n",param->nBond);
    fprintf(forfFile,"# ia\tib\ttype\tk\tr_0\tbeta\n");

    for(i=0; i<param->nBond; i++)
        fprintf(forfFile,_FORMAT2_,bond[i].a+1,bond[i].b+1,bond[i].type,bond[i].k,bond[i].r0,bond[i].beta);

    fprintf(forfFile,"Constraints\t%d\n",param->nConst);
    fprintf(forfFile,"# ia\tib\tr_0\n");

    for(i=0; i<param->nConst; i++)
        fprintf(forfFile,_FORMAT3_,constList[i].a+1,constList[i].b+1,constList[i].r0);

    fprintf(forfFile,"Urey-Bradley\t%d\n",param->nUb);
    fprintf(forfFile,"# ia\tib\ttype\tk r_0\tbeta\n");

    for(i=0; i<param->nUb; i++)
        fprintf(forfFile,_FORMAT4_,ub[i].a+1,ub[i].b+1,ub[i].k,ub[i].r0);

    fprintf(forfFile,"Angles\t%d\n",param->nAngle);
    fprintf(forfFile,"# ia\tib\tic\ttype\tk\ttheta_0\n");

    for(i=0; i<param->nAngle; i++)
        fprintf(forfFile,_FORMAT5_,angle[i].a+1,angle[i].b+1,angle[i].c+1,angle[i].type,angle[i].k,angle[i].theta0);

    nDihe=0;
    for(i=0; i<param->nDihedral; i++)
    {
        for(j=0; j<dihe[i].order; j++)
        {
            nDihe++;
        }
    }

    fprintf(forfFile,"Dihedrals\t%d\n",nDihe);
    fprintf(forfFile,"# ia\tib\tic\tid\ttype\torder\tk\tphi_0\tmult\n");

    for(i=0; i<param->nDihedral; i++)
    {
        for(j=0; j<dihe[i].order; j++)
            fprintf(forfFile,_FORMAT6_,dihe[i].a+1,dihe[i].b+1,dihe[i].c+1,dihe[i].d+1,dihe[i].type,dihe[i].order,dihe[i].k[j],dihe[i].phi0[j],dihe[i].mult[j]);
    }

    nImpr=0;
    for(i=0; i<param->nImproper; i++)
    {
        for(j=0; j<impr[i].order; j++)
        {
            nImpr++;
        }
    }

    fprintf(forfFile,"Impropers\t%d\n",nImpr);
    fprintf(forfFile,"# ia\tib\tic\tid\ttype\torder\tk\tphi_0\tmult\n");

    for(i=0; i<param->nImproper; i++)
    {
        for(j=0; j<impr[i].order; j++)
            fprintf(forfFile,_FORMAT6_,impr[i].a+1,impr[i].b+1,impr[i].c+1,impr[i].d+1,impr[i].type,impr[i].order,impr[i].k[j],impr[i].phi0[j],impr[i].mult[j]);
    }

    fprintf(forfFile,"vdw\t%d\n",param->nAtom);
    fprintf(forfFile,"# index\ttype\teps\tsig\tbeta\teps14\tsig14\tbeta14\n");

    for(i=0; i<param->nAtom; i++)
        fprintf(forfFile,_FORMAT7_,i+1,vdw[i].type,vdw[i].eps,vdw[i].sig,vdw[i].bet,vdw[i].eps14,vdw[i].sig14,vdw[i].bet14);

    fprintf(forfFile,"end");

    fclose(forfFile);
}

void read_PSF(INPUTS *inp,PARAM *param,ATOM **atom,BOND **bond,
              ANGLE **angle,DIHE **dihe,DIHE **impr,CONSTRAINT **constList)
{
    FILE *psfFile=NULL;
    char buff1[1024]="", *buff2=NULL, *buff3=NULL, *buff4=NULL, *buff5=NULL;
    int i,j,k,kk,kt,nalloc=100,nincr=100;
    int ia,ib,ic,id;

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
            param->nAtom=atoi(buff2);
            break;
        }
    }

    inp->typesNum=(int*)malloc(nalloc*sizeof(*(inp->typesNum)));

    *atom=(ATOM*)malloc(param->nAtom*sizeof(**atom));

    k=-1;
    for(i=0; i<param->nAtom; i++)
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
            for(j=0; j<k; j++)
            {
                if(kt==inp->typesNum[j])
                {
                    k--;
                    kk=j;
                    break;
                }
            }

            (*atom)[i].type=kk;
            (*atom)[i].q=atof(buff4);
            (*atom)[i].m=atof(buff5);

            buff2=strtok(NULL," \n\t");

            (*atom)[i].frozen=atoi(buff2);

            if( ( (*atom)[i].frozen!=0 ) && ( (*atom)[i].frozen!=1 ) )
                error(27);

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
            param->nBond=atoi(buff2);
            break;
        }
    }

    *bond=(BOND*)malloc(param->nBond*sizeof(BOND));

    i=0;
    while(fgets(buff1,1024,psfFile)!=NULL)
    {
        buff2=strtok(buff1," \n\t");
        buff3=strtok(NULL," \n\t");

        while(buff2!=NULL && buff3!=NULL)
        {
            ia=atoi(buff2)-1;
            ib=atoi(buff3)-1;

            if( ((*atom)[ia].frozen*(*atom)[ib].frozen) ) /** Test if both atoms are Frozen. */
            {
                buff2=strtok(NULL," \n\t");
                buff3=strtok(NULL," \n\t");

                param->nBond--;
            }
            else /** Test if both atoms are Frozen. */
            {
                (*bond)[i].a=ia;
                (*bond)[i].b=ib;

                buff2=strtok(NULL," \n\t");
                buff3=strtok(NULL," \n\t");
                i++;
            } /** Test if both atoms are Frozen. */

        }
        if(buff2!=NULL && buff3==NULL)
            error(23);
        if(i==param->nBond)
            break;
    }

    while(fgets(buff1,1024,psfFile)!=NULL)
    {
        if(strstr(buff1,"NTHETA")!=NULL)
        {
            buff2=strtok(buff1," \n\t");
            param->nAngle=atoi(buff2);
            break;
        }
    }

    *angle=(ANGLE*)malloc(param->nAngle*sizeof(ANGLE));

    i=0;
    while(fgets(buff1,1024,psfFile)!=NULL)
    {
        buff2=strtok(buff1," \n\t");
        buff3=strtok(NULL," \n\t");
        buff4=strtok(NULL," \n\t");

        while(buff2!=NULL && buff3!=NULL && buff4!=NULL)
        {
            ia=atoi(buff2)-1;
            ib=atoi(buff3)-1;
            ic=atoi(buff4)-1;

            if( ((*atom)[ia].frozen*(*atom)[ib].frozen*(*atom)[ic].frozen) ) /** Test if the three atoms are Frozen. */
            {
                buff2=strtok(NULL," \n\t");
                buff3=strtok(NULL," \n\t");
                buff4=strtok(NULL," \n\t");

                param->nAngle--;
            }
            else /** Test if the three atoms are Frozen. */
            {
                (*angle)[i].a=ia;
                (*angle)[i].b=ib;
                (*angle)[i].c=ic;

                buff2=strtok(NULL," \n\t");
                buff3=strtok(NULL," \n\t");
                buff4=strtok(NULL," \n\t");

                i++;
            } /** Test if the three atoms are Frozen. */

        }
        if((buff2!=NULL && buff3==NULL)||(buff2!=NULL && buff4==NULL))
            error(24);
        if(i==param->nAngle)
            break;
    }

    while(fgets(buff1,1024,psfFile)!=NULL)
    {
        if(strstr(buff1,"NPHI")!=NULL)
        {
            buff2=strtok(buff1," \n\t");
            param->nDihedral=atoi(buff2);
            break;
        }
    }

    *dihe=(DIHE*)malloc(param->nDihedral*sizeof(DIHE));

    i=0;
    while(fgets(buff1,1024,psfFile)!=NULL)
    {
        buff2=strtok(buff1," \n\t");
        buff3=strtok(NULL," \n\t");
        buff4=strtok(NULL," \n\t");
        buff5=strtok(NULL," \n\t");

        while(buff2!=NULL && buff3!=NULL && buff4!=NULL)
        {
            ia=atoi(buff2)-1;
            ib=atoi(buff3)-1;
            ic=atoi(buff4)-1;
            id=atoi(buff5)-1;

            if( ((*atom)[ia].frozen*(*atom)[ib].frozen*(*atom)[ic].frozen*(*atom)[id].frozen) ) /** Test if the four atoms are Frozen. */
            {
                buff2=strtok(NULL," \n\t");
                buff3=strtok(NULL," \n\t");
                buff4=strtok(NULL," \n\t");
                buff5=strtok(NULL," \n\t");

                param->nDihedral--;
            }
            else /** Test if the four atoms are Frozen. */
            {
                (*dihe)[i].a=ia;
                (*dihe)[i].b=ib;
                (*dihe)[i].c=ic;
                (*dihe)[i].d=id;

                buff2=strtok(NULL," \n\t");
                buff3=strtok(NULL," \n\t");
                buff4=strtok(NULL," \n\t");
                buff5=strtok(NULL," \n\t");
                i++;
            } /** Test if the four atoms are Frozen. */

        }
        if((buff2!=NULL && buff3==NULL)||(buff2!=NULL && buff4==NULL)||(buff2!=NULL && buff5==NULL))
            error(25);
        if(i==param->nDihedral)
            break;
    }

    while(fgets(buff1,1024,psfFile)!=NULL)
    {
        if(strstr(buff1,"NIMPHI")!=NULL)
        {
            buff2=strtok(buff1," \n\t");
            param->nImproper=atoi(buff2);
            break;
        }
    }

    *impr=(DIHE*)malloc(param->nImproper*sizeof(DIHE));

    i=0;
    while(fgets(buff1,1024,psfFile)!=NULL)
    {
        buff2=strtok(buff1," \n\t");
        buff3=strtok(NULL," \n\t");
        buff4=strtok(NULL," \n\t");
        buff5=strtok(NULL," \n\t");

        while(buff2!=NULL && buff3!=NULL && buff4!=NULL)
        {
            ia=atoi(buff2)-1;
            ib=atoi(buff3)-1;
            ic=atoi(buff4)-1;
            id=atoi(buff5)-1;

            if( ((*atom)[ia].frozen*(*atom)[ib].frozen*(*atom)[ic].frozen*(*atom)[id].frozen) ) /** Test if the four atoms are Frozen. */
            {
                buff2=strtok(NULL," \n\t");
                buff3=strtok(NULL," \n\t");
                buff4=strtok(NULL," \n\t");
                buff5=strtok(NULL," \n\t");
            }
            else /** Test if the four atoms are Frozen. */
            {
                (*impr)[i].a=ia;
                (*impr)[i].b=ib;
                (*impr)[i].c=ic;
                (*impr)[i].d=id;

                buff2=strtok(NULL," \n\t");
                buff3=strtok(NULL," \n\t");
                buff4=strtok(NULL," \n\t");
                buff5=strtok(NULL," \n\t");
                i++;
            } /** Test if the four atoms are Frozen. */

        }
        if((buff2!=NULL && buff3==NULL)||(buff2!=NULL && buff4==NULL)||(buff2!=NULL && buff5==NULL))
            error(26);
        if(i==param->nImproper)
            break;
    }

    param->nFrozen=0;
    for(i=0; i<param->nAtom; i++)
    {
        param->nFrozen+=(*atom)[i].frozen;

        if( ((*atom)[i].frozen) || ((*atom)[i].m<DBL_EPSILON) )
        {
            (*atom)[i].m=0.;
            (*atom)[i].rm=0.;
        }
        else
        {
            (*atom)[i].rm=1./(*atom)[i].m;
        }

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
    for(i=0; i<inp->nTypes; i++)
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
                for(i=0; i<inp->nTypes; i++)
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

    for(i=0; i<inp->nBondTypes; i++)
    {
        inp->bondTypes[i]=(int*)malloc(2*sizeof(**(inp->bondTypes)));
        inp->bondTypesParm[i]=(double*)malloc(2*sizeof(**(inp->bondTypesParm)));
    }

//       Angle types arrays allocation

    inp->angTypes=(int**)malloc(inp->nAngTypes*sizeof(*(inp->angTypes)));
    inp->angTypesParm=(double**)malloc(inp->nAngTypes*sizeof(*(inp->angTypesParm)));

    for(i=0; i<inp->nAngTypes; i++)
    {
        inp->angTypes[i]=(int*)malloc(3*sizeof(**(inp->angTypes)));
        inp->angTypesParm[i]=(double*)malloc(2*sizeof(**(inp->angTypesParm)));
    }

//       Uray-Bradley types arrays allocation

    inp->ubTypes=(int**)malloc(inp->nUbTypes*sizeof(*(inp->ubTypes)));
    inp->ubTypesParm=(double**)malloc(inp->nUbTypes*sizeof(*(inp->ubTypesParm)));

    for(i=0; i<inp->nUbTypes; i++)
    {
        inp->ubTypes[i]=(int*)malloc(3*sizeof(**(inp->ubTypes)));
        inp->ubTypesParm[i]=(double*)malloc(2*sizeof(**(inp->ubTypesParm)));
    }

//       Dihedral types arrays allocation

    inp->diheTypes=(int**)malloc(inp->nDiheTypes*sizeof(*(inp->diheTypes)));
    inp->diheTypesParm=(double**)malloc(inp->nDiheTypes*sizeof(*(inp->diheTypesParm)));
    inp->nDiheTypesParm=(int*)malloc(inp->nDiheTypes*sizeof(*(inp->nDiheTypesParm)));

    for(i=0; i<inp->nDiheTypes; i++)
    {
        inp->diheTypes[i]=(int*)malloc(4*sizeof(**(inp->diheTypes)));
        inp->diheTypesParm[i]=(double*)malloc(3*sizeof(**(inp->diheTypesParm)));
    }

    for(i=0; i<inp->nDiheTypes; i++)
        inp->nDiheTypesParm[i]=1;

//       Improper types arrays allocation

    inp->imprTypes=(int**)malloc(inp->nImprTypes*sizeof(*(inp->imprTypes)));
    inp->imprTypesParm=(double**)malloc(inp->nImprTypes*sizeof(*(inp->imprTypesParm)));

    for(i=0; i<inp->nImprTypes; i++)
    {
        inp->imprTypes[i]=(int*)malloc(4*sizeof(**(inp->imprTypes)));
        inp->imprTypesParm[i]=(double*)malloc(3*sizeof(**(inp->imprTypesParm)));
    }

//       Non-bonding types arrays allocation

    inp->nonBondedTypesParm=(double**)malloc(inp->nTypes*sizeof(*(inp->nonBondedTypesParm)));

    for(i=0; i<inp->nTypes; i++)
        inp->nonBondedTypesParm[i]=(double*)malloc(6*sizeof(**(inp->nonBondedTypesParm)));

    for(i=0; i<inp->nTypes; i++)
    {
        for(j=0; j<6; j++)
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
                for(l=0; l<inp->nTypes; l++)
                {
                    if(!strcmp(buff2,inp->types[l]))
                    {
                        ia=l;
                        break;
                    }
                }

                if(ia==-1)
                    continue;

                for(l=0; l<inp->nTypes; l++)
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

                for(l=0; l<inp->nTypes; l++)
                {
                    if(!strcmp(buff2,inp->types[l]))
                    {
                        ia=l;
                        break;
                    }
                }

                if(ia==-1)
                    continue;

                for(l=0; l<inp->nTypes; l++)
                {
                    if(!strcmp(buff3,inp->types[l]))
                    {
                        ib=l;
                        break;
                    }
                }

                if(ib==-1)
                    continue;

                for(l=0; l<inp->nTypes; l++)
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
                    for(l=0; l<inp->nTypes; l++)
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

                for(l=0; l<inp->nTypes; l++)
                {
                    if(!strcmp(buff3,inp->types[l]))
                    {
                        ib=l;
                        break;
                    }
                }

                if(ib==-1)
                    continue;

                for(l=0; l<inp->nTypes; l++)
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
                    for(l=0; l<inp->nTypes; l++)
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
                for(l=0; l<inp->nDiheTypes; l++)
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
                    for(l=0; l<inp->nTypes; l++)
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
                    for(l=0; l<inp->nTypes; l++)
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
                    for(l=0; l<inp->nTypes; l++)
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
                    for(l=0; l<inp->nTypes; l++)
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
                for(l=0; l<inp->nTypes; l++)
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

void read_CONF(PARAM *param,ATOM atom[])
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

    if(natomCheck!=param->nAtom)
        error(41);

    for(i=0; i<natomCheck; i++)
    {
        fscanf(confFile,"%d %d %s %s %lf %lf %lf %s %d %lf",&atn,&ire,ren,atl,&xx,&yy,&zz,sen,&res,&wei);
        atom[i].x=xx;
        atom[i].y=yy;
        atom[i].z=zz;

        atom[i].ires=ire;

    }

    fclose(confFile);
}

void setup(INPUTS *inp,PARAM *param,ATOM atom[],BOND bond[],
           ANGLE angle[],DIHE dihe[],DIHE impr[],CONSTRAINT constList[],
           BOND **ub,VDW **vdw)
{
    int i,j,k,ia,ib,ic,id,i0,i1,i2,i3,itype;

    if(param->nBond>0)
    {
        for(i=0; i<param->nBond; i++)
            bond[i].type=-1;
    }

    if(param->nAngle>0)
    {
        for(i=0; i<param->nAngle; i++)
            angle[i].type=-1;
    }

    if(param->nDihedral>0)
    {
        for(i=0; i<param->nDihedral; i++)
            dihe[i].type=-1;
    }

    if(param->nImproper>0)
    {
        for(i=0; i<param->nImproper; i++)
            impr[i].type=-1;
    }

    *vdw=(VDW*)malloc(param->nAtom*sizeof(VDW));

    for(i=0; i<param->nBond; i++)
    {
        ia=atom[bond[i].a].type;
        ib=atom[bond[i].b].type;

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
        for(j=0; j<inp->nBondTypes; j++)
        {
            if( (i0==inp->bondTypes[j][0]) && (i1==inp->bondTypes[j][1]) )
            {
                itype=j;
                break;
            }
        }

        if(itype==-1)
            error(71);

        bond[i].k=inp->bondTypesParm[itype][0]*2.;
        bond[i].r0=inp->bondTypesParm[itype][1];
        bond[i].type=0;
        bond[i].beta=0.;

    }

    if(param->nConst>0)
    {
        for(i=0; i<param->nConst; i++)
        {
            if(constList[i].r0<0)
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
                for(j=0; j<inp->nBondTypes; j++)
                {
                    if( (i0==inp->bondTypes[j][0]) && (i1==inp->bondTypes[j][1]) )
                    {
                        itype=j;
                        break;
                    }
                }

                if(itype==-1)
                    error(71);

                constList[i].r0=inp->bondTypesParm[itype][1];
            }
        }
    }

    param->nUb=0;
    for(i=0; i<param->nAngle; i++)
    {
        ia=atom[angle[i].a].type;
        ib=atom[angle[i].b].type;
        ic=atom[angle[i].c].type;

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
        for(j=0; j<inp->nAngTypes; j++)
        {
            if( (i0==inp->angTypes[j][0]) && (i1==inp->angTypes[j][1]) && (i2==inp->angTypes[j][2]) )
            {
                itype=j;
                break;
            }
        }

        if(itype==-1)
            error(72);

        angle[i].k=inp->angTypesParm[itype][0]*2.;
        angle[i].theta0=inp->angTypesParm[itype][1];
        angle[i].type=0;

        for(j=0; j<inp->nUbTypes; j++)
        {
            if( (i0==inp->ubTypes[j][0]) && (i1==inp->ubTypes[j][1]) && (i2==inp->ubTypes[j][2]) )
            {
                param->nUb++;
                break;
            }
        }
    }

    if(param->nUb>0)
    {
        *ub=(BOND*)malloc(param->nUb*sizeof(BOND));

        for(i=0; i<param->nUb; i++)
            (*ub)[i].type=-1;
    }

    k=0;
    for(i=0; i<param->nAngle; i++)
    {
        ia=atom[angle[i].a].type;
        ib=atom[angle[i].b].type;
        ic=atom[angle[i].c].type;

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
        for(j=0; j<inp->nUbTypes; j++)
        {
            if( (i0==inp->ubTypes[j][0]) && (i1==inp->ubTypes[j][1]) && (i2==inp->ubTypes[j][2]) )
            {
                itype=j;
                break;
            }
        }

        if(itype==-1)
            continue;

        (*ub)[k].a=angle[i].a;
        (*ub)[k].b=angle[i].c;
        (*ub)[k].k=inp->ubTypesParm[itype][0]*2.;
        (*ub)[k].r0=inp->ubTypesParm[itype][1];
        (*ub)[k].type=0;
        (*ub)[k].beta=0.;
        k++;
    }

    for(i=0; i<param->nDihedral; i++)
    {
        ia=atom[dihe[i].a].type;
        ib=atom[dihe[i].b].type;
        ic=atom[dihe[i].c].type;
        id=atom[dihe[i].d].type;

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
        for(j=0; j<inp->nDiheTypes; j++)
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

            for(j=0; j<inp->nDiheTypes; j++)
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

        dihe[i].order=inp->nDiheTypesParm[itype];

        dihe[i].k=(double*)malloc(inp->nDiheTypesParm[itype]*sizeof(double));
        dihe[i].mult=(double*)malloc(inp->nDiheTypesParm[itype]*sizeof(double));
        dihe[i].phi0=(double*)malloc(inp->nDiheTypesParm[itype]*sizeof(double));

        for(j=0; j<dihe[i].order; j++)
        {
            k=3*j;
            dihe[i].k[j]=inp->diheTypesParm[itype][k];
            dihe[i].mult[j]=inp->diheTypesParm[itype][k+1];
            dihe[i].phi0[j]=inp->diheTypesParm[itype][k+2];
        }

        if(dihe[i].order==1&&dihe[i].mult[0]<0.5)
        {
            dihe[i].k[0]=dihe[i].k[0]*2.;
            dihe[i].mult[0]=0.;
            dihe[i].type=2;
        }
        else if(dihe[i].order>1&&dihe[i].mult[0]<0.5)
            error(50);
        else
            dihe[i].type=1;
    }

    for(i=0; i<param->nImproper; i++)
    {
        ia=atom[impr[i].a].type;
        ib=atom[impr[i].b].type;
        ic=atom[impr[i].c].type;
        id=atom[impr[i].d].type;

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
        for(j=0; j<inp->nImprTypes; j++)
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

            for(j=0; j<inp->nImprTypes; j++)
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

            for(j=0; j<inp->nImprTypes; j++)
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

            for(j=0; j<inp->nImprTypes; j++)
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

            for(j=0; j<inp->nImprTypes; j++)
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

            for(j=0; j<inp->nImprTypes; j++)
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

            for(j=0; j<inp->nImprTypes; j++)
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

        impr[i].order=1;

        impr[i].k=(double*)malloc(sizeof(double));
        impr[i].mult=(double*)malloc(sizeof(double));
        impr[i].phi0=(double*)malloc(sizeof(double));

        impr[i].k[0]=inp->imprTypesParm[itype][0];
        impr[i].mult[0]=inp->imprTypesParm[itype][1];
        impr[i].phi0[0]=inp->imprTypesParm[itype][2];

        if(impr[i].mult[0]<0.5)
        {
            impr[i].k[0]=impr[i].k[0]*2.;
            impr[i].mult[0]=0.;
            impr[i].type=2;
        }
        else
            impr[i].type=1;
    }

    for(i=0; i<param->nAtom; i++)
    {

        (*vdw)[i].eps=-inp->nonBondedTypesParm[atom[i].type][1];
        (*vdw)[i].sig=inp->nonBondedTypesParm[atom[i].type][2]/sq6rt2;
        (*vdw)[i].bet=inp->nonBondedTypesParm[atom[i].type][0];

        if(inp->nonBondedTypesParm[atom[i].type][5]>=0.)
        {
            (*vdw)[i].eps14=-inp->nonBondedTypesParm[atom[i].type][4];
            (*vdw)[i].sig14=inp->nonBondedTypesParm[atom[i].type][5]/sq6rt2;
            (*vdw)[i].bet14=inp->nonBondedTypesParm[atom[i].type][3];
        }
        else
        {
            (*vdw)[i].eps14=(*vdw)[i].eps;
            (*vdw)[i].sig14=(*vdw)[i].sig;
            (*vdw)[i].bet14=(*vdw)[i].bet;
        }
    }

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
    case 47:
        printf("MDBas cannot open configuration file RESCONF.\n");
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
    default:
        printf("MDBas failed due to unknown error number: %d\n",errorNumber);
        printf("Reading the manual will not help you. You are by yourself.\n");
        printf("Errare humanum est.\n");
        break;
    }

    exit(errorNumber);

}

//for allocating an array of dimensions dim1*dim2 and of bytes size si
void** calloc_2D(int dim1, int dim2, size_t si)
{
    int i;
    void **array=NULL;

    //Here, allocation of the first dimension : rows
    array = calloc(dim1,si);
    assert(array!=NULL);

    //Then allocation of the second dimensions : columns
    for (i=0; i<dim1; i++)
    {
        array[i] = calloc(dim2,si);
        assert(array[i]!=NULL);
    }

    return array;
}

//for freeing one or more dynamically allocated 2D arrays of first dimension dim1
void free_2D(int dim1, ...)
{
    int i;
    void **array=NULL;
    va_list ap;

    va_start(ap,dim1);
    array=va_arg(ap,void**);

    while(array!=NULL)
    {
        for(i=0; i<dim1; i++)
            free(array[i]);
        free(array);
        array=NULL;
        array=va_arg(ap,void**);
    }
    va_end(ap);
}

//for allocating an array of dimensions dim1*dim2*dim3 and of bytes size si
void*** calloc_3D(int dim1, int dim2, int dim3, size_t si)
{
    int i,j;
    void ***array=NULL;

    array=calloc(dim1,si);
    assert(array!=NULL);

    for(i=0; i<dim1; i++)
    {
        array[i]=calloc(dim2,si);
        assert(array[i]!=NULL);

        for(j=0; j<dim2; j++)
        {
            array[i][j]=calloc(dim3,si);
            assert(array[i]!=NULL);
        }
    }
    return array;
}

//for freeing one or more dynamically allocated 3D arrays of dimensions dim1,dim2
void free_3D(int dim1, int dim2, ...)
{
    int i,j;
    void ***array=NULL;
    va_list ap;

    va_start(ap,dim2);
    array=va_arg(ap,void***);

    while(array!=NULL)
    {
        for(i=0; i<dim1; i++)
        {
            for(j=0; j<dim2; j++)
                free(array[i][j]);
            free(array[i]);
        }
        free(array);
        array=NULL;
        array=va_arg(ap,void***);
    }
    va_end(ap);
}
