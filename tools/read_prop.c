#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define _FORMATDS_ "%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\n"
#define _FORMATLS_ "# %11s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\n"

#define _FORMATDA_ "%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\t%#13.5le\n"
#define _FORMATLA_ "# %11s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\n"

enum OPT {
    SHORT,
    ALL,
};

int main(int argc, char* argv[])
{

    FILE *propFile=NULL;

    size_t count=29;
    size_t test;

    int i;

    double buffer[29]= {0.};

    enum OPT opt;

    char propName[512];

    char label[29][7]= { {"Step"} , {"Time"} , {"Temp"} , {"Press"} , {"Vol"} ,
        {"Etot"} , {"Ekin"}  , {"Epot"} , {"Ecoul"} , {"Evdw"} ,
        {"Ebond"} , {"Eangle"} , {"Eub"} , {"Edihe"} , {"Eimpr"} ,
        {"Virtot"} , {"Virpot"} , {"Virele"} , {"Virvdw"} ,
        {"Virbnd"} , {"Virub"} , {"Virshk"} , {"Consv"} ,
        {"box1"} , {"box2"} , {"box3"} , {"box4"} ,
        {"box5"} , {"box6"}
    };

    strcpy(propName,"PROP");
    opt=SHORT;

    i=1;
    while(i<argc)
    {

        if(!strcmp(argv[i],"-i"))
        {
            strcpy(propName,argv[++i]);
        }
        else if(!strcmp(argv[i],"-a"))
        {
            opt=ALL;
        }
        else
        {
            printf("%s prop_file\n",argv[0]);
            exit(0);
        }

        i++;
    }

    propFile=fopen(propName,"rb");

    switch(opt)
    {
    case SHORT:

        printf(_FORMATLS_,label[0],label[1],label[2],label[3],label[4],
               label[5],label[6],label[7],label[8],label[9],
               label[10],label[11],label[12],label[13],label[14]);

        while(!feof(propFile))
        {

            test=fread(buffer,sizeof(double),29,propFile);

            if(test!=count)
            {
                printf("Incomplete read statement.\n");
                exit(0);
            }

            printf(_FORMATDS_,buffer[0],buffer[1],buffer[2],buffer[3],buffer[4],
                   buffer[5],buffer[6],buffer[7],buffer[8],buffer[9],
                   buffer[10],buffer[11],buffer[12],buffer[13],buffer[14]);

        }
        break;

    case ALL:

        printf(_FORMATLA_,label[0],label[1],label[2],label[3],label[4],
               label[5],label[6],label[7],label[8],label[9],
               label[10],label[11],label[12],label[13],label[14],
               label[15],label[16],label[17],label[18],label[19],
               label[20],label[21],label[22],label[23],label[24],
               label[25],label[26],label[27],label[28]);

        while(!feof(propFile))
        {

            test=fread(buffer,sizeof(double),29,propFile);

            if(test!=count)
            {
                printf("Incomplete read statement.\n");
                exit(0);
            }

            printf(_FORMATDA_,buffer[0],buffer[1],buffer[2],buffer[3],buffer[4],
                   buffer[5],buffer[6],buffer[7],buffer[8],buffer[9],
                   buffer[10],buffer[11],buffer[12],buffer[13],buffer[14],
                   buffer[15],buffer[16],buffer[17],buffer[18],buffer[19],
                   buffer[20],buffer[21],buffer[22],buffer[23],buffer[24],
                   buffer[25],buffer[26],buffer[27],buffer[28]);

        }
        break;
    }

    fclose(propFile);

    return 0;
}
