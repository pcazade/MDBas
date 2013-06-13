/**
 * \file list.c
 * \brief Contains functions for building and updating the nonbonded interactions list.
 * \author Pierre-Andre Cazade and Florent Hedin
 * \version alpha-branch
 * \date 2012
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "global.h"
#include "utils.h"
#include "memory.h"
#include "io.h"
#include "list.h"
#include "errors.h"

#if (defined TIMING && defined __unix__ && !defined __STRICT_ANSI__)
#define TIMER
#include "timing.h"
#endif


/** Pointer to the output file. **/
extern FILE *outFile;

void makelist(CTRL *ctrl, PARAM *param, PBC *box, NEIGH *neigh, CONSTRAINT constList[],
        BOND bond[], ANGLE angle[], DIHE dihe[], DIHE impr[], double x[], double y[],
        double z[], int frozen[], int **neighList, int **neighPair, int **neighOrder,
        int **neighList14, int ***exclList, int **exclPair)
{
    int nlcx, nlcy, nlcz;
    double cutnb;

    if (ctrl->newjob == 1)
    {

        neigh->linkRatio = MIN(MAX(neigh->linkRatio, 1), 5);

        cutnb = param->cutOff + param->delr;

        ctrl->keyLink = 1;

        nlcx = (int) (box->pa * (double) neigh->linkRatio / cutnb);
        nlcy = (int) (box->pb * (double) neigh->linkRatio / cutnb);
        nlcz = (int) (box->pc * (double) neigh->linkRatio / cutnb);

        if ((nlcx < 3) || (nlcy < 3) || (nlcz < 3))
            ctrl->keyLink = 0;

        if ((nlcx * nlcy * nlcz) <= 27)
            ctrl->keyLink = 0;

        if (ctrl->noLink)
            ctrl->keyLink = 0;

        if (ctrl->keyLink)
        {
#ifdef TIMER
            create_new_timer(TIMER_LNKCEL_BUILD);
            create_new_timer(TIMER_LNKCEL_UPDATE);
#endif
            link_cell_exclude_list(ctrl, param, neigh, constList, bond, angle, dihe, impr,
                    neighList14, exclList, exclPair);

            link_cell_verlet_list(param, box, neigh, x, y, z, frozen, neighList, neighPair,
                    neighOrder, *exclList, *exclPair);
        }
        else
        {
#ifdef TIMER
            create_new_timer(TIMER_VERLET_BUILD);
            create_new_timer(TIMER_VERLET_UPDATE);
#endif
            exclude_list(ctrl, param, neigh, constList, bond, angle, dihe, impr,
                    neighList14, exclList, exclPair);

            verlet_list(param, box, neigh, x, y, z, frozen, neighList, neighPair,
                    neighOrder, *exclList, *exclPair);
        }

        ctrl->newjob = 0;

    }
    else if ((param->step % neigh->update == 0))
    {

        if (ctrl->keyLink)
        {
            //       fprintf(outFile,"\nCalling link_cell_verlet_list_update at step %d\n",param->step);
            link_cell_verlet_list_update(param, box, neigh, x, y, z, frozen, neighList, *neighPair,
                    *neighOrder, *exclList, *exclPair);
            //       fprintf(outFile,"List update done.\n\n");
        }
        else
        {
            //       fprintf(outFile,"\nCalling verlet_list_update at step %d\n",param->step);
            verlet_list_update(param, box, neigh, x, y, z, frozen, neighList, *neighPair,
                    *exclList, *exclPair);
            //       fprintf(outFile,"List update done.\n\n");
        }

    }
}

void exclude_list(CTRL *ctrl, PARAM *param, NEIGH *neigh, CONSTRAINT constList[],
        BOND bond[], ANGLE angle[], DIHE dihe[], DIHE impr[],
        int **neighList14, int ***exclList, int **exclPair)
{
    int i, j, k, l, ii, jj, kk, ia, ib, ic, id, exclude;
    int **tempAtom = NULL, **tempVer14 = NULL, **tempConnect = NULL, *tempConnectNum = NULL;
    int nAlloc = 16, nIncr = 16, nConnect = 16;

    *exclPair = (int*) my_malloc(param->nAtom * sizeof (**exclPair));
    for (i = 0; i < param->nAtom; i++)
        (*exclPair)[i] = 0;

    tempAtom = (int**) my_malloc(param->nAtom * sizeof (*(tempAtom)));
    for (i = 0; i < param->nAtom; i++)
        tempAtom[i] = (int*) my_malloc(nAlloc * sizeof (**(tempAtom)));

    /** The following temporary arrays are designed for chechking the
     * connectivity of system and make sure that all existing but not
     * parameterised angles and dihedrals are properly excluded. */

    tempConnectNum = (int*) my_malloc(param->nAtom * sizeof (*tempConnectNum));
    for (i = 0; i < param->nAtom; i++)
        tempConnectNum[i] = 0;

    tempConnect = (int**) my_malloc(param->nAtom * sizeof (*tempConnect));
    for (i = 0; i < param->nAtom; i++)
        tempConnect[i] = (int*) my_malloc(nConnect * sizeof (**tempConnect));


    for (i = 0; i < param->nBond; i++)
    {
        ia = bond[i].a;
        ib = bond[i].b;

        if (ib < ia)
        {
            ii = ib;
            jj = ia;
        }
        else
        {
            ii = ia;
            jj = ib;
        }

        if ((*exclPair)[ii] >= nAlloc)
        {
            nAlloc += nIncr;
            for (j = 0; j < param->nAtom; j++)
                tempAtom[j] = (int*) realloc(tempAtom[j], nAlloc * sizeof (**(tempAtom)));
        }

        if (tempConnectNum[ii] >= nConnect || tempConnectNum[jj] >= nConnect)
        {
            nConnect += nIncr;
            for (j = 0; j < param->nAtom; j++)
                tempConnect[j] = (int*) realloc(tempConnect[j], nConnect * sizeof (**tempConnect));
        }

        tempConnect[ii][tempConnectNum[ii]] = jj;
        tempConnect[jj][tempConnectNum[jj]] = ii;
        tempConnectNum[ii]++;
        tempConnectNum[jj]++;

        tempAtom[ii][(*exclPair)[ii]] = jj;
        (*exclPair)[ii]++;
    }

    if (ctrl->keyConstH)
    {
        for (i = 0; i < param->nConst; i++)
        {
            ia = constList[i].a;
            ib = constList[i].b;

            if (ib < ia)
            {
                ii = ib;
                jj = ia;
            }
            else
            {
                ii = ia;
                jj = ib;
            }

            if ((*exclPair)[ii] >= nAlloc)
            {
                nAlloc += nIncr;
                for (j = 0; j < param->nAtom; j++)
                    tempAtom[j] = (int*) realloc(tempAtom[j], nAlloc * sizeof (**(tempAtom)));
            }

            if (tempConnectNum[ii] >= nConnect || tempConnectNum[jj] >= nConnect)
            {
                nConnect += nIncr;
                for (j = 0; j < param->nAtom; j++)
                    tempConnect[j] = (int*) realloc(tempConnect[j], nConnect * sizeof (**tempConnect));
            }

            tempConnect[ii][tempConnectNum[ii]] = jj;
            tempConnect[jj][tempConnectNum[jj]] = ii;
            tempConnectNum[ii]++;
            tempConnectNum[jj]++;

            tempAtom[ii][(*exclPair)[ii]] = jj;
            (*exclPair)[ii]++;
        }
    }

    for (i = 0; i < param->nAngle; i++)
    {
        ia = angle[i].a;
        ib = angle[i].b;
        ic = angle[i].c;

        if (ib < ia)
        {
            ii = ib;
            jj = ia;
        }
        else
        {
            ii = ia;
            jj = ib;
        }

        if ((*exclPair)[ii] >= nAlloc)
        {
            nAlloc += nIncr;
            for (j = 0; j < param->nAtom; j++)
                tempAtom[j] = (int*) realloc(tempAtom[j], nAlloc * sizeof (**(tempAtom)));
        }

        exclude = 1;
        for (j = 0; j < (*exclPair)[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
        {
            tempAtom[ii][(*exclPair)[ii]] = jj;
            (*exclPair)[ii]++;
        }

        if (ic < ia)
        {
            ii = ic;
            jj = ia;
        }
        else
        {
            ii = ia;
            jj = ic;
        }

        if ((*exclPair)[ii] >= nAlloc)
        {
            nAlloc += nIncr;
            for (j = 0; j < param->nAtom; j++)
                tempAtom[j] = (int*) realloc(tempAtom[j], nAlloc * sizeof (**(tempAtom)));
        }

        exclude = 1;
        for (j = 0; j < (*exclPair)[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
        {
            tempAtom[ii][(*exclPair)[ii]] = jj;
            (*exclPair)[ii]++;
        }

        if (ic < ib)
        {
            ii = ic;
            jj = ib;
        }
        else
        {
            ii = ib;
            jj = ic;
        }

        if ((*exclPair)[ii] >= nAlloc)
        {
            nAlloc += nIncr;
            for (j = 0; j < param->nAtom; j++)
                tempAtom[j] = (int*) realloc(tempAtom[j], nAlloc * sizeof (**(tempAtom)));
        }

        exclude = 1;
        for (j = 0; j < (*exclPair)[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
        {
            tempAtom[ii][(*exclPair)[ii]] = jj;
            (*exclPair)[ii]++;
        }

    }

    tempVer14 = (int**) my_malloc(5 * param->nDihedral * sizeof (*tempVer14));
    for (i = 0; i < 5 * param->nDihedral; i++)
        tempVer14[i] = (int*) my_malloc(2 * sizeof (**tempVer14));

    neigh->nPair14 = 0;

    for (i = 0; i < param->nDihedral; i++)
    {
        ia = dihe[i].a;
        ib = dihe[i].b;
        ic = dihe[i].c;
        id = dihe[i].d;

        if (ib < ia)
        {
            ii = ib;
            jj = ia;
        }
        else
        {
            ii = ia;
            jj = ib;
        }

        if ((*exclPair)[ii] >= nAlloc)
        {
            nAlloc += nIncr;
            for (j = 0; j < param->nAtom; j++)
                tempAtom[j] = (int*) realloc(tempAtom[j], nAlloc * sizeof (**(tempAtom)));
        }

        exclude = 1;
        for (j = 0; j < (*exclPair)[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
        {
            tempAtom[ii][(*exclPair)[ii]] = jj;
            (*exclPair)[ii]++;
        }

        if (ic < ia)
        {
            ii = ic;
            jj = ia;
        }
        else
        {
            ii = ia;
            jj = ic;
        }

        if ((*exclPair)[ii] >= nAlloc)
        {
            nAlloc += nIncr;
            for (j = 0; j < param->nAtom; j++)
                tempAtom[j] = (int*) realloc(tempAtom[j], nAlloc * sizeof (**(tempAtom)));
        }

        exclude = 1;
        for (j = 0; j < (*exclPair)[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
        {
            tempAtom[ii][(*exclPair)[ii]] = jj;
            (*exclPair)[ii]++;
        }

        if (id < ia)
        {
            ii = id;
            jj = ia;
        }
        else
        {
            ii = ia;
            jj = id;
        }

        if ((*exclPair)[ii] >= nAlloc)
        {
            nAlloc += nIncr;
            for (j = 0; j < param->nAtom; j++)
                tempAtom[j] = (int*) realloc(tempAtom[j], nAlloc * sizeof (**(tempAtom)));
        }

        exclude = 1;
        for (j = 0; j < (*exclPair)[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
        {
            tempVer14[neigh->nPair14][0] = ia;
            tempVer14[neigh->nPair14][1] = id;
            neigh->nPair14++;

            tempAtom[ii][(*exclPair)[ii]] = jj;
            (*exclPair)[ii]++;
        }

        if (ic < ib)
        {
            ii = ic;
            jj = ib;
        }
        else
        {
            ii = ib;
            jj = ic;
        }

        if ((*exclPair)[ii] >= nAlloc)
        {
            nAlloc += nIncr;
            for (j = 0; j < param->nAtom; j++)
                tempAtom[j] = (int*) realloc(tempAtom[j], nAlloc * sizeof (**(tempAtom)));
        }

        exclude = 1;
        for (j = 0; j < (*exclPair)[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
        {
            tempAtom[ii][(*exclPair)[ii]] = jj;
            (*exclPair)[ii]++;
        }

        if (id < ib)
        {
            ii = id;
            jj = ib;
        }
        else
        {
            ii = ib;
            jj = id;
        }

        if ((*exclPair)[ii] >= nAlloc)
        {
            nAlloc += nIncr;
            for (j = 0; j < param->nAtom; j++)
                tempAtom[j] = (int*) realloc(tempAtom[j], nAlloc * sizeof (**(tempAtom)));
        }

        exclude = 1;
        for (j = 0; j < (*exclPair)[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
        {
            tempAtom[ii][(*exclPair)[ii]] = jj;
            (*exclPair)[ii]++;
        }

        if (id < ic)
        {
            ii = id;
            jj = ic;
        }
        else
        {
            ii = ic;
            jj = id;
        }

        if ((*exclPair)[ii] >= nAlloc)
        {
            nAlloc += nIncr;
            for (j = 0; j < param->nAtom; j++)
                tempAtom[j] = (int*) realloc(tempAtom[j], nAlloc * sizeof (**(tempAtom)));
        }

        exclude = 1;
        for (j = 0; j < (*exclPair)[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
        {
            tempAtom[ii][(*exclPair)[ii]] = jj;
            (*exclPair)[ii]++;
        }

    }

    for (i = 0; i < param->nImproper; i++)
    {
        ia = impr[i].a;
        ib = impr[i].b;
        ic = impr[i].c;
        id = impr[i].d;

        if (ib < ia)
        {
            ii = ib;
            jj = ia;
        }
        else
        {
            ii = ia;
            jj = ib;
        }

        if ((*exclPair)[ii] >= nAlloc)
        {
            nAlloc += nIncr;
            for (j = 0; j < param->nAtom; j++)
                tempAtom[j] = (int*) realloc(tempAtom[j], nAlloc * sizeof (**(tempAtom)));
        }

        exclude = 1;
        for (j = 0; j < (*exclPair)[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
        {
            tempAtom[ii][(*exclPair)[ii]] = jj;
            (*exclPair)[ii]++;
        }

        if (ic < ia)
        {
            ii = ic;
            jj = ia;
        }
        else
        {
            ii = ia;
            jj = ic;
        }

        if ((*exclPair)[ii] >= nAlloc)
        {
            nAlloc += nIncr;
            for (j = 0; j < param->nAtom; j++)
                tempAtom[j] = (int*) realloc(tempAtom[j], nAlloc * sizeof (**(tempAtom)));
        }

        exclude = 1;
        for (j = 0; j < (*exclPair)[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
        {
            tempAtom[ii][(*exclPair)[ii]] = jj;
            (*exclPair)[ii]++;
        }

        if (id < ia)
        {
            ii = id;
            jj = ia;
        }
        else
        {
            ii = ia;
            jj = id;
        }

        if ((*exclPair)[ii] >= nAlloc)
        {
            nAlloc += nIncr;
            for (j = 0; j < param->nAtom; j++)
                tempAtom[j] = (int*) realloc(tempAtom[j], nAlloc * sizeof (**(tempAtom)));
        }

        exclude = 1;
        for (j = 0; j < (*exclPair)[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
        {
            tempAtom[ii][(*exclPair)[ii]] = jj;
            (*exclPair)[ii]++;
        }

        if (ic < ib)
        {
            ii = ic;
            jj = ib;
        }
        else
        {
            ii = ib;
            jj = ic;
        }

        if ((*exclPair)[ii] >= nAlloc)
        {
            nAlloc += nIncr;
            for (j = 0; j < param->nAtom; j++)
                tempAtom[j] = (int*) realloc(tempAtom[j], nAlloc * sizeof (**(tempAtom)));
        }

        exclude = 1;
        for (j = 0; j < (*exclPair)[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
        {
            tempAtom[ii][(*exclPair)[ii]] = jj;
            (*exclPair)[ii]++;
        }

        if (id < ib)
        {
            ii = id;
            jj = ib;
        }
        else
        {
            ii = ib;
            jj = id;
        }

        if ((*exclPair)[ii] >= nAlloc)
        {
            nAlloc += nIncr;
            for (j = 0; j < param->nAtom; j++)
                tempAtom[j] = (int*) realloc(tempAtom[j], nAlloc * sizeof (**(tempAtom)));
        }

        exclude = 1;
        for (j = 0; j < (*exclPair)[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
        {
            tempAtom[ii][(*exclPair)[ii]] = jj;
            (*exclPair)[ii]++;
        }

        if (id < ic)
        {
            ii = id;
            jj = ic;
        }
        else
        {
            ii = ic;
            jj = id;
        }

        if ((*exclPair)[ii] >= nAlloc)
        {
            nAlloc += nIncr;
            for (j = 0; j < param->nAtom; j++)
                tempAtom[j] = (int*) realloc(tempAtom[j], nAlloc * sizeof (**(tempAtom)));
        }

        exclude = 1;
        for (j = 0; j < (*exclPair)[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
        {
            tempAtom[ii][(*exclPair)[ii]] = jj;
            (*exclPair)[ii]++;
        }

    }

    for (ia = 0; ia < param->nAtom; ia++)
    {
        for (j = 0; j < tempConnectNum[ia]; j++)
        {
            ib = tempConnect[ia][j];
            for (k = 0; k < tempConnectNum[ib]; k++)
            {
                ic = tempConnect[ib][k];

                if (ic == ia)
                    continue;

                if (ic < ia)
                {
                    ii = ic;
                    jj = ia;
                }
                else
                {
                    ii = ia;
                    jj = ic;
                }

                if ((*exclPair)[ii] >= nAlloc)
                {
                    nAlloc += nIncr;
                    for (kk = 0; kk < param->nAtom; kk++)
                        tempAtom[kk] = (int*) realloc(tempAtom[kk], nAlloc * sizeof (**(tempAtom)));
                }

                exclude = 1;
                for (kk = 0; kk < (*exclPair)[ii]; kk++)
                {
                    if (tempAtom[ii][kk] == jj)
                    {
                        exclude = 0;
                        break;
                    }
                }

                if (exclude)
                {
                    tempAtom[ii][(*exclPair)[ii]] = jj;
                    (*exclPair)[ii]++;
                }

                for (l = 0; l < tempConnectNum[ic]; l++)
                {
                    id = tempConnect[ic][l];

                    if (id == ia || id == ib)
                        continue;

                    if (id < ia)
                    {
                        ii = id;
                        jj = ia;
                    }
                    else
                    {
                        ii = ia;
                        jj = id;
                    }

                    if ((*exclPair)[ii] >= nAlloc)
                    {
                        nAlloc += nIncr;
                        for (kk = 0; kk < param->nAtom; kk++)
                            tempAtom[kk] = (int*) realloc(tempAtom[kk], nAlloc * sizeof (**(tempAtom)));
                    }

                    if (neigh->nPair14 >= 5 * param->nDihedral)
                        my_error(DIHE_NONPARAM_ERROR, __FILE__, __LINE__, 0);

                    exclude = 1;
                    for (kk = 0; kk < (*exclPair)[ii]; kk++)
                    {
                        if (tempAtom[ii][kk] == jj)
                        {
                            exclude = 0;
                            break;
                        }
                    }

                    if (exclude)
                    {

                        tempVer14[neigh->nPair14][0] = ia;
                        tempVer14[neigh->nPair14][1] = id;
                        neigh->nPair14++;

                        tempAtom[ii][(*exclPair)[ii]] = jj;
                        (*exclPair)[ii]++;
                    }

                }
            }
        }
    }

    *neighList14 = (int*) my_malloc(2 * neigh->nPair14 * sizeof (**neighList14));

    for (i = 0; i < neigh->nPair14; i++)
    {
        (*neighList14)[2 * i] = tempVer14[i][0];
        (*neighList14)[2 * i + 1] = tempVer14[i][1];
    }

    free_2D(5 * param->nDihedral, tempVer14, NULL);

    free_2D(param->nAtom, tempConnect, NULL);
    free(tempConnectNum);

    if ((*exclPair)[param->nAtom - 1] != 0)
        my_error(EXCLLIST_LASTATOM_ERROR, __FILE__, __LINE__, 0);

    /** Brode-Ahlrichs sorted exclude list **/

    *exclList = (int**) my_malloc((param->nAtom - 1) * sizeof (**exclList));

    for (i = 0; i < param->nAtom - 1; i++)
    {
        (*exclList)[i] = (int*) my_malloc((*exclPair)[i] * sizeof (***exclList));

        for (j = 0; j < (*exclPair)[i]; j++)
            (*exclList)[i][j] = tempAtom[i][j];
    }

    free_2D(param->nAtom, tempAtom, NULL);

}

void verlet_list(PARAM *param, PBC *box, NEIGH *neigh, double x[], double y[], double z[],
        int frozen[], int **neighList, int **neighPair, int **neighOrder,
        int **exclList, int exclPair[])
{
    int i, j, k, offset, exclude;
    double r2, cutnb2, delta[3];

    cutnb2 = param->cutOff + param->delr;
    cutnb2 = X2(cutnb2);

    neigh->sizeList = 0;

#ifdef TIMER
    update_timer_begin(TIMER_VERLET_BUILD, __func__);
#endif

    for (i = 0; i < param->nAtom - 1; i++)
    {
        for (j = i + 1; j < param->nAtom; j++)
        {

            delta[0] = x[j] - x[i];
            delta[1] = y[j] - y[i];
            delta[2] = z[j] - z[i];

            r2 = dist(box, delta);

            if (r2 <= cutnb2)
            {
                exclude = 0;

                if ((frozen[i] * frozen[j]))
                {
                    exclude = 1;
                }
                else
                {
                    for (k = 0; k < exclPair[i]; k++)
                    {
                        if (exclList[i][k] == j)
                        {
                            exclude = 1;
                            break;
                        }
                    }
                }

                if (!exclude)
                {

                    neigh->sizeList++;

                }
            }
        }
    }

    neigh->sizeList = (int) (neigh->sizeList * (1. + 2. * TOLLIST)) + 1;

    *neighPair = (int*) my_malloc(param->nAtom * sizeof (**neighPair));

    for (i = 0; i < param->nAtom; i++)
        (*neighPair)[i] = 0;

    *neighOrder = (int*) my_malloc(param->nAtom * sizeof (**neighOrder));
    for (i = 0; i < param->nAtom; i++)
        (*neighOrder)[i] = i;

    *neighList = (int*) my_malloc(neigh->sizeList * sizeof (**neighList));

    offset = 0;

    for (i = 0; i < param->nAtom - 1; i++)
    {

        for (j = i + 1; j < param->nAtom; j++)
        {

            delta[0] = x[j] - x[i];
            delta[1] = y[j] - y[i];
            delta[2] = z[j] - z[i];

            r2 = dist(box, delta);

            if (r2 <= cutnb2)
            {
                exclude = 0;

                if ((frozen[i] * frozen[j]))
                {
                    exclude = 1;
                }
                else
                {
                    for (k = 0; k < exclPair[i]; k++)
                    {
                        if (exclList[i][k] == j)
                        {
                            exclude = 1;
                            break;
                        }
                    }
                }

                if (!exclude)
                {

                    if (offset >= neigh->sizeList)
                    {
                        fprintf(outFile, "WARNING: List larger than estimated. Size increased from %d", neigh->sizeList);
                        neigh->sizeList = (int) (neigh->sizeList * (1. + TOLLIST)) + 1;
                        fprintf(outFile, " to %d.\n", neigh->sizeList);

                        *neighList = (int*) realloc(*neighList, neigh->sizeList * sizeof (**neighList));
                    }

                    (*neighList)[offset] = j;
                    (*neighPair)[i]++;
                    offset++;

                }
            }
        }
    }

#ifdef TIMER
    update_timer_end(TIMER_VERLET_BUILD, __func__);
#endif

    /*printf("%d %d\n",offset,neigh->sizeList);
    FILE *toto=fopen("list.new","wt");
    for(i=0;i<offset;i++)
    {
      fprintf(toto,"%d\n",(*neighList)[i]);
    }
    fflush(toto);*/
    /*printf("%d %d\n",offset,neigh->sizeList);
    FILE *toto=fopen("list.con","wt");
    for(i=0;i<param->nAtom;i++)
    {
      fprintf(toto,"%d %d\n",(*neighPair)[i],(*neighOrder)[i]);
    }
    fflush(toto);*/

}

void verlet_list_update(PARAM *param, PBC *box, NEIGH *neigh, double x[], double y[], double z[],
        int frozen[], int **neighList, int neighPair[], int **exclList,
        int exclPair[])
{
    int i, j, k, exclude, offset;
    double r2, cutnb2, delta[3];

    cutnb2 = param->cutOff + param->delr;
    cutnb2 = X2(cutnb2);

    offset = 0;

#ifdef TIMER
    update_timer_begin(TIMER_VERLET_UPDATE, __func__);
#endif

#ifdef _OPENMP
#pragma omp parallel default(none) shared(param,neigh,atom,box,cutnb2,neigh->sizeList) private(i,j,k,r,delta,exclude)
    {
#pragma omp for schedule(dynamic,1000) nowait
#endif
        for (i = 0; i < param->nAtom - 1; i++)
        {
            neighPair[i] = 0;
            for (j = i + 1; j < param->nAtom; j++)
            {

                delta[0] = x[j] - x[i];
                delta[1] = y[j] - y[i];
                delta[2] = z[j] - z[i];

                r2 = dist(box, delta);

                if (r2 <= cutnb2)
                {
                    exclude = 0;

                    if ((frozen[i] * frozen[j]))
                    {
                        exclude = 1;
                    }
                    else
                    {
                        for (k = 0; k < exclPair[i]; k++)
                        {
                            if (exclList[i][k] == j)
                            {
                                exclude = 1;
                                break;
                            }
                        }
                    }

                    if (!exclude)
                    {

                        if (offset >= neigh->sizeList)
                        {
                            fprintf(outFile, "WARNING: List larger than estimated. Size increased from %d", neigh->sizeList);
                            neigh->sizeList = (int) (neigh->sizeList * (1. + TOLLIST)) + 1;
                            fprintf(outFile, " to %d.\n", neigh->sizeList);

                            *neighList = (int*) realloc(*neighList, neigh->sizeList * sizeof (**neighList));
                        }

                        (*neighList)[offset] = j;
                        neighPair[i]++;
                        offset++;

                    }
                }
            }
        }
#ifdef _OPENMP
    } // END OF parallel zone
#endif

#ifdef TIMER
    update_timer_end(TIMER_VERLET_UPDATE, __func__);
#endif

}

void link_cell_exclude_list(CTRL *ctrl, PARAM *param, NEIGH *neigh, CONSTRAINT constList[],
        BOND bond[], ANGLE angle[], DIHE dihe[], DIHE impr[],
        int **neighList14, int ***exclList, int **exclPair)
{
    int i, j, k, l, kk, ii, jj, ia, ib, ic, id, exclude;
    int **tempAtom = NULL, **tempVer14 = NULL, **tempConnect = NULL, *tempConnectNum = NULL;
    int nAlloc = 32, nIncr = 16, nConnect = 16;

    *exclPair = (int*) my_malloc(param->nAtom * sizeof (**exclPair));
    for (i = 0; i < param->nAtom; i++)
        (*exclPair)[i] = 0;

    tempAtom = (int**) my_malloc(param->nAtom * sizeof (*(tempAtom)));
    for (i = 0; i < param->nAtom; i++)
        tempAtom[i] = (int*) my_malloc(nAlloc * sizeof (**(tempAtom)));

    /* The following temporary arrays are designed for chechking the
     * connectivity of system and make sure that all existing but not
     * parameterised angles and dihedrals are properly excluded. */

    tempConnectNum = (int*) my_malloc(param->nAtom * sizeof (*tempConnectNum));
    for (i = 0; i < param->nAtom; i++)
        tempConnectNum[i] = 0;

    tempConnect = (int**) my_malloc(param->nAtom * sizeof (*tempConnect));
    for (i = 0; i < param->nAtom; i++)
        tempConnect[i] = (int*) my_malloc(nConnect * sizeof (**tempConnect));


    for (i = 0; i < param->nBond; i++)
    {
        ia = bond[i].a;
        ib = bond[i].b;

        ii = ia;
        jj = ib;

        if ((*exclPair)[ii] >= nAlloc || (*exclPair)[jj] >= nAlloc)
        {
            nAlloc += nIncr;
            for (j = 0; j < param->nAtom; j++)
                tempAtom[j] = (int*) realloc(tempAtom[j], nAlloc * sizeof (**(tempAtom)));
        }

        if (tempConnectNum[ii] >= nConnect || tempConnectNum[jj] >= nConnect)
        {
            nConnect += nIncr;
            for (j = 0; j < param->nAtom; j++)
                tempConnect[j] = (int*) realloc(tempConnect[j], nConnect * sizeof (**tempConnect));
        }

        tempConnect[ii][tempConnectNum[ii]] = jj;
        tempConnect[jj][tempConnectNum[jj]] = ii;
        tempConnectNum[ii]++;
        tempConnectNum[jj]++;

        tempAtom[ii][(*exclPair)[ii]] = jj;
        tempAtom[jj][(*exclPair)[jj]] = ii;
        (*exclPair)[ii]++;
        (*exclPair)[jj]++;
    }

    if (ctrl->keyConstH)
    {
        for (i = 0; i < param->nConst; i++)
        {
            ia = constList[i].a;
            ib = constList[i].b;

            ii = ia;
            jj = ib;

            if ((*exclPair)[ii] >= nAlloc || (*exclPair)[jj] >= nAlloc)
            {
                nAlloc += nIncr;
                for (j = 0; j < param->nAtom; j++)
                    tempAtom[j] = (int*) realloc(tempAtom[j], nAlloc * sizeof (**(tempAtom)));
            }

            if (tempConnectNum[ii] >= nConnect || tempConnectNum[jj] >= nConnect)
            {
                nConnect += nIncr;
                for (j = 0; j < param->nAtom; j++)
                    tempConnect[j] = (int*) realloc(tempConnect[j], nConnect * sizeof (**tempConnect));
            }

            tempConnect[ii][tempConnectNum[ii]] = jj;
            tempConnect[jj][tempConnectNum[jj]] = ii;
            tempConnectNum[ii]++;
            tempConnectNum[jj]++;

            tempAtom[ii][(*exclPair)[ii]] = jj;
            tempAtom[jj][(*exclPair)[jj]] = ii;
            (*exclPair)[ii]++;
            (*exclPair)[jj]++;
        }
    }

    for (i = 0; i < param->nAngle; i++)
    {
        ia = angle[i].a;
        ib = angle[i].b;
        ic = angle[i].c;

        ii = ia;
        jj = ib;

        if ((*exclPair)[ii] >= nAlloc || (*exclPair)[jj] >= nAlloc)
        {
            nAlloc += nIncr;
            for (j = 0; j < param->nAtom; j++)
                tempAtom[j] = (int*) realloc(tempAtom[j], nAlloc * sizeof (**(tempAtom)));
        }

        exclude = 1;
        for (j = 0; j < (*exclPair)[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
        {
            tempAtom[ii][(*exclPair)[ii]] = jj;
            tempAtom[jj][(*exclPair)[jj]] = ii;
            (*exclPair)[ii]++;
            (*exclPair)[jj]++;
        }

        ii = ia;
        jj = ic;

        if ((*exclPair)[ii] >= nAlloc || (*exclPair)[jj] >= nAlloc)
        {
            nAlloc += nIncr;
            for (j = 0; j < param->nAtom; j++)
                tempAtom[j] = (int*) realloc(tempAtom[j], nAlloc * sizeof (**(tempAtom)));
        }

        exclude = 1;
        for (j = 0; j < (*exclPair)[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
        {
            tempAtom[ii][(*exclPair)[ii]] = jj;
            tempAtom[jj][(*exclPair)[jj]] = ii;
            (*exclPair)[ii]++;
            (*exclPair)[jj]++;
        }

        ii = ib;
        jj = ic;

        if ((*exclPair)[ii] >= nAlloc || (*exclPair)[jj] >= nAlloc)
        {
            nAlloc += nIncr;
            for (j = 0; j < param->nAtom; j++)
                tempAtom[j] = (int*) realloc(tempAtom[j], nAlloc * sizeof (**(tempAtom)));
        }

        exclude = 1;
        for (j = 0; j < (*exclPair)[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
        {
            tempAtom[ii][(*exclPair)[ii]] = jj;
            tempAtom[jj][(*exclPair)[jj]] = ii;
            (*exclPair)[ii]++;
            (*exclPair)[jj]++;
        }

    }

    tempVer14 = (int**) my_malloc(5 * param->nDihedral * sizeof (*tempVer14));
    for (i = 0; i < 5 * param->nDihedral; i++)
        tempVer14[i] = (int*) my_malloc(2 * sizeof (**tempVer14));

    neigh->nPair14 = 0;

    for (i = 0; i < param->nDihedral; i++)
    {
        ia = dihe[i].a;
        ib = dihe[i].b;
        ic = dihe[i].c;
        id = dihe[i].d;

        ii = ia;
        jj = ib;

        if ((*exclPair)[ii] >= nAlloc || (*exclPair)[jj] >= nAlloc)
        {
            nAlloc += nIncr;
            for (j = 0; j < param->nAtom; j++)
                tempAtom[j] = (int*) realloc(tempAtom[j], nAlloc * sizeof (**(tempAtom)));
        }

        exclude = 1;
        for (j = 0; j < (*exclPair)[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
        {
            tempAtom[ii][(*exclPair)[ii]] = jj;
            tempAtom[jj][(*exclPair)[jj]] = ii;
            (*exclPair)[ii]++;
            (*exclPair)[jj]++;
        }

        ii = ia;
        jj = ic;

        if ((*exclPair)[ii] >= nAlloc || (*exclPair)[jj] >= nAlloc)
        {
            nAlloc += nIncr;
            for (j = 0; j < param->nAtom; j++)
                tempAtom[j] = (int*) realloc(tempAtom[j], nAlloc * sizeof (**(tempAtom)));
        }

        exclude = 1;
        for (j = 0; j < (*exclPair)[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
        {
            tempAtom[ii][(*exclPair)[ii]] = jj;
            tempAtom[jj][(*exclPair)[jj]] = ii;
            (*exclPair)[ii]++;
            (*exclPair)[jj]++;
        }

        ii = ia;
        jj = id;

        if ((*exclPair)[ii] >= nAlloc || (*exclPair)[jj] >= nAlloc)
        {
            nAlloc += nIncr;
            for (j = 0; j < param->nAtom; j++)
                tempAtom[j] = (int*) realloc(tempAtom[j], nAlloc * sizeof (**(tempAtom)));
        }

        exclude = 1;
        for (j = 0; j < (*exclPair)[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
        {
            tempVer14[neigh->nPair14][0] = ia;
            tempVer14[neigh->nPair14][1] = id;
            neigh->nPair14++;

            tempAtom[ii][(*exclPair)[ii]] = jj;
            tempAtom[jj][(*exclPair)[jj]] = ii;
            (*exclPair)[ii]++;
            (*exclPair)[jj]++;
        }

        ii = ib;
        jj = ic;

        if ((*exclPair)[ii] >= nAlloc || (*exclPair)[jj] >= nAlloc)
        {
            nAlloc += nIncr;
            for (j = 0; j < param->nAtom; j++)
                tempAtom[j] = (int*) realloc(tempAtom[j], nAlloc * sizeof (**(tempAtom)));
        }

        exclude = 1;
        for (j = 0; j < (*exclPair)[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
        {
            tempAtom[ii][(*exclPair)[ii]] = jj;
            tempAtom[jj][(*exclPair)[jj]] = ii;
            (*exclPair)[ii]++;
            (*exclPair)[jj]++;
        }

        ii = ib;
        jj = id;

        if ((*exclPair)[ii] >= nAlloc || (*exclPair)[jj] >= nAlloc)
        {
            nAlloc += nIncr;
            for (j = 0; j < param->nAtom; j++)
                tempAtom[j] = (int*) realloc(tempAtom[j], nAlloc * sizeof (**(tempAtom)));
        }

        exclude = 1;
        for (j = 0; j < (*exclPair)[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
        {
            tempAtom[ii][(*exclPair)[ii]] = jj;
            tempAtom[jj][(*exclPair)[jj]] = ii;
            (*exclPair)[ii]++;
            (*exclPair)[jj]++;
        }

        ii = ic;
        jj = id;

        if ((*exclPair)[ii] >= nAlloc || (*exclPair)[jj] >= nAlloc)
        {
            nAlloc += nIncr;
            for (j = 0; j < param->nAtom; j++)
                tempAtom[j] = (int*) realloc(tempAtom[j], nAlloc * sizeof (**(tempAtom)));
        }

        exclude = 1;
        for (j = 0; j < (*exclPair)[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
        {
            tempAtom[ii][(*exclPair)[ii]] = jj;
            tempAtom[jj][(*exclPair)[jj]] = ii;
            (*exclPair)[ii]++;
            (*exclPair)[jj]++;
        }

    }

    for (i = 0; i < param->nImproper; i++)
    {
        ia = impr[i].a;
        ib = impr[i].b;
        ic = impr[i].c;
        id = impr[i].d;

        ii = ia;
        jj = ib;

        if ((*exclPair)[ii] >= nAlloc || (*exclPair)[jj] >= nAlloc)
        {
            nAlloc += nIncr;
            for (j = 0; j < param->nAtom; j++)
                tempAtom[j] = (int*) realloc(tempAtom[j], nAlloc * sizeof (**(tempAtom)));
        }

        exclude = 1;
        for (j = 0; j < (*exclPair)[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
        {
            tempAtom[ii][(*exclPair)[ii]] = jj;
            tempAtom[jj][(*exclPair)[jj]] = ii;
            (*exclPair)[ii]++;
            (*exclPair)[jj]++;
        }

        ii = ia;
        jj = ic;

        if ((*exclPair)[ii] >= nAlloc || (*exclPair)[jj] >= nAlloc)
        {
            nAlloc += nIncr;
            for (j = 0; j < param->nAtom; j++)
                tempAtom[j] = (int*) realloc(tempAtom[j], nAlloc * sizeof (**(tempAtom)));
        }

        exclude = 1;
        for (j = 0; j < (*exclPair)[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
        {
            tempAtom[ii][(*exclPair)[ii]] = jj;
            tempAtom[jj][(*exclPair)[jj]] = ii;
            (*exclPair)[ii]++;
            (*exclPair)[jj]++;
        }

        ii = ia;
        jj = id;

        if ((*exclPair)[ii] >= nAlloc || (*exclPair)[jj] >= nAlloc)
        {
            nAlloc += nIncr;
            for (j = 0; j < param->nAtom; j++)
                tempAtom[j] = (int*) realloc(tempAtom[j], nAlloc * sizeof (**(tempAtom)));
        }

        exclude = 1;
        for (j = 0; j < (*exclPair)[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
        {
            tempAtom[ii][(*exclPair)[ii]] = jj;
            tempAtom[jj][(*exclPair)[jj]] = ii;
            (*exclPair)[ii]++;
            (*exclPair)[jj]++;
        }

        ii = ib;
        jj = ic;

        if ((*exclPair)[ii] >= nAlloc || (*exclPair)[jj] >= nAlloc)
        {
            nAlloc += nIncr;
            for (j = 0; j < param->nAtom; j++)
                tempAtom[j] = (int*) realloc(tempAtom[j], nAlloc * sizeof (**(tempAtom)));
        }

        exclude = 1;
        for (j = 0; j < (*exclPair)[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
        {
            tempAtom[ii][(*exclPair)[ii]] = jj;
            tempAtom[jj][(*exclPair)[jj]] = ii;
            (*exclPair)[ii]++;
            (*exclPair)[jj]++;
        }

        ii = ib;
        jj = id;

        if ((*exclPair)[ii] >= nAlloc || (*exclPair)[jj] >= nAlloc)
        {
            nAlloc += nIncr;
            for (j = 0; j < param->nAtom; j++)
                tempAtom[j] = (int*) realloc(tempAtom[j], nAlloc * sizeof (**(tempAtom)));
        }

        exclude = 1;
        for (j = 0; j < (*exclPair)[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
        {
            tempAtom[ii][(*exclPair)[ii]] = jj;
            tempAtom[jj][(*exclPair)[jj]] = ii;
            (*exclPair)[ii]++;
            (*exclPair)[jj]++;
        }

        ii = ic;
        jj = id;

        if ((*exclPair)[ii] >= nAlloc || (*exclPair)[jj] >= nAlloc)
        {
            nAlloc += nIncr;
            for (j = 0; j < param->nAtom; j++)
                tempAtom[j] = (int*) realloc(tempAtom[j], nAlloc * sizeof (**(tempAtom)));
        }

        exclude = 1;
        for (j = 0; j < (*exclPair)[ii]; j++)
        {
            if (tempAtom[ii][j] == jj)
            {
                exclude = 0;
                break;
            }
        }

        if (exclude)
        {
            tempAtom[ii][(*exclPair)[ii]] = jj;
            tempAtom[jj][(*exclPair)[jj]] = ii;
            (*exclPair)[ii]++;
            (*exclPair)[jj]++;
        }

    }

    for (ia = 0; ia < param->nAtom; ia++)
    {
        for (j = 0; j < tempConnectNum[ia]; j++)
        {
            ib = tempConnect[ia][j];
            for (k = 0; k < tempConnectNum[ib]; k++)
            {
                ic = tempConnect[ib][k];

                if (ic == ia)
                    continue;

                ii = ia;
                jj = ic;

                if ((*exclPair)[ii] >= nAlloc || (*exclPair)[jj] >= nAlloc)
                {
                    nAlloc += nIncr;
                    for (kk = 0; kk < param->nAtom; kk++)
                        tempAtom[kk] = (int*) realloc(tempAtom[kk], nAlloc * sizeof (**(tempAtom)));
                }

                exclude = 1;
                for (kk = 0; kk < (*exclPair)[ii]; kk++)
                {
                    if (tempAtom[ii][kk] == jj)
                    {
                        exclude = 0;
                        break;
                    }
                }

                if (exclude)
                {
                    tempAtom[ii][(*exclPair)[ii]] = jj;
                    tempAtom[jj][(*exclPair)[jj]] = ii;
                    (*exclPair)[ii]++;
                    (*exclPair)[jj]++;
                }

                for (l = 0; l < tempConnectNum[ic]; l++)
                {
                    id = tempConnect[ic][l];

                    if (id == ia || id == ib)
                        continue;

                    ii = ia;
                    jj = id;

                    if ((*exclPair)[ii] >= nAlloc || (*exclPair)[jj] >= nAlloc)
                    {
                        nAlloc += nIncr;
                        for (kk = 0; kk < param->nAtom; kk++)
                            tempAtom[kk] = (int*) realloc(tempAtom[kk], nAlloc * sizeof (**(tempAtom)));
                    }

                    if (neigh->nPair14 >= 5 * param->nDihedral)
                        my_error(DIHE_NONPARAM_ERROR, __FILE__, __LINE__, 0);

                    exclude = 1;
                    for (kk = 0; kk < (*exclPair)[ii]; kk++)
                    {
                        if (tempAtom[ii][kk] == jj)
                        {
                            exclude = 0;
                            break;
                        }
                    }

                    if (exclude)
                    {
                        tempVer14[neigh->nPair14][0] = ia;
                        tempVer14[neigh->nPair14][1] = id;
                        neigh->nPair14++;

                        tempAtom[ii][(*exclPair)[ii]] = jj;
                        tempAtom[jj][(*exclPair)[jj]] = ii;
                        (*exclPair)[ii]++;
                        (*exclPair)[jj]++;
                    }

                }
            }
        }
    }

    *neighList14 = (int*) my_malloc(2 * neigh->nPair14 * sizeof (**neighList14));

    for (i = 0; i < neigh->nPair14; i++)
    {
        (*neighList14)[2 * i] = tempVer14[i][0];
        (*neighList14)[2 * i + 1] = tempVer14[i][1];
    }

    free_2D(5 * param->nDihedral, tempVer14, NULL);

    free_2D(param->nAtom, tempConnect, NULL);
    free(tempConnectNum);

    *exclList = (int**) my_malloc(param->nAtom * sizeof (**exclList));
    for (i = 0; i < param->nAtom; i++)
        (*exclList)[i] = (int*) my_malloc((*exclPair)[i] * sizeof (***exclList));

    for (i = 0; i < param->nAtom; i++)
    {
        for (j = 0; j < (*exclPair)[i]; j++)
            (*exclList)[i][j] = tempAtom[i][j];
    }

    free_2D(param->nAtom, tempAtom, NULL);

}

void link_cell_verlet_list(PARAM *param, PBC *box, NEIGH *neigh, double x[], double y[], double z[],
        int frozen[], int **neighList, int **neighPair, int **neighOrder,
        int **exclList, int exclPair[])
{

    int transx[508] = {0, 1, 0, 0, -1, 1, 0, -1, 1, 0, -1, 1, -1, 1, 2, 0, 0, -2, 2, -1, 1, 0, -2, 2, 0,
        0, -1, 1, 0, -1, 1, -2, 2, -2, 2, -1, 1, -1, 1, -1, 1, -2, 2, 0, -2, 2, 0, -2, 2, -2, 2,
        -1, 1, -2, 2, -2, 2, -1, 1, -2, 2, -2, 2, 3, 0, 0, -3, 3, -1, 1, 0, -3, 3, 0, 0, -1, 1, 0,
        -1, 1, -3, 3, -3, 3, -1, 1, -1, 1, -1, 1, -3, 3, -2, 2, 0, -3, 3, 0, 0, -2, 2, 0, -2, 2,
        -3, 3, -3, 3, -2, 2, -1, 1, -3, 3, -3, 3, -1, 1, -1, 1, -2, 2, -2, 2, -1, 1, -2, 2, -3, 3,
        -3, 3, -2, 2, -2, 2, -2, 2, -3, 3, 0, -3, 3, 0, -3, 3, -3, 3, -1, 1, -3, 3, -3, 3, -1, 1,
        -3, 3, -3, 3, -2, 2, -3, 3, -3, 3, -2, 2, -3, 3, -3, 3, 4, 0, 0, -4, 4, -1, 1, 0, -4, 4, 0,
        0, -1, 1, 0, -1, 1, -4, 4, -4, 4, -1, 1, -1, 1, -1, 1, -4, 4, -2, 2, 0, -4, 4, 0, 0, -2, 2,
        0, -2, 2, -4, 4, -4, 4, -2, 2, -1, 1, -4, 4, -4, 4, -1, 1, -1, 1, -2, 2, -2, 2, -1, 1, -2,
        2, -4, 4, -4, 4, -2, 2, -2, 2, -2, 2, -4, 4, -3, 3, 0, -4, 4, 0, 0, -3, 3, 0, -3, 3, -4, 4,
        -4, 4, -3, 3, -1, 1, -4, 4, -4, 4, -1, 1, -1, 1, -3, 3, -3, 3, -1, 1, -3, 3, -4, 4, -4, 4,
        -3, 3, -2, 2, -4, 4, -4, 4, -2, 2, -2, 2, -3, 3, -3, 3, -2, 2, -3, 3, -4, 4, -4, 4, -3, 3,
        -3, 3, -3, 3, -4, 4, 0, -4, 4, 0, -4, 4, -4, 4, -1, 1, -4, 4, -4, 4, -1, 1, -4, 4, -4, 4,
        -2, 2, -4, 4, -4, 4, -2, 2, -4, 4, -4, 4, -3, 3, -4, 4, -4, 4, -3, 3, 5, 0, 0, -5, 5, -1,
        1, 0, -5, 5, 0, 0, -1, 1, 0, -1, 1, -5, 5, -5, 5, -1, 1, -1, 1, -1, 1, -5, 5, -2, 2, 0, -5,
        5, 0, 0, -2, 2, 0, -2, 2, -5, 5, -5, 5, -2, 2, -1, 1, -5, 5, -5, 5, -1, 1, -1, 1, -2, 2,
        -2, 2, -1, 1, -2, 2, -5, 5, -5, 5, -2, 2, -2, 2, -2, 2, -5, 5, -3, 3, 0, -5, 5, 0, 0, -3,
        3, 0, -3, 3, -5, 5, -5, 5, -3, 3, -1, 1, -5, 5, -5, 5, -1, 1, -1, 1, -3, 3, -3, 3, -1, 1,
        -3, 3, -5, 5, -5, 5, -3, 3, -2, 2, -5, 5, -5, 5, -2, 2, -2, 2, -3, 3, -3, 3, -2, 2, -3, 3,
        -5, 5, -5, 5, -3, 3, -3, 3, -3, 3};

    int transy[508] = {0, 0, 1, 0, 1, 1, -1, 0, 0, 1, -1, -1, 1, 1, 0, 2, 0, 1, 1, 2, 2, -2, 0, 0, 2,
        -1, 0, 0, 1, -2, -2, -1, -1, 1, 1, 2, 2, -1, -1, 1, 1, 2, 2, -2, 0, 0, 2, -2, -2, 2, 2, -2,
        -2, -1, -1, 1, 1, 2, 2, -2, -2, 2, 2, 0, 3, 0, 1, 1, 3, 3, -3, 0, 0, 3, -1, 0, 0, 1, -3, -3,
        -1, -1, 1, 1, 3, 3, -1, -1, 1, 1, 2, 2, 3, 3, -3, 0, 0, 3, -2, 0, 0, 2, -3, -3, -2, -2, 2,
        2, 3, 3, -3, -3, -1, -1, 1, 1, 3, 3, -2, -2, -1, -1, 1, 1, 2, 2, -3, -3, -2, -2, 2, 2, 3,
        3, -2, -2, 2, 2, 3, 3, -3, 0, 0, 3, -3, -3, 3, 3, -3, -3, -1, -1, 1, 1, 3, 3, -3, -3, 3, 3,
        -3, -3, -2, -2, 2, 2, 3, 3, -3, -3, 3, 3, 0, 4, 0, 1, 1, 4, 4, -4, 0, 0, 4, -1, 0, 0, 1, -4,
        -4, -1, -1, 1, 1, 4, 4, -1, -1, 1, 1, 2, 2, 4, 4, -4, 0, 0, 4, -2, 0, 0, 2, -4, -4, -2, -2,
        2, 2, 4, 4, -4, -4, -1, -1, 1, 1, 4, 4, -2, -2, -1, -1, 1, 1, 2, 2, -4, -4, -2, -2, 2, 2,
        4, 4, -2, -2, 2, 2, 3, 3, 4, 4, -4, 0, 0, 4, -3, 0, 0, 3, -4, -4, -3, -3, 3, 3, 4, 4, -4,
        -4, -1, -1, 1, 1, 4, 4, -3, -3, -1, -1, 1, 1, 3, 3, -4, -4, -3, -3, 3, 3, 4, 4, -4, -4,
        -2, -2, 2, 2, 4, 4, -3, -3, -2, -2, 2, 2, 3, 3, -4, -4, -3, -3, 3, 3, 4, 4, -3, -3, 3, 3,
        4, 4, -4, 0, 0, 4, -4, -4, 4, 4, -4, -4, -1, -1, 1, 1, 4, 4, -4, -4, 4, 4, -4, -4, -2, -2,
        2, 2, 4, 4, -4, -4, 4, 4, -4, -4, -3, -3, 3, 3, 4, 4, 0, 5, 0, 1, 1, 5, 5, -5, 0, 0, 5, -1,
        0, 0, 1, -5, -5, -1, -1, 1, 1, 5, 5, -1, -1, 1, 1, 2, 2, 5, 5, -5, 0, 0, 5, -2, 0, 0, 2, -5,
        -5, -2, -2, 2, 2, 5, 5, -5, -5, -1, -1, 1, 1, 5, 5, -2, -2, -1, -1, 1, 1, 2, 2, -5, -5,
        -2, -2, 2, 2, 5, 5, -2, -2, 2, 2, 3, 3, 5, 5, -5, 0, 0, 5, -3, 0, 0, 3, -5, -5, -3, -3, 3,
        3, 5, 5, -5, -5, -1, -1, 1, 1, 5, 5, -3, -3, -1, -1, 1, 1, 3, 3, -5, -5, -3, -3, 3, 3, 5,
        5, -5, -5, -2, -2, 2, 2, 5, 5, -3, -3, -2, -2, 2, 2, 3, 3, -5, -5, -3, -3, 3, 3, 5, 5, -3,
        -3, 3, 3};

    int transz[508] = {0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 2, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2,
        2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 0, 0, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 0, 0, 3, 0, 0, 0, 0, 1, 1, 1, 1, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 0, 0,
        0, 0, 2, 2, 2, 2, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3,
        3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 0, 0, 3, 3, 3, 3, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3,
        2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 0, 4, 0, 0, 0, 0, 1, 1, 1, 1, 4, 4, 4, 4, 1,
        1, 1, 1, 1, 1, 1, 1, 4, 4, 4, 4, 0, 0, 0, 0, 2, 2, 2, 2, 4, 4, 4, 4, 1, 1, 1, 1, 1, 1, 1, 1, 2,
        2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 4, 4, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 0, 0, 0, 0, 3,
        3, 3, 3, 4, 4, 4, 4, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 2,
        2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 4,
        4, 4, 4, 0, 0, 4, 4, 4, 4, 1, 1, 1, 1, 4, 4, 4, 4, 4, 4, 4, 4, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 4,
        4, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 0, 0, 5, 0, 0, 0, 0, 1, 1, 1, 1, 5, 5, 5, 5, 1, 1, 1, 1,
        1, 1, 1, 1, 5, 5, 5, 5, 0, 0, 0, 0, 2, 2, 2, 2, 5, 5, 5, 5, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2,
        2, 2, 2, 2, 5, 5, 5, 5, 5, 5, 5, 5, 2, 2, 2, 2, 2, 2, 2, 2, 5, 5, 5, 5, 0, 0, 0, 0, 3, 3, 3, 3,
        5, 5, 5, 5, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3, 5, 5, 5, 5, 5, 5, 5, 5, 2, 2, 2, 2,
        2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 5, 5, 5, 5, 5, 5, 5, 5, 3, 3, 3, 3, 3, 3, 3, 3, 5, 5, 5, 5};

    int cellCheck, enoughLinkCells, nlcx, nlcy, nlcz;
    int i, j, k, l, kk, ll, icell, exclude, offset, iorder;
    int ix, iy, iz, jx, jy, jz;
    double r2, cutnb, cutnb2, dnlcx, dnlcy, dnlcz;
    double cx, cy, cz, xt, yt, zt, xd, yd, zd, *xu, *yu, *zu;

#ifdef TIMER
    update_timer_begin(TIMER_LNKCEL_BUILD, __func__);
#endif

    xu = (double*) my_malloc(param->nAtom * sizeof (*xu));
    yu = (double*) my_malloc(param->nAtom * sizeof (*yu));
    zu = (double*) my_malloc(param->nAtom * sizeof (*zu));

    cutnb = param->cutOff + param->delr;
    cutnb2 = X2(cutnb);

    neigh->linkRatio = MIN(MAX(neigh->linkRatio, 1), 5);

    if (neigh->linkRatio == 1)
        cellCheck = 14;
    else if (neigh->linkRatio == 2)
        cellCheck = 63;
    else if (neigh->linkRatio == 3)
        cellCheck = 156;
    else if (neigh->linkRatio == 4)
        cellCheck = 307;
    else if (neigh->linkRatio == 5)
        cellCheck = 508;
    else
        my_error(LNKCELL_CUTOFF_ERROR, __FILE__, __LINE__, 0);

    nlcx = (int) (box->pa * (double) neigh->linkRatio / cutnb);
    nlcy = (int) (box->pb * (double) neigh->linkRatio / cutnb);
    nlcz = (int) (box->pc * (double) neigh->linkRatio / cutnb);

    dnlcx = (double) nlcx;
    dnlcy = (double) nlcy;
    dnlcz = (double) nlcz;

    enoughLinkCells = 1;
    if (nlcx < 2 * neigh->linkRatio + 1)
        enoughLinkCells = 0;
    if (nlcy < 2 * neigh->linkRatio + 1)
        enoughLinkCells = 0;
    if (nlcz < 2 * neigh->linkRatio + 1)
        enoughLinkCells = 0;
    if (!enoughLinkCells)
        my_error(LNKCELL_NCELLS_ERROR, __FILE__, __LINE__, 0);

    neigh->nCells = nlcx * nlcy*nlcz;

    int* link = (int*) my_malloc(param->nAtom * sizeof (*link));

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(param,neigh->sizeList,link) private(i)
#endif
    for (i = 0; i < param->nAtom; i++)
    {
        link[i] = -1;
    }

    int* head = (int*) my_malloc(neigh->nCells * sizeof (*head));

    for (i = 0; i < neigh->nCells; i++)
    {
        head[i] = -1;
    }

    //   image_update(atom,simulCond,box); //if new job.
    image_update(param, box, x, y, z); //if new job.

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(param,atom,box,xu,yu,zu) private(i)
#endif
    for (i = 0; i < param->nAtom; i++)
    {
        xu[i] = (x[i] * box->u1 + y[i] * box->u2 + z[i] * box->u3) + 0.5;
        yu[i] = (x[i] * box->v1 + y[i] * box->v2 + z[i] * box->v3) + 0.5;
        zu[i] = (x[i] * box->w1 + y[i] * box->w2 + z[i] * box->w3) + 0.5;
    }

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(param,xu,yu,zu,nlcx,nlcy,nlcz,dnlcx,dnlcy,dnlcz,head,link) private(i,ix,iy,iz,icell,j)
#endif
    for (i = 0; i < param->nAtom; i++)
    {

        ix = (int) dnlcx * xu[i];
        ix = MIN(ix, nlcx - 1);

        iy = (int) dnlcy * yu[i];
        iy = MIN(iy, nlcy - 1);

        iz = (int) dnlcz * zu[i];
        iz = MIN(iz, nlcz - 1);

        icell = ix + nlcx * (iy + nlcy * iz);

#ifdef _OPENMP
#pragma omp critical
        {
#endif
            j = head[icell];
            head[icell] = i;
            link[i] = j;
#ifdef _OPENMP
        }
#endif

    }//end for

    *neighOrder = (int*) my_malloc(param->nAtom * sizeof (**neighOrder));

    iorder = 0;
    for (k = 0; k < neigh->nCells; k++)
    {
        i = head[k];
        while (i != -1)
        {
            (*neighOrder)[iorder++] = i;
            i = link[i];
        }
    }

    ix = 0;
    iy = 0;
    iz = 0;

    neigh->sizeList = 0;

    for (k = 0; k < neigh->nCells; k++)
    {
        i = head[k];

        while (i != -1)
        {

            for (l = 0; l < cellCheck; l++)
            {

                cx = 0.;
                cy = 0.;
                cz = 0.;

                jx = ix + transx[l];
                jy = iy + transy[l];
                jz = iz + transz[l];

                if (jx > nlcx - 1)
                {
                    jx = jx - nlcx;
                    cx = 1.;
                }
                else if (jx < 0)
                {
                    jx = jx + nlcx;
                    cx = -1.;
                }

                if (jy > nlcy - 1)
                {
                    jy = jy - nlcy;
                    cy = 1.;
                }
                else if (jy < 0)
                {
                    jy = jy + nlcy;
                    cy = -1.;
                }

                if (jz > nlcz - 1)
                {
                    jz = jz - nlcz;
                    cz = 1.;
                }
                else if (jz < 0)
                {
                    jz = jz + nlcz;
                    cz = -1.;
                }

                ll = jx + nlcx * (jy + nlcy * jz);

                j = head[ll];

                if (k == ll)j = link[i];

                while (j != -1)
                {
                    xt = xu[j] - xu[i] + cx;
                    yt = yu[j] - yu[i] + cy;
                    zt = zu[j] - zu[i] + cz;

                    xd = xt * box->a1 + yt * box->b1 + zt * box->c1;
                    yd = xt * box->a2 + yt * box->b2 + zt * box->c2;
                    zd = xt * box->a3 + yt * box->b3 + zt * box->c3;

                    r2 = X2(xd) + X2(yd) + X2(zd);

                    if (r2 <= cutnb2)
                    {
                        exclude = 0;

                        if ((frozen[i] * frozen[j]))
                        {
                            exclude = 1;
                        }
                        else
                        {
                            for (kk = 0; kk < exclPair[i]; kk++)
                            {
                                if (exclList[i][kk] == j)
                                {
                                    exclude = 1;
                                    break;
                                }
                            }
                        }

                        if (!exclude)
                        {

                            neigh->sizeList++;

                        } //exclude
                    }//cutnb

                    j = link[j];

                }//while j
            }//for checkcell

            i = link[i];

        }// while i

        ix = ix + 1;
        if (ix > nlcx - 1)
        {
            ix = 0;
            iy = iy + 1;

            if (iy > nlcy - 1)
            {
                iy = 0;
                iz = iz + 1;
            }

        }

    }//for ncells

    neigh->sizeList = (int) ((double) neigh->sizeList * (1. + 2. * TOLLIST)) + 1;

    *neighList = (int*) my_malloc(neigh->sizeList * sizeof (**neighList));

    *neighPair = (int*) my_malloc(param->nAtom * sizeof (**neighPair));

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(param,neigh,neigh->sizeList) private(i)
#endif
    for (i = 0; i < param->nAtom; i++)
    {
        (*neighPair)[i] = 0;
    }

    ix = 0;
    iy = 0;
    iz = 0;

    offset = 0;

    for (k = 0; k < neigh->nCells; k++)
    {
        i = head[k];

        while (i != -1)
        {

            for (l = 0; l < cellCheck; l++)
            {

                cx = 0.;
                cy = 0.;
                cz = 0.;

                jx = ix + transx[l];
                jy = iy + transy[l];
                jz = iz + transz[l];

                if (jx > nlcx - 1)
                {
                    jx = jx - nlcx;
                    cx = 1.;
                }
                else if (jx < 0)
                {
                    jx = jx + nlcx;
                    cx = -1.;
                }

                if (jy > nlcy - 1)
                {
                    jy = jy - nlcy;
                    cy = 1.;
                }
                else if (jy < 0)
                {
                    jy = jy + nlcy;
                    cy = -1.;
                }

                if (jz > nlcz - 1)
                {
                    jz = jz - nlcz;
                    cz = 1.;
                }
                else if (jz < 0)
                {
                    jz = jz + nlcz;
                    cz = -1.;
                }

                ll = jx + nlcx * (jy + nlcy * jz);

                j = head[ll];

                if (k == ll)j = link[i];

                while (j != -1)
                {
                    xt = xu[j] - xu[i] + cx;
                    yt = yu[j] - yu[i] + cy;
                    zt = zu[j] - zu[i] + cz;

                    xd = xt * box->a1 + yt * box->b1 + zt * box->c1;
                    yd = xt * box->a2 + yt * box->b2 + zt * box->c2;
                    zd = xt * box->a3 + yt * box->b3 + zt * box->c3;

                    r2 = X2(xd) + X2(yd) + X2(zd);

                    if (r2 <= cutnb2)
                    {
                        exclude = 0;

                        if ((frozen[i] * frozen[j]))
                        {
                            exclude = 1;
                        }
                        else
                        {
                            for (kk = 0; kk < exclPair[i]; kk++)
                            {
                                if (exclList[i][kk] == j)
                                {
                                    exclude = 1;
                                    break;
                                }
                            }
                        }

                        if (!exclude)
                        {

                            if (offset >= neigh->sizeList)
                            {
                                fprintf(outFile, "WARNING: List larger than estimated. Size increased from %d", neigh->sizeList);
                                neigh->sizeList = (int) ((double) neigh->sizeList * (1. + TOLLIST)) + 1;
                                fprintf(outFile, " to %d.\n", neigh->sizeList);

                                *neighList = (int*) realloc(*neighList, neigh->sizeList * sizeof (**neighList));
                            }

                            (*neighList)[offset] = j;
                            (*neighPair)[i]++;
                            offset++;

                        } //exclude
                    }//cutnb

                    j = link[j];

                }//while j
            }//for checkcell

            i = link[i];

        }// while i

        ix = ix + 1;
        if (ix > nlcx - 1)
        {
            ix = 0;
            iy = iy + 1;

            if (iy > nlcy - 1)
            {
                iy = 0;
                iz = iz + 1;
            }

        }

    }//for ncells

    free(xu);
    free(yu);
    free(zu);

    free(link);
    free(head);

#ifdef TIMER
    update_timer_end(TIMER_LNKCEL_BUILD, __func__);
#endif

}

void link_cell_verlet_list_update(PARAM *param, PBC *box, NEIGH *neigh, double x[], double y[], double z[],
        int frozen[], int **neighList, int neighPair[], int neighOrder[],
        int **exclList, int exclPair[])
{

    int transx[508] = {0, 1, 0, 0, -1, 1, 0, -1, 1, 0, -1, 1, -1, 1, 2, 0, 0, -2, 2, -1, 1, 0, -2, 2, 0,
        0, -1, 1, 0, -1, 1, -2, 2, -2, 2, -1, 1, -1, 1, -1, 1, -2, 2, 0, -2, 2, 0, -2, 2, -2, 2,
        -1, 1, -2, 2, -2, 2, -1, 1, -2, 2, -2, 2, 3, 0, 0, -3, 3, -1, 1, 0, -3, 3, 0, 0, -1, 1, 0,
        -1, 1, -3, 3, -3, 3, -1, 1, -1, 1, -1, 1, -3, 3, -2, 2, 0, -3, 3, 0, 0, -2, 2, 0, -2, 2,
        -3, 3, -3, 3, -2, 2, -1, 1, -3, 3, -3, 3, -1, 1, -1, 1, -2, 2, -2, 2, -1, 1, -2, 2, -3, 3,
        -3, 3, -2, 2, -2, 2, -2, 2, -3, 3, 0, -3, 3, 0, -3, 3, -3, 3, -1, 1, -3, 3, -3, 3, -1, 1,
        -3, 3, -3, 3, -2, 2, -3, 3, -3, 3, -2, 2, -3, 3, -3, 3, 4, 0, 0, -4, 4, -1, 1, 0, -4, 4, 0,
        0, -1, 1, 0, -1, 1, -4, 4, -4, 4, -1, 1, -1, 1, -1, 1, -4, 4, -2, 2, 0, -4, 4, 0, 0, -2, 2,
        0, -2, 2, -4, 4, -4, 4, -2, 2, -1, 1, -4, 4, -4, 4, -1, 1, -1, 1, -2, 2, -2, 2, -1, 1, -2,
        2, -4, 4, -4, 4, -2, 2, -2, 2, -2, 2, -4, 4, -3, 3, 0, -4, 4, 0, 0, -3, 3, 0, -3, 3, -4, 4,
        -4, 4, -3, 3, -1, 1, -4, 4, -4, 4, -1, 1, -1, 1, -3, 3, -3, 3, -1, 1, -3, 3, -4, 4, -4, 4,
        -3, 3, -2, 2, -4, 4, -4, 4, -2, 2, -2, 2, -3, 3, -3, 3, -2, 2, -3, 3, -4, 4, -4, 4, -3, 3,
        -3, 3, -3, 3, -4, 4, 0, -4, 4, 0, -4, 4, -4, 4, -1, 1, -4, 4, -4, 4, -1, 1, -4, 4, -4, 4,
        -2, 2, -4, 4, -4, 4, -2, 2, -4, 4, -4, 4, -3, 3, -4, 4, -4, 4, -3, 3, 5, 0, 0, -5, 5, -1,
        1, 0, -5, 5, 0, 0, -1, 1, 0, -1, 1, -5, 5, -5, 5, -1, 1, -1, 1, -1, 1, -5, 5, -2, 2, 0, -5,
        5, 0, 0, -2, 2, 0, -2, 2, -5, 5, -5, 5, -2, 2, -1, 1, -5, 5, -5, 5, -1, 1, -1, 1, -2, 2,
        -2, 2, -1, 1, -2, 2, -5, 5, -5, 5, -2, 2, -2, 2, -2, 2, -5, 5, -3, 3, 0, -5, 5, 0, 0, -3,
        3, 0, -3, 3, -5, 5, -5, 5, -3, 3, -1, 1, -5, 5, -5, 5, -1, 1, -1, 1, -3, 3, -3, 3, -1, 1,
        -3, 3, -5, 5, -5, 5, -3, 3, -2, 2, -5, 5, -5, 5, -2, 2, -2, 2, -3, 3, -3, 3, -2, 2, -3, 3,
        -5, 5, -5, 5, -3, 3, -3, 3, -3, 3};

    int transy[508] = {0, 0, 1, 0, 1, 1, -1, 0, 0, 1, -1, -1, 1, 1, 0, 2, 0, 1, 1, 2, 2, -2, 0, 0, 2,
        -1, 0, 0, 1, -2, -2, -1, -1, 1, 1, 2, 2, -1, -1, 1, 1, 2, 2, -2, 0, 0, 2, -2, -2, 2, 2, -2,
        -2, -1, -1, 1, 1, 2, 2, -2, -2, 2, 2, 0, 3, 0, 1, 1, 3, 3, -3, 0, 0, 3, -1, 0, 0, 1, -3, -3,
        -1, -1, 1, 1, 3, 3, -1, -1, 1, 1, 2, 2, 3, 3, -3, 0, 0, 3, -2, 0, 0, 2, -3, -3, -2, -2, 2,
        2, 3, 3, -3, -3, -1, -1, 1, 1, 3, 3, -2, -2, -1, -1, 1, 1, 2, 2, -3, -3, -2, -2, 2, 2, 3,
        3, -2, -2, 2, 2, 3, 3, -3, 0, 0, 3, -3, -3, 3, 3, -3, -3, -1, -1, 1, 1, 3, 3, -3, -3, 3, 3,
        -3, -3, -2, -2, 2, 2, 3, 3, -3, -3, 3, 3, 0, 4, 0, 1, 1, 4, 4, -4, 0, 0, 4, -1, 0, 0, 1, -4,
        -4, -1, -1, 1, 1, 4, 4, -1, -1, 1, 1, 2, 2, 4, 4, -4, 0, 0, 4, -2, 0, 0, 2, -4, -4, -2, -2,
        2, 2, 4, 4, -4, -4, -1, -1, 1, 1, 4, 4, -2, -2, -1, -1, 1, 1, 2, 2, -4, -4, -2, -2, 2, 2,
        4, 4, -2, -2, 2, 2, 3, 3, 4, 4, -4, 0, 0, 4, -3, 0, 0, 3, -4, -4, -3, -3, 3, 3, 4, 4, -4,
        -4, -1, -1, 1, 1, 4, 4, -3, -3, -1, -1, 1, 1, 3, 3, -4, -4, -3, -3, 3, 3, 4, 4, -4, -4,
        -2, -2, 2, 2, 4, 4, -3, -3, -2, -2, 2, 2, 3, 3, -4, -4, -3, -3, 3, 3, 4, 4, -3, -3, 3, 3,
        4, 4, -4, 0, 0, 4, -4, -4, 4, 4, -4, -4, -1, -1, 1, 1, 4, 4, -4, -4, 4, 4, -4, -4, -2, -2,
        2, 2, 4, 4, -4, -4, 4, 4, -4, -4, -3, -3, 3, 3, 4, 4, 0, 5, 0, 1, 1, 5, 5, -5, 0, 0, 5, -1,
        0, 0, 1, -5, -5, -1, -1, 1, 1, 5, 5, -1, -1, 1, 1, 2, 2, 5, 5, -5, 0, 0, 5, -2, 0, 0, 2, -5,
        -5, -2, -2, 2, 2, 5, 5, -5, -5, -1, -1, 1, 1, 5, 5, -2, -2, -1, -1, 1, 1, 2, 2, -5, -5,
        -2, -2, 2, 2, 5, 5, -2, -2, 2, 2, 3, 3, 5, 5, -5, 0, 0, 5, -3, 0, 0, 3, -5, -5, -3, -3, 3,
        3, 5, 5, -5, -5, -1, -1, 1, 1, 5, 5, -3, -3, -1, -1, 1, 1, 3, 3, -5, -5, -3, -3, 3, 3, 5,
        5, -5, -5, -2, -2, 2, 2, 5, 5, -3, -3, -2, -2, 2, 2, 3, 3, -5, -5, -3, -3, 3, 3, 5, 5, -3,
        -3, 3, 3};

    int transz[508] = {0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 2, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2,
        2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 0, 0, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 0, 0, 3, 0, 0, 0, 0, 1, 1, 1, 1, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 0, 0,
        0, 0, 2, 2, 2, 2, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3,
        3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 0, 0, 3, 3, 3, 3, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3,
        2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 0, 4, 0, 0, 0, 0, 1, 1, 1, 1, 4, 4, 4, 4, 1,
        1, 1, 1, 1, 1, 1, 1, 4, 4, 4, 4, 0, 0, 0, 0, 2, 2, 2, 2, 4, 4, 4, 4, 1, 1, 1, 1, 1, 1, 1, 1, 2,
        2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 4, 4, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 0, 0, 0, 0, 3,
        3, 3, 3, 4, 4, 4, 4, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 2,
        2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 4,
        4, 4, 4, 0, 0, 4, 4, 4, 4, 1, 1, 1, 1, 4, 4, 4, 4, 4, 4, 4, 4, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 4,
        4, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 0, 0, 5, 0, 0, 0, 0, 1, 1, 1, 1, 5, 5, 5, 5, 1, 1, 1, 1,
        1, 1, 1, 1, 5, 5, 5, 5, 0, 0, 0, 0, 2, 2, 2, 2, 5, 5, 5, 5, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2,
        2, 2, 2, 2, 5, 5, 5, 5, 5, 5, 5, 5, 2, 2, 2, 2, 2, 2, 2, 2, 5, 5, 5, 5, 0, 0, 0, 0, 3, 3, 3, 3,
        5, 5, 5, 5, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3, 5, 5, 5, 5, 5, 5, 5, 5, 2, 2, 2, 2,
        2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 5, 5, 5, 5, 5, 5, 5, 5, 3, 3, 3, 3, 3, 3, 3, 3, 5, 5, 5, 5};

    int cellCheck, enoughLinkCells, nlcx, nlcy, nlcz;
    int i, j, k, l, ll, kk, icell, exclude, offset, iorder;
    int ix, iy, iz, jx, jy, jz;
    double r2, cutnb, cutnb2, dnlcx, dnlcy, dnlcz;
    double cx, cy, cz, xt, yt, zt, xd, yd, zd, *xu, *yu, *zu;

#ifdef TIMER
    update_timer_begin(TIMER_LNKCEL_UPDATE, __func__);
#endif

    xu = (double*) my_malloc(param->nAtom * sizeof (*xu));
    yu = (double*) my_malloc(param->nAtom * sizeof (*yu));
    zu = (double*) my_malloc(param->nAtom * sizeof (*zu));

    cutnb = param->cutOff + param->delr;
    cutnb2 = X2(cutnb);

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(param,neigh,neigh->sizeList) private(i)
#endif
    for (i = 0; i < param->nAtom; i++)
    {
        neighPair[i] = 0;
    }

    neigh->linkRatio = MIN(MAX(neigh->linkRatio, 1), 5);

    if (neigh->linkRatio == 1)
        cellCheck = 14;
    else if (neigh->linkRatio == 2)
        cellCheck = 63;
    else if (neigh->linkRatio == 3)
        cellCheck = 156;
    else if (neigh->linkRatio == 4)
        cellCheck = 307;
    else if (neigh->linkRatio == 5)
        cellCheck = 508;
    else
        my_error(LNKCELL_CUTOFF_ERROR, __FILE__, __LINE__, 0);

    nlcx = (int) (box->pa * (double) neigh->linkRatio / cutnb);
    nlcy = (int) (box->pb * (double) neigh->linkRatio / cutnb);
    nlcz = (int) (box->pc * (double) neigh->linkRatio / cutnb);

    dnlcx = (double) nlcx;
    dnlcy = (double) nlcy;
    dnlcz = (double) nlcz;

    enoughLinkCells = 1;
    if (nlcx < 2 * neigh->linkRatio + 1)
        enoughLinkCells = 0;
    if (nlcy < 2 * neigh->linkRatio + 1)
        enoughLinkCells = 0;
    if (nlcz < 2 * neigh->linkRatio + 1)
        enoughLinkCells = 0;
    if (!enoughLinkCells)
        my_error(LNKCELL_NCELLS_ERROR, __FILE__, __LINE__, 0);

    neigh->nCells = nlcx * nlcy*nlcz;

    int* link = (int*) my_malloc(param->nAtom * sizeof (*link));

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(param,neigh->sizeList,link) private(i)
#endif
    for (i = 0; i < param->nAtom; i++)
    {
        link[i] = -1;
    }

    int* head = (int*) my_malloc(neigh->nCells * sizeof (*head));

    for (i = 0; i < neigh->nCells; i++)
    {
        head[i] = -1;
    }


#ifdef _OPENMP
#pragma omp parallel for default(none) shared(param,atom,box,xu,yu,zu) private(i)
#endif
    for (i = 0; i < param->nAtom; i++)
    {
        xu[i] = (x[i] * box->u1 + y[i] * box->u2 + z[i] * box->u3) + 0.5;
        yu[i] = (x[i] * box->v1 + y[i] * box->v2 + z[i] * box->v3) + 0.5;
        zu[i] = (x[i] * box->w1 + y[i] * box->w2 + z[i] * box->w3) + 0.5;
    }

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(param,xu,yu,zu,nlcx,nlcy,nlcz,dnlcx,dnlcy,dnlcz,head,link) private(i,ix,iy,iz,icell,j)
#endif
    for (i = 0; i < param->nAtom; i++)
    {

        ix = (int) dnlcx * xu[i];
        ix = MIN(ix, nlcx - 1);

        iy = (int) dnlcy * yu[i];
        iy = MIN(iy, nlcy - 1);

        iz = (int) dnlcz * zu[i];
        iz = MIN(iz, nlcz - 1);

        icell = ix + nlcx * (iy + nlcy * iz);

#ifdef _OPENMP
#pragma omp critical
        {
#endif
            j = head[icell];
            head[icell] = i;
            link[i] = j;
#ifdef _OPENMP
        }
#endif

    } //end for

    iorder = 0;
    for (k = 0; k < neigh->nCells; k++)
    {
        i = head[k];
        while (i != -1)
        {
            neighOrder[iorder++] = i;
            i = link[i];
        }
    }

    ix = 0;
    iy = 0;
    iz = 0;

    offset = 0;

    for (k = 0; k < neigh->nCells; k++)
    {
        i = head[k];

        while (i != -1)
        {

            for (l = 0; l < cellCheck; l++)
            {

                cx = 0.;
                cy = 0.;
                cz = 0.;

                jx = ix + transx[l];
                jy = iy + transy[l];
                jz = iz + transz[l];

                if (jx > nlcx - 1)
                {
                    jx = jx - nlcx;
                    cx = 1.;
                }
                else if (jx < 0)
                {
                    jx = jx + nlcx;
                    cx = -1.;
                }

                if (jy > nlcy - 1)
                {
                    jy = jy - nlcy;
                    cy = 1.;
                }
                else if (jy < 0)
                {
                    jy = jy + nlcy;
                    cy = -1.;
                }

                if (jz > nlcz - 1)
                {
                    jz = jz - nlcz;
                    cz = 1.;
                }
                else if (jz < 0)
                {
                    jz = jz + nlcz;
                    cz = -1.;
                }

                ll = jx + nlcx * (jy + nlcy * jz);

                j = head[ll];

                if (k == ll)j = link[i];

                while (j != -1)
                {
                    xt = xu[j] - xu[i] + cx;
                    yt = yu[j] - yu[i] + cy;
                    zt = zu[j] - zu[i] + cz;

                    xd = xt * box->a1 + yt * box->b1 + zt * box->c1;
                    yd = xt * box->a2 + yt * box->b2 + zt * box->c2;
                    zd = xt * box->a3 + yt * box->b3 + zt * box->c3;

                    r2 = X2(xd) + X2(yd) + X2(zd);

                    if (r2 <= cutnb2)
                    {
                        exclude = 0;

                        if ((frozen[i] * frozen[j]))
                        {
                            exclude = 1;
                        }
                        else
                        {
                            for (kk = 0; kk < exclPair[i]; kk++)
                            {
                                if (exclList[i][kk] == j)
                                {
                                    exclude = 1;
                                    break;
                                }
                            }
                        }

                        if (!exclude)
                        {

                            if (offset >= neigh->sizeList)
                            {
                                fprintf(outFile, "WARNING: List larger than estimated. Size increased from %d", neigh->sizeList);
                                neigh->sizeList = (int) ((double) neigh->sizeList * (1. + TOLLIST)) + 1;
                                fprintf(outFile, " to %d.\n", neigh->sizeList);

                                *neighList = (int*) realloc(*neighList, neigh->sizeList * sizeof (**neighList));
                            }

                            (*neighList)[offset] = j;
                            neighPair[i]++;
                            offset++;

                        } //exclude
                    }//cutnb

                    j = link[j];

                }//while j
            }//for checkcell

            i = link[i];

        }// while i

        ix = ix + 1;
        if (ix > nlcx - 1)
        {
            ix = 0;
            iy = iy + 1;

            if (iy > nlcy - 1)
            {
                iy = 0;
                iz = iz + 1;
            }

        }

    }//for ncells

    free(xu);
    free(yu);
    free(zu);

    free(link);
    free(head);

#ifdef TIMER
    update_timer_end(TIMER_LNKCEL_UPDATE, __func__);
#endif

}

void fast_verlet_list(PARAM *param, PBC *box, NEIGH *neigh, double x[], double y[], double z[],
        int frozen[], int **neighList, int **neighPair, int **neighOrder,
        int **exclList, int exclPair[])
{
    /**
     * This routine is derived from the algorithm described
     * by Tim. N. Heinz and Philippe H. Hunenberger
     * J. Comput. Chem., Vol. 25, No. 12, 1474--1486 (2004)
     * DOI: 10.1002/jcc.20071
     */

    int i, j, k, kk, m, dm, s, m1, m2, stripes, atpercell;
    int nlcx, nlcy, nlcz;
    int dmx, dmy, dmz, dnx, dny, dnz;
    int t1ny, t2ny, t1nz, t2nz;
    int exclude, offset, iorder;

    double r2, rlcx, rlcy, rlcz, rlcx2, rlcy2, rlcz2, cutnb, cutnb2;
    double dnlcx, dnlcy, dnlcz;

    int ix, iy, iz;

    double delta[3];

    int *ptrmask = NULL, *ptrcell = NULL, *cell = NULL, *tempcell = NULL;
    double *xu = NULL, *yu = NULL, *zu = NULL;

    xu = (double*) my_malloc(param->nAtom * sizeof (*xu));
    yu = (double*) my_malloc(param->nAtom * sizeof (*yu));
    zu = (double*) my_malloc(param->nAtom * sizeof (*zu));
    ;

    cutnb = param->cutOff + param->delr;
    cutnb2 = X2(cutnb);

    nlcx = (int) (box->pa * (double) neigh->linkRatio / cutnb);
    nlcy = (int) (box->pb * (double) neigh->linkRatio / cutnb);
    nlcz = (int) (box->pc * (double) neigh->linkRatio / cutnb);

    neigh->nCells = nlcx * nlcy*nlcz;

    dnlcx = (double) nlcx;
    dnlcy = (double) nlcy;
    dnlcz = (double) nlcz;

    rlcx = box->pa*dnlcx;
    rlcy = box->pb*dnlcy;
    rlcz = box->pc*dnlcz;

    rlcx2 = X2(rlcx);
    rlcy2 = X2(rlcy);
    rlcz2 = X2(rlcz);

    stripes = (int) (1.5 * PI * cutnb2 * dnlcy * dnlcz / (box->pb * box->pc));
    atpercell = (int) (1.5 * param->nAtom / neigh->nCells);

    ptrmask = (int*) my_malloc(2 * stripes * sizeof (*ptrmask));
    for (i = 0; i < 2 * stripes; i++)
        ptrmask[i] = -1;

    cell = (int*) my_malloc(param->nAtom * sizeof (*cell));
    tempcell = (int*) my_malloc(neigh->nCells * atpercell * sizeof (*tempcell));

    ptrcell = (int*) my_malloc((neigh->nCells + 1) * sizeof (*ptrcell));
    for (i = 0; i < neigh->nCells + 1; i++)
        ptrcell[i] = i * atpercell;

    ptrcell[neigh->nCells] = param->nAtom;

    k = 0;
    for (dm = 1; dm < neigh->nCells - 1; dm++)
    {

        dmz = (int) (dm / (nlcx * nlcy));
        dmy = (int) ((dm % (nlcx * nlcy)) / nlcx);
        dmx = dm % nlcx;

        dnx = abs(dmx - nlcx * nint((double) dmx / (double) nlcx));

        if ((dmx == 0) || ((dmy == nlcy - 1) && (dmz == nlcz - 1)))
            dny = dmy - nlcy * nint((double) dmy / (double) nlcy);
        else
        {
            t1ny = abs(dmy - nlcy * nint((double) dmy / (double) nlcy));
            t2ny = abs((dmy + 1) - nlcy * nint((double) (dmy + 1) / (double) nlcy));
            dny = MIN(t1ny, t2ny);
        }

        if ((dmz == nlcz - 1) || ((dmx == 0) && (dmy == 0)))
            dnz = dmz - nlcz * nint((double) dmz / (double) nlcz);
        else
        {
            t1nz = abs(dmz - nlcz * nint((double) dmz / (double) nlcz));
            t2nz = abs((dmz + 1) - nlcz * nint((double) (dmz + 1) / (double) nlcz));
            dnz = MIN(t1nz, t2nz);
        }

        r2 = X2(MAX(dnx, 1) - 1) * rlcx2 + X2(MAX(dny, 1) - 1) * rlcy2 + X2(MAX(dnz, 1) - 1) * rlcz2;

        if (r2 <= cutnb2)
        {

            if (ptrmask[k] == -1)
                ptrmask[k] = dm;

            ptrmask[k + 1] = dm;

        }
        else
        {
            if (ptrmask[k] != -1)
                k += 2;
        }

    }

    for (i = 0; i < param->nAtom; i++)
    {
        xu[i] = (x[i] * box->u1 + y[i] * box->u2 + z[i] * box->u3) + 0.5;
        yu[i] = (x[i] * box->v1 + y[i] * box->v2 + z[i] * box->v3) + 0.5;
        zu[i] = (x[i] * box->w1 + y[i] * box->w2 + z[i] * box->w3) + 0.5;
    }

    for (i = 0; i < param->nAtom; i++)
    {

        ix = (int) dnlcx * xu[i];
        ix = MIN(ix, nlcx - 1);

        iy = (int) dnlcy * yu[i];
        iy = MIN(iy, nlcy - 1);

        iz = (int) dnlcz * zu[i];
        iz = MIN(iz, nlcz - 1);

        m = ix + nlcx * (iy + nlcy * iz);

        cell[ptrcell[m]] = i;
        if (ptrcell[m]<(m + 1) * atpercell)
            ptrcell[m]++;
        else
            my_error(130, __FILE__, __LINE__, 0);

    }

    k = ptrcell[0];
    ptrcell[0] = 0;
    for (m = 1; m < neigh->nCells; m++)
    {
        j = ptrcell[m];
        ptrcell[m] = k;
        for (i = m * atpercell; i < j; i++)
        {
            cell[k] = cell[i];
            k++;
        }
    }

    cell = (int*) realloc(cell, param->nAtom * sizeof (*cell));

    neigh->sizeList = 0;

    for (m = 0; m < neigh->nCells; m++)
    {
        for (i = ptrcell[m]; i < ptrcell[m + 1]; i++)
        {
            for (j = i + 1; i < ptrcell[m + 1]; i++)
            {
                exclude = 0;

                if ((frozen[i] * frozen[j]))
                {
                    exclude = 1;
                }
                else
                {
                    for (kk = 0; kk < exclPair[i]; kk++)
                    {
                        if (exclList[i][kk] == j)
                        {
                            exclude = 1;
                            break;
                        }
                    }
                }

                if (!exclude)
                {
                    neigh->sizeList++;
                }
            }

            for (s = 0; s < stripes; s++)
            {
                m1 = m + ptrmask[2 * s];
                if (m1 < neigh->nCells)
                {
                    m2 = MIN(m + ptrmask[2 * s + 1], neigh->nCells - 1);
                    for (j = ptrcell[m1]; j < ptrcell[m2 + 1]; j++)
                    {

                        delta[0] = x[j] - x[i];
                        delta[1] = y[j] - y[i];
                        delta[2] = z[j] - z[i];

                        r2 = dist(box, delta);

                        if (r2 <= cutnb2)
                        {
                            exclude = 0;

                            if ((frozen[i] * frozen[j]))
                            {
                                exclude = 1;
                            }
                            else
                            {
                                for (kk = 0; kk < exclPair[i]; kk++)
                                {
                                    if (exclList[i][kk] == j)
                                    {
                                        exclude = 1;
                                        break;
                                    }
                                }
                            }

                            if (!exclude)
                            {
                                neigh->sizeList++;
                            }
                        }
                    }
                }
            }
        }
    }

    *neighList = (int*) my_malloc(neigh->sizeList * sizeof (**neighList));

    *neighPair = (int*) my_malloc(param->nAtom * sizeof (**neighPair));

    for (i = 0; i < param->nAtom; i++)
    {
        (*neighPair)[i] = 0;
    }

    *neighOrder = (int*) my_malloc(param->nAtom * sizeof (**neighOrder));

    iorder = 0;
    offset = 0;

    for (m = 0; m < neigh->nCells; m++)
    {
        for (i = ptrcell[m]; i < ptrcell[m + 1]; i++)
        {
            (*neighOrder)[iorder++] = i;

            for (j = i + 1; i < ptrcell[m + 1]; i++)
            {
                exclude = 0;

                if ((frozen[i] * frozen[j]))
                {
                    exclude = 1;
                }
                else
                {
                    for (kk = 0; kk < exclPair[i]; kk++)
                    {
                        if (exclList[i][kk] == j)
                        {
                            exclude = 1;
                            break;
                        }
                    }
                }

                if (!exclude)
                {
                    (*neighList)[offset] = j;
                    (*neighPair)[i]++;
                    offset++;
                }
            }

            for (s = 0; s < stripes; s++)
            {
                m1 = m + ptrmask[2 * s];
                if (m1 < neigh->nCells)
                {
                    m2 = MIN(m + ptrmask[2 * s + 1], neigh->nCells - 1);
                    for (j = ptrcell[m1]; j < ptrcell[m2 + 1]; j++)
                    {
                        delta[0] = x[j] - x[i];
                        delta[1] = y[j] - y[i];
                        delta[2] = z[j] - z[i];

                        r2 = dist(box, delta);

                        if (r2 <= cutnb2)
                        {
                            exclude = 0;

                            if ((frozen[i] * frozen[j]))
                            {
                                exclude = 1;
                            }
                            else
                            {
                                for (kk = 0; kk < exclPair[i]; kk++)
                                {
                                    if (exclList[i][kk] == j)
                                    {
                                        exclude = 1;
                                        break;
                                    }
                                }
                            }

                            if (!exclude)
                            {

                                if (offset >= neigh->sizeList)
                                {
                                    fprintf(outFile, "WARNING: List larger than estimated. Size increased from %d", neigh->sizeList);
                                    neigh->sizeList = (int) (neigh->sizeList * (1. + TOLLIST)) + 1;
                                    fprintf(outFile, " to %d.\n", neigh->sizeList);

                                    *neighList = (int*) realloc(*neighList, neigh->sizeList * sizeof (**neighList));
                                }

                                (*neighList)[offset] = j;
                                (*neighPair)[i]++;
                                offset++;
                            }
                        }
                    }
                }
            }
        }
    }


}
