/**
 * \file list.h
 * \brief Prototypes for file list.c
 * \author Pierre-Andre Cazade and Florent Hedin
 * \version alpha-branch
 * \date 2012
 */

#ifndef LISTH_INCLUDED
#define LISTH_INCLUDED

void makelist(SIMULPARAMS *simulCond,ATOM atom[],FORCEFIELD *ff,CONSTRAINT *constList,PBC *box);
void exclude_list(SIMULPARAMS *simulCond,ATOM atom[],FORCEFIELD *ff,CONSTRAINT *constList);
void verlet_list(SIMULPARAMS *simulCond,ATOM atom[],FORCEFIELD *ff,PBC *box);
void verlet_list_update(SIMULPARAMS *simulCond,ATOM atom[],FORCEFIELD *ff,PBC *box);
void link_cell_exclude_list(SIMULPARAMS *simulCond,ATOM atom[],FORCEFIELD *ff,CONSTRAINT *constList);
void link_cell_verlet_list(SIMULPARAMS *simulCond,ATOM atom[],FORCEFIELD *ff,PBC *box);
void link_cell_verlet_list_update(SIMULPARAMS *simulCond,ATOM atom[],FORCEFIELD *ff,PBC *box);

#endif