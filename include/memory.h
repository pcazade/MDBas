/**
 * \file memory.h
 * \brief Prototypes for file memory.c
 * \author Pierre-Andre Cazade and Florent Hedin
 * \version alpha-branch
 * \date 2012
 */

#ifndef MEMORY_H_INCLUDED
#define MEMORY_H_INCLUDED

#include <stdlib.h>

void** calloc_2D(int dim1, int dim2, size_t si);
void free_2D(int dim1, ...);

void*** calloc_3D(int dim1, int dim2, int dim3, size_t si);
void free_3D(int dim1, int dim2, ...);


#endif // MEMORY_H_INCLUDED
