/*
 * Copyright (c) 2013 Pierre-Andre Cazade
 * Copyright (c) 2013 Florent hedin
 * 
 * This file is part of MDBas.
 *
 * MDBas is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MDBas is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MDBas.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * \file memory.c
 * \brief Contains functions for managing allocations and freeing of memory dynamically.
 * \author Pierre-Andre Cazade and Florent Hedin
 * \version alpha-branch
 * \date 2012
 */

#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <stdarg.h>

#include "memory.h"
#include "errors.h"

//totally allocated memory
//static size_t totSize = 0;

// redefinitions of malloc and calloc with failure check ; totSize used for keeping trace of totally allocated memory
void* my_malloc(size_t size)
{
  void *ptr = NULL;  
  ptr=(void*)malloc(size);
  
  if(ptr==NULL)
    my_error(MEM_ALLOC_ERROR,__FILE__,__LINE__,0);
  
  //totSize += size;
  
  return ptr;
}

void* my_realloc(void *ptr1,size_t size)
{
  void *ptr2 = NULL;  
  ptr2=(void*)realloc(ptr1,size);
  
  if(ptr2==NULL)
    my_error(MEM_ALLOC_ERROR,__FILE__,__LINE__,0);
  
  //totSize += size;
  
  return ptr2;
}

void* my_calloc(int dim, size_t size)
{
  void *ptr = NULL;  
  ptr=(void*)calloc(dim,size);
  
  if(ptr==NULL)
    my_error(MEM_ALLOC_ERROR,__FILE__,__LINE__,0);
  
  //totSize += dim*size;
  
  return ptr;
}

//for allocating an array of dimensions dim1*dim2 and of bytes size si
void** calloc_2D(int dim1, int dim2, size_t si)
{
    int i;
    void **array=NULL;

    //Here, allocation of the first dimension : rows
    array = my_calloc(dim1,si);

    //Then allocation of the second dimensions : columns
    for (i=0; i<dim1; i++)
        array[i] = my_calloc(dim2,si);

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

    array=my_calloc(dim1,si);

    for(i=0; i<dim1; i++)
    {
        array[i]=my_calloc(dim2,si);

        for(j=0; j<dim2; j++)
            array[i][j]=my_calloc(dim3,si);
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