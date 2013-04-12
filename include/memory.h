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

#ifndef MEMORY_H_INCLUDED
#define MEMORY_H_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

#include <stdlib.h>

void* my_malloc(size_t size);
void* my_realloc(void *ptr1,size_t size);
void* my_calloc(int dim, size_t size);

void** calloc_2D(int dim1, int dim2, size_t si);
void free_2D(int dim1, ...);

void*** calloc_3D(int dim1, int dim2, int dim3, size_t si);
void free_3D(int dim1, int dim2, ...);

#ifdef	__cplusplus
}
#endif

#endif // MEMORY_H_INCLUDED
