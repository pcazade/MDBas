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
