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
 * \file cuda_utils.c
 * \brief Contains various and general utilitary functions.
 * \author Pierre-Andre Cazade
 * \version alpha-branch
 * \date 2014
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <float.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <cuComplex.h>

#include "global.h"

# ifdef DOUBLE_CUDA

__host__ __device__ __inline__ cplx cudaAdd(cplx za,cplx zb)
{
  return cuCadd(za,zb);
}

__host__ __device__ __inline__ cplx cudaSub(cplx za,cplx zb)
{
  return cuCsub(za,zb);
}

__host__ __device__ __inline__ cplx cudaMul(cplx za,cplx zb)
{
  return cuCmul(za,zb);
}

__host__ __device__ __inline__ cplx cudaDiv(cplx za,cplx zb)
{
  return cuCdiv(za,zb);
}

__host__ __device__ __inline__ cplx cudaConj(cplx z)
{
  return cuConj(z);
}

__host__ __device__ __inline__ cplx cudaExpC(cplx z)
{

  cplx res;

  real t = exp(z.x);

  res.x=cos(z.y);
  res.y=sin(z.y);

  res.x *= t;
  res.y *= t;

  return res;

}

__host__ __device__ __inline__ cplx cudaExp(real z)
{

  cplx res;

  res.x=exp(z);
  res.y=0.0;

  return res;

}

__host__ __device__ __inline__ cplx cudaExpI(real z)
{

  cplx res;

  res.x=cos(z);
  res.y=sin(z);

  return res;

}

__host__ __device__ __inline__ real cudaRe(cplx z)
{
  return z.x;
}

__host__ __device__ __inline__ real cudaIm(cplx z)
{
  return z.y;
}

__host__ __device__ __inline__ real cudaComp(real z)
{
  
  cplx res;
  
  res.x=z;
  res.y=0.0;
  
  return res;
  
}

#else

__host__ __device__ __inline__ cplx cudaAdd(cplx za,cplx zb)
{
  return cuCaddf(za,zb);
}

__host__ __device__ __inline__ cplx cudaSub(cplx za,cplx zb)
{
  return cuCsubf(za,zb);
}

__host__ __device__ __inline__ cplx cudaMul(cplx za,cplx zb)
{
  return cuCmulf(za,zb);
}

__host__ __device__ __inline__ cplx cudaDiv(cplx za,cplx zb)
{
  return cuCdivf(za,zb);
}

__host__ __device__ __inline__ cplx cudaConj(cplx z)
{
  return cuConjf(z);
}

__host__ __device__ __inline__ cplx cudaExpC(cplx z)
{

  cplx res;

  float t = expf(z.x);

  res.x=cosf(z.y);
  res.y=sinf(z.y);

  res.x *= t;
  res.y *= t;

  return res;

}

__host__ __device__ __inline__ cplx cudaExp(float z)
{

  cplx res;

  res.x=expf(z);
  res.y=0.f;

  return res;

}

__host__ __device__ __inline__ cplx cudaExpI(float z)
{

  cplx res;

  res.x=cosf(z);
  res.y=sinf(z);

  return res;

}

__host__ __device__ __inline__ real cudaRe(cplx z)
{
  return z.x;
}

__host__ __device__ __inline__ real cudaIm(cplx z)
{
  return z.y;
}

__host__ __device__ __inline__ real cudaComp(real z)
{
  
  cplx res;
  
  res.x=z;
  res.y=0.0;
  
  return res;
  
}

#endif