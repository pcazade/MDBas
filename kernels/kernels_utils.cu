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

# ifdef DOUBLE_CUDA

__host__ __device__ __inline__ cuDoubleComplex cudaAdd(cuDoubleComplex za,cuDoubleComplex zb)
{
  return cuCadd(za,zb);
}

__host__ __device__ __inline__ cuDoubleComplex cudaSub(cuDoubleComplex za,cuDoubleComplex zb)
{
  return cuCsub(za,zb);
}

__host__ __device__ __inline__ cuDoubleComplex cudaMul(cuDoubleComplex za,cuDoubleComplex zb)
{
  return cuCmul(za,zb);
}

__host__ __device__ __inline__ cuDoubleComplex cudaDiv(cuDoubleComplex za,cuDoubleComplex zb)
{
  return cuCdiv(za,zb);
}

__host__ __device__ __inline__ cuDoubleComplex cudaConj(cuDoubleComplex z)
{
  return cuConj(z);
}

__host__ __device__ __inline__ cuDoubleComplex cudaExpC(cuDoubleComplex z)
{

  cuDoubleComplex res;

  double t = exp(z.x);

  res.x=cos(z.y);
  res.y=sin(z.y);

  res.x *= t;
  res.y *= t;

  return res;

}

__host__ __device__ __inline__ cuDoubleComplex cudaExp(double z)
{

  cuDoubleComplex res;

  res.x=exp(z);
  res.y=0.0;

  return res;

}

__host__ __device__ __inline__ cuDoubleComplex cudaExpI(double z)
{

  cuDoubleComplex res;

  res.x=cos(z);
  res.y=sin(z);

  return res;

}

#else

__host__ __device__ __inline__ cuComplex cudaAdd(cuComplex za,cuComplex zb)
{
  return cuCaddf(za,zb);
}

__host__ __device__ __inline__ cuComplex cudaSub(cuComplex za,cuComplex zb)
{
  return cuCsubf(za,zb);
}

__host__ __device__ __inline__ cuComplex cudaMul(cuComplex za,cuComplex zb)
{
  return cuCmulf(za,zb);
}

__host__ __device__ __inline__ cuComplex cudaDiv(cuComplex za,cuComplex zb)
{
  return cuCdivf(za,zb);
}

__host__ __device__ __inline__ cuComplex cudaConj(cuComplex z)
{
  return cuConjf(z);
}

__host__ __device__ __inline__ cuComplex cudaExpC(cuComplex z)
{

  cuComplex res;

  float t = expf(z.x);

  res.x=cosf(z.y);
  res.y=sinf(z.y);

  res.x *= t;
  res.y *= t;

  return res;

}

__host__ __device__ __inline__ cuComplex cudaExp(float z)
{

  cuComplex res;

  res.x=expf(z);
  res.y=0.f;

  return res;

}

__host__ __device__ __inline__ cuComplex cudaExpI(float z)
{

  cuComplex res;

  res.x=cosf(z);
  res.y=sinf(z);

  return res;

}

#endif