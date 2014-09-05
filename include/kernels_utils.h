#ifndef CUDA_UTILSH_INCLUDED
#define CUDA_UTILSH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif
  
# ifdef DOUBLE_CUDA

  __host__ __device__ __inline__ cplx cudaAdd(cplx za,cplx zb);
  __host__ __device__ __inline__ cplx cudaSub(cplx za,cplx zb);
  __host__ __device__ __inline__ cplx cudaMul(cplx za,cplx zb);
  __host__ __device__ __inline__ cplx cudaDiv(cplx za,cplx zb);
  __host__ __device__ __inline__ cplx cudaConj(cplx z)
  __host__ __device__ __inline__ cplx cudaExpC(cplx z);
  __host__ __device__ __inline__ cplx cudaExp(real z);
  __host__ __device__ __inline__ cplx cudaExpI(real z);
  __host__ __device__ __inline__ real cudaRe(cplx z);
  __host__ __device__ __inline__ real cudaIm(cplx z);
  __host__ __device__ __inline__ real cudaComp(real z);

#else
  
  __host__ __device__ __inline__ cplx cudaAdd(cplx za,cplx zb);
  __host__ __device__ __inline__ cplx cudaSub(cplx za,cplx zb);
  __host__ __device__ __inline__ cplx cudaMul(cplx za,cplx zb);
  __host__ __device__ __inline__ cplx cudaDiv(cplx za,cplx zb);
  __host__ __device__ __inline__ cplx cudaConj(cplx z);
  __host__ __device__ __inline__ cplx cudaExpC(cplx z);
  __host__ __device__ __inline__ cplx cudaExp(real z);
  __host__ __device__ __inline__ cplx cudaExpI(real z);
  __host__ __device__ __inline__ real cudaRe(cplx z);
  __host__ __device__ __inline__ real cudaIm(cplx z);
  __host__ __device__ __inline__ real cudaComp(real z);

#ifdef	__cplusplus
}
#endif

#endif
