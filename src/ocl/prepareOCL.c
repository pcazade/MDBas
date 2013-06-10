#ifdef OPENCL

#include <stdlib.h>
#include <stdio.h>

#include "ocl/ocl_top.h"

/**
 * The first action is to find an openCL platform, i.e. an implementation of the OpenCL standard by a vendor.
 * For example : 
 * 	* NVIDIA OpenCL SDK for cuda
 * 	* AMD APP SDK
 * 	* INTEL OpenCL SDK
 * 	* other SDK from various vendors (IBM, ... ? )
 * For the moment only 1 platform is selected.
 * Informations printed in a log file.
 */
void ocl_get_platform()
{
  cl_int errCode;
  char buff[4096];

  cl_uint num_platforms = 1;
  
  errCode = clGetPlatformIDs(num_platforms, &platform, NULL);
  if (errCode == CL_SUCCESS) ; else { fprintf(ocl_info,"Error when obtaining platforms : error code : %d\n",errCode);  }

  errCode = clGetPlatformInfo(platform, CL_PLATFORM_NAME, sizeof(buff), &buff, NULL);
  if (errCode == CL_SUCCESS) fprintf(ocl_info,"Platform name is : %s\n", buff) ; else { fprintf(ocl_info,"Error for Platform name  : error code %d\n", errCode);  }

  errCode = clGetPlatformInfo(platform, CL_PLATFORM_VENDOR, sizeof(buff), &buff, NULL);
  if (errCode == CL_SUCCESS) fprintf(ocl_info,"Platform vendor is : %s\n", buff) ; else { fprintf(ocl_info,"Error for Platform vendor  : error code %d\n", errCode);  }

  errCode = clGetPlatformInfo(platform, CL_PLATFORM_VERSION, sizeof(buff), &buff, NULL);
  if (errCode == CL_SUCCESS) fprintf(ocl_info,"Platform version is  : %s\n", buff) ; else { fprintf(ocl_info,"Error for Platform version : error code %d\n", errCode);  }

  errCode = clGetPlatformInfo(platform, CL_PLATFORM_PROFILE, sizeof(buff), &buff, NULL);
  if (errCode == CL_SUCCESS) fprintf(ocl_info,"Platform profile is  : %s\n", buff) ; else { fprintf(ocl_info,"Error for Platform profile : error code %d\n", errCode);  }

  errCode = clGetPlatformInfo(platform, CL_PLATFORM_EXTENSIONS, sizeof(buff), &buff, NULL);
  if (errCode == CL_SUCCESS) fprintf(ocl_info,"Platform extensions are  : %s\n\n", buff) ; else { fprintf(ocl_info,"Error for Platform extensions  : error code %d\n", errCode);  }

}

/**
 * The second action is to find devices on the selected platform.
 * For example :
 * 	* A gaming GPU (NVIDIA GEFORCE or AMD RADEON HD ...) or a computing card (NVIDIA TESLA ...)
 * 	* A compatible CPU (most of the CPUS supporting SSE2 instructions are compatible)
 * 	* other dedicated devices (cards from IBM, and others ... ?)
 * For the moment only 1 device is selected.
 * Informations printed in a log file.
 */
void ocl_get_device()
{
  cl_int errCode;
  char buff[4096];

  cl_uint num_devices = 1;
  cl_ulong mem_info;
  
  errCode = clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, num_devices, &device, NULL);
  if (errCode == CL_SUCCESS) ; else { fprintf(ocl_info,"Error when obtaining devices for Platform : error code : %d\n", errCode);  }
  
  errCode = clGetDeviceInfo(device, CL_DEVICE_NAME, sizeof(buff), &buff, NULL);
  if (errCode == CL_SUCCESS) fprintf(ocl_info,"Platform Device name is : %s\n", buff) ; else { fprintf(ocl_info,"Error for Platform Device name : error code %d\n", errCode);  }
  
  errCode = clGetDeviceInfo(device, CL_DEVICE_VENDOR, sizeof(buff), &buff, NULL);
  if (errCode == CL_SUCCESS) fprintf(ocl_info,"Platform Device vendor is : %s\n", buff) ; else { fprintf(ocl_info,"Error for Platform Device vendor : error code %d\n", errCode);  }
  
  errCode = clGetDeviceInfo(device, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(mem_info), &mem_info, NULL);
  if (errCode == CL_SUCCESS) fprintf(ocl_info,"Platform Device global memory is (MB) : %lf\n", mem_info/TO_MB) ; else { fprintf(ocl_info,"Error for Platform Device global memory : error code %d\n", errCode);  }
  
  errCode = clGetDeviceInfo(device, CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(mem_info), &mem_info, NULL);
  if (errCode == CL_SUCCESS) fprintf(ocl_info,"Platform Device max alloc memory is (MB) : %lf\n", mem_info/TO_MB) ; else { fprintf(ocl_info,"Error for Platform Device max alloc memory : error code %d\n", errCode);  }

  errCode = clGetDeviceInfo(device, CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, sizeof(mem_info), &mem_info, NULL);
  if (errCode == CL_SUCCESS) fprintf(ocl_info,"Platform Device constant memory is (KB) : %lf\n", mem_info/TO_KB) ; else { fprintf(ocl_info,"Error for Platform Device constant memory  : error code %d\n", errCode);  }

  errCode = clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(mem_info), &mem_info, NULL);
  if (errCode == CL_SUCCESS) fprintf(ocl_info,"Platform Device local memory is (KB) : %lf\n", mem_info/TO_KB) ; else { fprintf(ocl_info,"Error for Platform Device local memory  : error code %d\n", errCode);  }

  errCode = clGetDeviceInfo(device, CL_DEVICE_EXTENSIONS, sizeof(buff), &buff, NULL);
  if (errCode == CL_SUCCESS) fprintf(ocl_info,"Platform Device extensions are : %s\n\n", buff) ; else { fprintf(ocl_info,"Error for Platform Device extensions  : error code %d\n", errCode);  }

}

/**
 * Creates a context : a context is an ensemble on devices on which calculations occurs.
 * For the moment we create a basic context made of only 1 device selected previously
 */
void ocl_get_context()
{
    cl_int errCode;
    
    context = clCreateContext(NULL, 1, &device, NULL, NULL, &errCode);
    if (errCode == CL_SUCCESS) ; else { fprintf(ocl_info,"Error when adding a device to context : error code : %d\n", errCode); }
}

/**
 * Creates a command queue : it manages all the instructions sent to the context and its devices
 */
void ocl_get_cmdQueue()
{
  cl_int errCode;
  
  command_queue = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &errCode);
  if (errCode == CL_SUCCESS) ; else { fprintf(ocl_info,"Error when creating command queue : error code : %d\n", errCode); }

}

#endif

