
#ifdef	__cplusplus
extern "C" {
#endif
    
#ifdef OPENCL

#ifndef OCL_INCLUDED
#define OCL_INCLUDED

#include <CL/cl.h>

#define TO_MB	1048576.0
#define TO_KB	1024.0

//file for info log
extern FILE *ocl_info;

//common CL_STRUCTURES
extern cl_platform_id 	platform;
extern cl_device_id 	device;
extern cl_context	context;
extern cl_command_queue	command_queue;

//in initOCL.c
void ocl_init_all();
void ocl_end_all();

//in prepareOCL.c
void ocl_get_platform();
void ocl_get_device();
void ocl_get_context();
void ocl_get_cmdQueue();

#endif

#ifdef	__cplusplus
}
#endif

#endif

