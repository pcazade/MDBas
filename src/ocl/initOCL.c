#ifdef OPENCL

#include <stdlib.h>
#include <stdio.h>

#include "ocl/ocl_top.h"

//file for info log
FILE *ocl_info;

//common CL_STRUCTURES
cl_platform_id 		platform;
cl_device_id 		device;
cl_context		context;
cl_command_queue	command_queue;

// a call to this function will initialise anything required by openCL
void ocl_init_all()
{
  ocl_info = fopen("openCL.log","wt");

  // Prepare ocl parameters (cf prepareOCL.c for informations)
  ocl_get_platform();
  ocl_get_device();
  ocl_get_context();
  ocl_get_cmdQueue();

}

//close files and free memory and more ...
void ocl_end_all()
{
  fclose(ocl_info);
  
  clFlush(command_queue);
  clFinish(command_queue);
    
  clReleaseCommandQueue(command_queue);
  clReleaseContext(context);
}

#endif

