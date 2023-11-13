#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <CL/opencl.h>
#include <CL/cl_ext.h>
#include "opencl_lib.h"
#include "params.h"


//TODO Probably not needed to be defined and should be pruned
#define CL_HPP_CL_1_2_DEFAULT_BUILD
#define CL_HPP_TARGET_OPENCL_VERSION 120
#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#define CL_HPP_ENABLE_PROGRAM_CONSTRUCTION_FROM_ARRAY_COMPATIBILITY 1
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS

cl_uint load_file_to_memory(const char *filename, char **result)
{
    cl_uint size = 0;
    FILE *f = fopen(filename, "rb");
    if (f == NULL) {
        *result = NULL;
        return -1; // -1 means file opening fail
    }
    fseek(f, 0, SEEK_END);
    size = ftell(f);
    fseek(f, 0, SEEK_SET);
    *result = (char *)malloc(size+1);
    if (size != fread(*result, sizeof(char), size, f)) {
        free(*result);
        return -2; // -2 means file reading fail
    }
    fclose(f);
    (*result)[size] = 0;
    return size;
}


void platform_init(cl_program *program, cl_context *context, cl_command_queue *commands, char *xclbin){

 
    cl_int err;
    cl_platform_id platform_id;
    cl_device_id device_id;
	const char* target_device_name;
	const char* zcu_target_device_name = "zcu102_base";
	const char* alveo_target_device_name = "xilinx_u280_xdma_201920_3";
    char cl_platform_vendor[1001];

	#ifdef ALVEO
    {
		target_device_name = alveo_target_device_name;
    }
	#endif
	#ifdef ZCU
	{
		target_device_name = zcu_target_device_name;
	}
	#endif

    cl_platform_id platforms[16];
    cl_uint platform_count;
    cl_uint platform_found = 0;
    cl_uint num_devices;
    cl_uint device_found = 0;
    cl_device_id devices[16];
    char cl_device_name[1001];
    cl_int status;

    //OpenCL Initialization code STARTS
	// ------------------------------------------------------------------------------------
	// Step 1: Get All PLATFORMS, then search for Target_Platform_Vendor (CL_PLATFORM_VENDOR)
	// ------------------------------------------------------------------------------------

	// Get the number of platforms
	// ..................................................

    err = clGetPlatformIDs(16, platforms, &platform_count);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		printf("Error: Failed to find an OpenCL platform!\n");
		printf("Test failed\n");
		return EXIT_FAILURE;
	}
	printf("INFO: Found %d platforms\n", platform_count);
	#endif



	  // ....................................................................................
	  // step 1:  Search for Platform (ex: Xilinx) using: CL_PLATFORM_VENDOR = Target_Platform_Vendor
	  // Check if the current platform matches Target_Platform_Vendor
	  // ................device_....................................................................

	for (cl_uint iplat=0; iplat<platform_count; iplat++) {
		err = clGetPlatformInfo(platforms[iplat], CL_PLATFORM_VENDOR, 1000, (void *)cl_platform_vendor,NULL);
		#ifdef OCL_API_DEBUG
		if (err != CL_SUCCESS) {
			printf("Error: clGetPlatformInfo(CL_PLATFORM_VENDOR) failed!\n");
			printf("Test failed\n");
			return EXIT_FAILURE;
		}
		if (strcmp(cl_platform_vendor, "Xilinx") == 0) {
			printf("INFO: Selected platform %d from %s\n", iplat, cl_platform_vendor);
			platform_id = platforms[iplat];
			platform_found = 1;
		}
		#endif
	}
	#ifdef OCL_API_DEBUG
	if (!platform_found) {
		printf("ERROR: Platform Xilinx not found. Exit.\n");
		return EXIT_FAILURE;
	}
	#endif

	   // ------------------------------------------------------------------------------------
	   // Step 1:  Get All Devices for selected platform Target_Platform_ID
	   //            then search for Xilinx platform (CL_DEVICE_TYPE_ACCELERATOR = Target_Device_Name)
	   // ------------------------------------------------------------------------------------


    err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_ACCELERATOR, 16, devices, &num_devices);
	#ifdef OCL_API_DEBUG
	printf("INFO: Found %d devices\n", num_devices);
	if (err != CL_SUCCESS) {
		printf("ERROR: Failed to create a device group!\n");
		printf("ERROR: Test failed\n");
		return -1;
	}
	#endif
	 // ------------------------------------------------------------------------------------
	 // Step 1:  Search for CL_DEVICE_NAME = Target_Device_Name
	 // ............................................................................

   for (cl_uint i=0; i<num_devices; i++) {
		err = clGetDeviceInfo(devices[i], CL_DEVICE_NAME, 1024, cl_device_name, 0);
		#ifdef OCL_API_DEBUG
		if (err != CL_SUCCESS) {
			printf("Error: Failed to get device name for device %d!\n", i);
			printf("Test failed\n");
			return EXIT_FAILURE;
		}
		printf("CL_DEVICE_NAME %s\n", cl_device_name);
		#endif



	   // ............................................................................
	   // Step 1: Check if the current device matches Target_Device_Name
	   // ............................................................................

	   if(strcmp(cl_device_name, target_device_name) == 0) {
			device_id = devices[i];
			device_found = 1;
			#ifdef OCL_API_DEBUG
			printf("Selected %s as the target device\n", cl_device_name);
			#endif
	   }
	}


	// ------------------------------------------------------------------------------------
	// Step 1: Create Context
	// ------------------------------------------------------------------------------------
	*context = clCreateContext(0, 1, &device_id, NULL, NULL, &err);
	#ifdef OCL_API_DEBUG
	if (!(*context)) {
		printf("Error: Failed to create a compute context!\n");
		printf("Test failed\n");
		return EXIT_FAILURE;
	}
	#endif

	// ------------------------------------------------------------------------------------
	// Step 1: Create Command Queue
	// ------------------------------------------------------------------------------------
	*commands = clCreateCommandQueue(*context, device_id, CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE | CL_QUEUE_PROFILING_ENABLE, &err);
//	commands = clCreateCommandQueue(context, device_id, NULL, &err);
	#ifdef OCL_API_DEBUG
	if (!(*commands)) {
		printf("Error: Failed to create a command commands!\n");
		printf("Error: code %i\n",err);
		printf("Test failed\n");
		return EXIT_FAILURE;
	}
	#endif



   unsigned char *kernelbinary;
   //TODO that has been given as a cml line argument. Dele if it works!
//    char *xclbin = argv[1];

   // ------------------------------------------------------------------
   // Step 1: Load Binary File from a disk to Memory
   // ------------------------------------------------------------------
	#ifdef OCL_API_DEBUG
    printf("INFO: loading xclbin %s\n", xclbin);
	#endif
    //TODO Should I dereference the xclbin here?
    cl_uint n_i0 = load_file_to_memory(xclbin, (char **) &kernelbinary);
	#ifdef OCL_API_DEBUG
    if (n_i0 < 0) {
	    printf("failed to load kernel from xclbin: %s\n", xclbin);
	    printf("Test failed\n");
	    return EXIT_FAILURE;
    }
	#endif



	// ------------------------------------------------------------
	// Step 1: Create a program using a Binary File
	// ------------------------------------------------------------
	size_t n0 = n_i0;
	*program = clCreateProgramWithBinary(*context, 1, &device_id, &n0, (const unsigned char **) &kernelbinary, &status, &err);
	free(kernelbinary);

	// ============================================================================
	// Step 2: Create Program and Kernels
	// ============================================================================
	//   o) Build a Program from a Binary File
	//   o) Create Kernels
	// ============================================================================
	#ifdef OCL_API_DEBUG
	if ((!*program) || (err!=CL_SUCCESS)) {
		printf("Error: Failed to create compute program from binary %d!\n", err);
		printf("Test failed\n");
		return EXIT_FAILURE;
	}
	#endif

	// -------------------------------------------------------------
	// Step 2: Build (compiles and links) a program executable from binary
	// -------------------------------------------------------------
	err = clBuildProgram(*program, 0, NULL, NULL, NULL, NULL);
	#ifdef OCL_API_DEBUG
	if (err != CL_SUCCESS) {
		size_t len;
		char buffer[2048];

		printf("Error: Failed to build program executable!\n");
		clGetProgramBuildInfo(*program, device_id, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
		printf("%s\n", buffer);
		printf("Test failed\n");
		return EXIT_FAILURE;
	}
	#endif

	return;
}
