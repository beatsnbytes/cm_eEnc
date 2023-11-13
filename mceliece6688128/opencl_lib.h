#include <CL/opencl.h>
#include <CL/cl_ext.h>

cl_uint load_file_to_memory(const char *, char **);
void platform_init(cl_program *, cl_context *, cl_command_queue *, char *);

