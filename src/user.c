#include <stdio.h>
#include <stdlib.h>

// dynamic loading of code
#ifdef __unix__
#include <dlfcn.h>
#endif

#include "global.h"
#include "errors.h"
#include "user.h"

//extern FILE *outFile;

static void *plugin;
static void *fnc;

UserEnergyPtr loadUserPlugin(const char pluginName[], const char funcName[])
{
    plugin = NULL;
    fnc = NULL;
    
#ifdef __unix__
    plugin = dlopen(pluginName, RTLD_LAZY);

    if(plugin==NULL)
        my_error(PLUGINS_DLOPEN_ERROR,__FILE__,__LINE__,2,__func__,pluginName);
    
    fnc = dlsym(plugin,funcName);
    
    if(fnc==NULL)
        my_error(PLUGINS_DLSYM_ERROR,__FILE__,__LINE__,3,__func__,pluginName,funcName);
#endif

    UserEnergyPtr eneCustom = (UserEnergyPtr) fnc;
    
    return eneCustom;
    
}

void closeUserPlugin(void)
{
    int err = dlclose(plugin);
}

