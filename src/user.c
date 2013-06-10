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

