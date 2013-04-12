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

#ifndef TIMINGH_INCLUDED
#define TIMINGH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

#if (defined TIMING && defined __unix__ && !defined __STRICT_ANSI__)

#include <stdint.h>

enum TIMER_LOCATION{
  TIMER_ALL,
  TIMER_ENERGY_TOT,
  TIMER_ENERGY_NB,
  TIMER_ENERGY_NB14,
  TIMER_ENERGY_BOND,
  TIMER_ENERGY_ANGL,
  TIMER_ENERGY_UB,
  TIMER_ENERGY_DIHE,
  
  TIMER_VERLET_BUILD,
  TIMER_VERLET_UPDATE,
  
  TIMER_LNKCEL_BUILD,
  TIMER_LNKCEL_UPDATE,
  
  TIMER_SHAKE,
  TIMER_INTEGRATE,
  
  /* 
   * It is easy to add new type of timers here ...
   * One can also use the generic TIMER_OTHER
   */
  
  TIMER_OTHER
};

void init_timers();

void create_new_timer(enum TIMER_LOCATION loc);

void update_timer_begin(enum TIMER_LOCATION loc, const char *function_name);

void update_timer_end(enum TIMER_LOCATION loc, const char *function_name);

void print_timers();

void free_timers();

#endif

#ifdef	__cplusplus
}
#endif

#endif