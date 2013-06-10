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

#if (defined TIMING && defined __unix__ && !defined __STRICT_ANSI__)

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

#include "memory.h"
#include "timing.h"

extern FILE *outFile;

typedef struct elem{
  
  struct timespec t_ini;
  struct timespec t_end;
  
  uint64_t n_update;
  uint64_t total_time_nsec;
  
  enum TIMER_LOCATION location;
  const char *function_name;
  
  struct elem *prev;
  struct elem *next;
  
}TIMER_ELEMENT;

typedef struct{
  
  TIMER_ELEMENT *first;
  TIMER_ELEMENT *last;
  
  uint32_t nTimers;
  
}TIMER_LIST;

static TIMER_LIST *list_of_timers;

void init_timers()
{
  list_of_timers=NULL;
  list_of_timers=(TIMER_LIST*)my_malloc(sizeof(TIMER_LIST));
  
  list_of_timers->nTimers=0;
  list_of_timers->first=NULL;
  list_of_timers->last=NULL;
}

void create_new_timer(enum TIMER_LOCATION loc)
{
  TIMER_ELEMENT *new_timer=NULL;
  new_timer=(TIMER_ELEMENT*)my_malloc(sizeof(TIMER_ELEMENT));

//   clock_gettime(CLOCK_REALTIME,&(new_timer->t_ini));
  new_timer->location=loc;
  new_timer->n_update=0;
  new_timer->total_time_nsec=0;
  
  if(list_of_timers->nTimers==0)
  {
    list_of_timers->first = new_timer;
    list_of_timers->last = new_timer;
    
    new_timer->prev = NULL;
    new_timer->next = NULL;
    
    (list_of_timers->nTimers)++;
  }
  else
  {
    list_of_timers->last->next=new_timer;
    
    new_timer->prev = list_of_timers->last;
    new_timer->next = NULL;
    
    list_of_timers->last=new_timer;
    (list_of_timers->nTimers)++;
  }
}

void update_timer_begin(enum TIMER_LOCATION loc, const char *function_name)
{
  struct timespec tmp;
  clock_gettime(CLOCK_REALTIME,&tmp);
  
  TIMER_ELEMENT *ptr=list_of_timers->first;
  while(ptr->location != loc)
  {
    ptr = ptr->next;
  }
  
  ptr->t_ini = tmp;
  ptr->function_name=function_name;
}

void update_timer_end(enum TIMER_LOCATION loc, const char *function_name)
{
  struct timespec tmp;
  clock_gettime(CLOCK_REALTIME,&tmp);
  
  TIMER_ELEMENT *ptr=list_of_timers->first;
  while(ptr->location != loc)
  {
    ptr = ptr->next;
  }
  
  ptr->t_end = tmp;
  
  ptr->total_time_nsec += (ptr->t_end.tv_sec - ptr->t_ini.tv_sec)*1000000000 + (ptr->t_end.tv_nsec - ptr->t_ini.tv_nsec);
  (ptr->n_update)++;
}

void free_timers()
{
  TIMER_ELEMENT *ptr=list_of_timers->first;
  
  while(list_of_timers->nTimers!=0)
  {
    if(ptr->next==NULL)
    {
      free(ptr);
      (list_of_timers->nTimers)--;
    }
    else
    {
      ptr = ptr->next;
      free(ptr->prev);
      (list_of_timers->nTimers)--;
    }
  }
  
  free(list_of_timers);
    
}

void print_timers()
{
  TIMER_ELEMENT *ptr=list_of_timers->first;
  
  while(ptr->next != NULL)
  {
    fprintf(outFile,"Time for timer\t%30s\t:\t%lf\tseconds\n",ptr->function_name,(double)ptr->total_time_nsec/1e9);
    ptr = ptr->next;
  }
  fprintf(outFile,"Time for timer\t%30s\t:\t%lf\tseconds\n",ptr->function_name,(double)ptr->total_time_nsec/1e9);
}

#endif












