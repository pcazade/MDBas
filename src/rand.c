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

/**
 * \file rand.c
 * \brief File for generating random numbers : dSFMT, an external high quality generator is used.
 * \author Pierre-Andre Cazade and Florent Hedin
 * \version alpha-branch
 * \date 2012
 */

#include <math.h>

#include "global.h"
#include "dSFMT.h"

static dsfmt_t dsfmt;

void init_rand(unsigned int seed)
{
    dsfmt_init_gen_rand(&dsfmt,seed);
}

real get_rand()
{
    return dsfmt_genrand_open_open(&dsfmt);
}

void get_BoxMuller(real *u, real *v)
{
    real a,b,s;

    do
    {

        a = 2. * get_rand() - 1.;
        b = 2. * get_rand() - 1.;
        s = a*a + b*b;

    } while (s >= 1.);

    *u= a*sqrt(-2.*log(s)/s);
    *v= b*sqrt(-2.*log(s)/s);

}
