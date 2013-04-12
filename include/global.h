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
 * \file global.h
 * \brief The main header file containing most of the defines, plus the patterns of enum and structs.
 * \author Pierre-Andre Cazade and Florent Hedin
 * \version alpha-branch
 * \date 2012
 */

#ifndef GLOBALH_INCLUDED
#define GLOBALH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

#include <stdio.h>

#define elemchg (1.602176565e-19)       /*!< The elementary charge in C (coulomb). */
#define angstr  (1.e-10)                /*!< SI definition of one Angstroem in m (meters). */
#define calory  (4.184)                 /*!< SI definition of one calory in J (joules). */
#define kcaltoiu (418.4)                /*!< Conversion factor : internal units (10 Joules) to kcal (kilo-calory). */
#define clight  (299792458.)            /*!< The speed of light (c) in \f$ m.s^{-1} \f$. */
#define NA	(6.02214129e+23)        /*!< The Avogadro constant in \f$ mol^{-1} \f$. */
#define bartoiu	(6.02214129e-3)         /*!< For converting pressure from bar \f$ 10^5.J.m^{-3} \f$ to internal units \f$ 10 J.mol^{-1}.A^{-3} \f$. */
#define kboltz  (1.3806488e-23)         /*!< The Boltzmann Constant in \f$ J.K^{-1} \f$. */
#define rboltz  (8.3144621)             /*!< The Gas constant in \f$ J.K^{-1}.mol^{-1} \f$ : it is kboltz divided by NA. */
#define rboltzui  (0.83144621)          /*!< The Gas constant in internal units,  \f$ 10 J.K^{-1}.mol^{-1} \f$. */
#define mu0     (1.e-7)                 /*!< The vacuum permeability in \f$ 4\pi  \frac{V·s}{A·m}  \f$. */
#define chgcharmm (332.0716)            /*!< For converting electrostatic energy from internal units to \f$ kcal.mol^{-1} \f$ for CHARMM. */
#define chgnamd   (332.0636)            /*!< For converting electrostatic energy from internal units to \f$ kcal.mol^{-1} \f$ for NAMD. */
#define chgdlpolyiu (138935.4835)       /*!< For converting electrostatic energy from internal units to \f$ kcal.mol^{-1} \f$ for DL_POLY. */
#define sq6rt2  (1.122462048309373)     /*!< This constant is \f$ \sqrt[6]{2} \f$. */
#define PI      (3.141592653589793)     /*!< This is \f$ \pi \f$. */
#define TWOPI   (6.283185307179586)     /*!< This is \f$ 2\pi \f$. */
#define SQRTPI  (1.772453850905516)     /*!< This is \f$ \sqrt\pi \f$. */
#define watercomp (0.007372)		/*!< This the compressibility of water in internal units. */

#define MAXLIST 2048                    /*!< Maximal number of pairs per atom for neighbours list. */
#define TOLLIST 0.01

#define FINAMELEN 512                   /*!< Default number of chracters (string length) for files name. */

#define X2(x) ((x)*(x))                 /*!< Raise any value at power 2. */
#define X3(x) (X2(x)*(x))               /*!< Raise any value at power 3. */
#define X4(x) (X2(x)*X2(x))             /*!< Raise any value at power 4. */
#define X6(x) (X3(x)*X3(x))             /*!< Raise any value at power 6. */
#define X12(x) (X6(x)*X6(x))            /*!< Raise any value at power 12. */

#define MAX(x,y) ((x)>=(y)?(x):(y))     /*!< Macro providing the maximum of (x,y). */
#define MIN(x,y) ((x)<=(y)?(x):(y))     /*!< Macro providing the minimum of (x,y). */

//#define COSDIH  1                       /*!< Cosine potential type for dihedral and improper angles. */
//#define HARMDIH 2                       /*!< Harmonic potential type for dihedral and improper angles. */

/*!
 * \enum ELEC_TYPE
 * \brief Enumeration of the possible ways of evaluating electrostatic potential and force.
 */
enum ELEC_TYPE{
  NOELEC = 0,       /*!< No electrostatic evaluation. */
  FULL   = 1,       /*!< Full electrostatic evaluation. */
  SHIFT1 = 2,       /*!< Electrostatic evaluation with SHIFT_1 cutoff. */
  SHIFT2 = 3,       /*!< Electrostatic evaluation with SHIFT_2 cutoff. */
  SWITCH = 4        /*!< Electrostatic evaluation with SWITCH cutoff. */
};

/*!
 * \enum VDW_TYPE
 * \brief Enumeration of the possible ways of evaluating Van der Waals potential and force.
 */
enum VDW_TYPE{
  NOVDW   = 0,      /*!< No Van der Waals evaluation. */
  VFULL   = 1,      /*!< Full Van der Waals evaluation. */
  VSWITCH = 2       /*!< Van der Waals evaluation with SWITCH cutoff. */
};

/*!
 * \enum BOX_TYPE
 * \brief Enumeration of the available simulation boxes.
 */
enum BOX_TYPE{
  NOBOX = 0,        /*!< No box for this simulation, so no Periodic Boundaries Conditions applied. */
  CUBIC = 1,        /*!< A cubic box is used. */
  ORBIC = 2,        /*!< An orthorhombic box is used. */
  TCLIN = 3         /*!< A triclinic box is used. */
};

/*!
 * \enum INTEGRATOR_TYPE
 * \brief Enumeration of the different integrators available.
 */
enum INTEGRATOR_TYPE{
  LEAPFROG = 0,
  VELOCITY = 1
};

/*!
 * \enum ENSEMBLE_TYPE
 * \brief Enumeration of the available emsembles.
 */
enum ENSEMBLE_TYPE{
  NVE = 0,
  NVT_B = 1,
  NPT_B = 2,
  NVT_H = 3,
  NPT_H = 4
};

/*!
 * \enum BOND_TYPE
 * \brief Enumeration of the possible bond potentials.
 */
enum BOND_TYPE{
  BHARM  = 0,       /*!< The bond is described by a harmonic potential: \f$ E=k\left(r-r_0\right)^2 \f$ */
  BMORSE = 1,       /*!< The bond is described by a Morse potential: \f$ E=k\left[1-\exp\left(-\beta\left(r-r_0\right)\right)\right]^2 \f$ */
};

/*!
 * \enum ANGLE_TYPE
 * \brief Enumeration of the possible angle potentials.
 */
enum ANGLE_TYPE{
  AHARM = 0,       /*!< The angle is described by a harmonic potential: \f$ E=k\left(\theta-\theta_0\right)^2 \f$ */
};

/*!
 * \enum DIHE_TYPE
 * \brief Enumeration of the possible dihedral angle potentials.
 */
enum DIHE_TYPE{
  DCOS  = 1,       /*!< The dihedral is described by a cosine potential: \f$ E=k\left[1+\cos\left(m*\phi-phi_0\right)\right] \f$ */
  DHARM = 2,       /*!< The dihedral is described by a harmonic potential: \f$ E=k\left(\phi-\phi_0\right)^2 \f$ */
};

/*!
 * \struct ATOM
 * \brief This structure represents an atom of the simulation and all its properties.
 *
 * First, labelling informations as 3 char[] (atom label, segment label, residue label) and 3 int (atom type, residue number, number of constraints this atom is involved in).
 *
 * Then, spatial informations : coordinates, speeds, forces associated to this atom.
 *
 * Finally, the mass and charge.
 */
typedef struct
{
  char label[5],segn[5],resn[5];
  int type,resi,ires;
}ATOM;

/*!
 * \struct PBC
 * \brief This structure contains parameters describing the simulation box.
 */
typedef struct
{
  enum BOX_TYPE type; /*<  */
  double a,a1,a2,a3,b,b1,b2,b3,c,c1,c2,c3;
  double u,u1,u2,u3,v,v1,v2,v3,w,w1,w2,w3;
  double pa,pb,pc,det,vol,vol0;
  double stress1,stress2,stress3;
  double stress4,stress5,stress6;
  double stress7,stress8,stress9;
}PBC;

typedef struct
{
  FILE *confFile,*forfFile,*propFile,*restFile,*rconFile,*trajFile;
  FILE *simuFile;
  char confName[FINAMELEN],forfName[FINAMELEN],propName[FINAMELEN];
  char restName[FINAMELEN],rconName[FINAMELEN],trajName[FINAMELEN];
  char simuName[FINAMELEN];
}IO;

typedef struct
{
  int newjob;
  
  int keyRand,keyTraj,keyProp,keyForF,keyRest;
  int printOut,printProp,printTraj,printRest;
  
  int keyMd,mdType;
  int seed;
  
  int keyMinim;
  int keyLink,noLink;
  
  int keyConstH;
  int keyNb14,keyNumForce;
  
  int keyEwald,keyAlpha,keyMmax;
  
  enum ELEC_TYPE elecType;
  enum VDW_TYPE  vdwType;
  enum INTEGRATOR_TYPE integrator;
  enum ENSEMBLE_TYPE ens;
}CTRL;

typedef struct
{
  int step,nSteps;
  double timeStep,rTimeStep;
  
  int nAtom,nDegFree,nFrozen;
  
  int nBond,nAngle,nDihedral,nImproper,nUb;
  
  int nConst,maxCycle;
  double tolShake;
  
  double chargeConst,delr,scal14;
  double cutOff,rcutOff,cutOff2,rcutOff2;
  double cutOn,cutOn2,switch2;
  
  double temp0,press0,kinTemp0;
  
  double tolMinim,maxminsiz,maxminst;
}PARAM;

typedef struct
{
  double tauT,chiT;
  
  double tauP,chiP,compress;
}BATH;

typedef struct
{
  enum BOND_TYPE type;
  int a,b;
  double k,r0,beta;
}BOND;

typedef struct
{
  enum ANGLE_TYPE type;
  int a,b,c;
  double k,theta0;
}ANGLE;

typedef struct
{
  enum DIHE_TYPE type;
  int order;
  int a,b,c,d;
  double k,phi0,mult;
}DIHE;

typedef struct
{
  int update;
  int linkRatio;
  int nCells;
  int sizeList;
  int nPair14;
}NEIGH;

/*!
 * \struct ENERGY
 * \brief This structure contains the total energy but also all its components.
 */
typedef struct
{
  double tot,pot,kin;
  double elec,vdw;
  double bond,ang,ub,dihe,impr;
  double conint,consv;
  double virelec,virvdw,virbond,virub;
  double virshake,virpot,virtot;
}ENERGY;

/*!
 * \struct CONSTRAINT
 * \brief This structure is used for storing parameters of the constraints.
 */
typedef struct
{
  int a,b;
  double rc2;
}CONSTRAINT;

/*!
 * \struct DELTA
 * \brief This structure represents a distance vector between two atoms.
 */
typedef struct
{
  double x,y,z;
}DELTA;

typedef struct
{
  int nbsp;
  int mmax,m1max,m2max,m3max;
  double prec,tol,tol1,alpha;
}EWALD;

#ifdef	__cplusplus
}
#endif

#endif
