/**
 * \file global.h
 * \brief The main header file containing most of the defines, plus the patterns of enum and structs.
 * \author Pierre-Andre Cazade and Florent Hedin
 * \version alpha-branch
 * \date 2012
 */

#ifndef GLOBALH_INCLUDED
#define GLOBALH_INCLUDED

#define elemchg (1.602176565e-19)       /*!< The elementary charge in C (coulomb). */
#define angstr  (1.e-10)                /*!< SI definition of one Angstroem in m (meters). */
#define calory  (4.184)                 /*!< SI definition of one calory in J (joules). */
#define kcaltoiu (418.4)                /*!< Conversion factor : internal units (10 Joules) to kcal (kilo-calory). */
#define clight  (299792458.)            /*!< The speed of light (c) in \f$ m.s^{-1} \f$. */
#define NA	(6.02214129e+23)        /*!< The Avogadro constant in \f$ mol^{-1} \f$. */
#define bartoiu	(6.02214129e+3)         /*!< For converting pressure from bar \f$ 10^5.J.m^{-3} \f$ to internal units \f$ 10 J.mol^{-1}.A^{-3} \f$. */
#define kboltz  (1.3806488e-23)         /*!< The Boltzmann Constant in \f$ J.K^{-1} \f$. */
#define rboltz  (8.3144621)             /*!< The Gas constant in \f$ J.K^{-1}.mol^{-1} \f$ : it is kboltz divided by NA. */
#define rboltzui  (0.83144621)          /*!< The Gas constant in internal units,  \f$ 10 J.K^{-1}.mol^{-1} \f$. */
#define mu0     (1.e-7)                 /*!< The vacuum permeability in \f$ 4\pi  \frac{V·s}{A·m}  \f$. */
#define chgcharmm (332.0716)            /*!< For converting electrostatic energy from internal units to \f$ kcal.mol^{-1} \f$ for CHARMM. */
#define chgnamd   (332.0636)            /*!< For converting electrostatic energy from internal units to \f$ kcal.mol^{-1} \f$ for NAMD. */
#define sq6rt2  (1.12246204830937)      /*!< This constant is \f$ \sqrt[6]{2} \f$. */
#define PI      (3.14159265358979)      /*!< This is \f$ \pi \f$. */
#define watercomp (0.007372)		/*!< This the compressibility of water in internal units. */

#define MAXLIST 2048                    /*!< Maximal number of pairs per atom for neighbours list. */

#define X2(x) ((x)*(x))                 /*!< Raise any value at power 2. */
#define X3(x) (X2(x)*(x))               /*!< Raise any value at power 3. */
#define X4(x) (X2(x)*X2(x))             /*!< Raise any value at power 4. */
#define X6(x) (X3(x)*X3(x))             /*!< Raise any value at power 6. */
#define X12(x) (X6(x)*X6(x))            /*!< Raise any value at power 12. */

#define MAX(x,y) ((x)>=(y)?(x):(y))     /*!< Macro providing the maximum of (x,y). */
#define MIN(x,y) ((x)<=(y)?(x):(y))     /*!< Macro providing the minimum of (x,y). */

#define COSDIH  1                       /*!< Cosine potential type for dihedral and improper angles. */
#define HARMDIH 2                       /*!< Harmonic potential type for dihedral and improper angles. */

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
 * \brief Enumeration of the availables simulation boxes.
 */
enum BOX_TYPE{
  NOBOX = 0,        /*!< No box for this simulation, so no Periodic Boundaries Conditions applied. */
  CUBIC = 1,        /*< A cubic box is used. */
  ORBIC = 2,        /*< An orthorhombic box is used. */
  TCLIN = 3         /*< A triclinic box is used. */
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
  char label[5],segi[5],resi[5];
  int type,resn,inconst;

  double x,y,z;
  double vx,vy,vz;
  double fx,fy,fz;
  
  double m,q;
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
  double stress1,stress2,stress3,stress4,stress5,stress6,stress7,stress8,stress9;
}PBC;

/*!
 * \struct INPUTS
 * \brief This structure contains data parsed from input files.
 *
 * Most of the arrays defined in this structure are allocated when parsing CHARMM files and the freed.
 * 
 */
typedef struct
{
//   FILE *cntrFile,*topFile,*psfFile,*parmFile,*corFile;
  char **types;
  int *typesNum;
  int nTypes,nBondTypes,nAngTypes,nUbTypes,nDiheTypes,nImprTypes,nNonBonded;
  int **bondTypes,**angTypes,**ubTypes,*nDiheTypesParm,**diheTypes,**imprTypes;
  double **bondTypesParm,**angTypesParm,**ubTypesParm;
  double **diheTypesParm,**imprTypesParm,**nonBondedTypesParm;
}INPUTS;

/*!
 * \struct FORCEFIELD
 * \brief This structure contains most of the arrays used for storing parameters for the forcefield.
 */
typedef struct
{
  int nBond,nAngle,nDihedral,nImproper,nUb,ncells;
  int **verList,*verPair,**ver14,npr14,*nParmDihe;
  double **parmVdw,scal14;
  double **parmBond,**parmUb,**parmAngle,**parmDihe;
  double **parmImpr;
}FORCEFIELD;

/*!
 * \struct SIMULPARAMS
 * \brief This structure contains parameters of the simulation.
 */
typedef struct
{
  int natom,nb14,step,nsteps,degfree,firstener,listupdate;
  int printo,printtr,integrator,ens,nconst,maxcycle;
  int keyrand,seed,keytraj,keyener,keyforf,keymd,keyconsth;
  int keyminim,maxminst,mdNature,numDeriv;
  int linkRatio,keylink,nolink;
  int *excludeNum,**excludeAtom;
  int *bondType,*ubType,*angleType,*diheType,*imprType;
  int **iBond,**iUb,**iAngle,**iDihedral,**iImproper;
  double chargeConst,cutoff,cuton,delr,tolshake;
  double lambdat,lambdap,kintemp0,taut;
  double tolminim,maxminsiz,temp,timeStep;
  double press,tempStep,pressStep,taup;
  enum ELEC_TYPE elecType;
  enum VDW_TYPE  vdwType;
}SIMULPARAMS;

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

#endif
