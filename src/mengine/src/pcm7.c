#define EXTERN extern

#include "pcwin.h"
#include "pcmod.h"

#include "angles.h"
#include "attached.h"
#include "energies.h"
#include "bonds_ff.h"
#include "derivs.h"
#include "hess.h"
#include "utility.h"
#include "cutoffs.h"
#include "nonbond.h"
#include "torsions.h"
#include "fix.h"

EXTERN struct t_strbnd {
        int nstrbnd, **isb;
        float *ksb1, *ksb2;
        } strbnd;
EXTERN struct t_improp {
        int nimptors, **iiprop;
        float *v1, *v2, *v3;
        int   *ph1, *ph2, *ph3;        
        } improp; 
            
EXTERN struct t_minim_control {
        int type, method, field, added_const;
        char added_path[256],added_name[256];
        } minim_control;
        
EXTERN struct t_minim_values {
        int iprint, ndc, nconst;
        float dielc;
        } minim_values;
struct t_solvent {
    int type;
    double EPSin, EPSsolv;
    double doffset, p1,p2,p3,p4,p5;
    double *shct,*asolv,*rsolv,*vsolv,*gpol,*rborn;
    } solvent;

int Missing_constants;
EXTERN double hesscut;
    
void message_alert(char *,char *);
void pcm7_min(void);
int initialise(void);
void attach(void);
void get_bonds(void);
void get_angles(void);
void get_torsions(int natom,int *type,int **iat,int **bo);
void set_field(int);
double energy(void);
double get_total_energy(void);
double get_total_deriv_x(int i);
double get_total_deriv_y(int i);
double get_total_deriv_z(int i);
double print_energy(void);
void print_energy_totals(void);
void set_active(int natom,int *use,int natom_fix,int *katm_fix);
int setup_calculation(void);
void end_calculation(void);

void ebond(int nbnd,int **,int *use,double *x,double *y,double *z,double *bl,double *bk,double *estr);
void ebond1(int natom,int nbnd,int **,int *use,double *x,double *y,double *z,double *bl,double *bk,double *estr,double **deb);
void eangle(int nang,int **,int *use,double *x,double *y,double *z,int *angin,int *angtype,float *anat,float *acon,double *ebend);
void eangle1(int nang,int natom,int **,int *use,double *x,double *y,double *z,int *angin,int *angtype,float *anat,float *acon,double *ebend,double **dea);
void ebufcharge(int natom,int *use,int **skip,double *x,double *y,double *z,double *charge,double chrgcut,double *eu);
void ebufcharge1(int natom,int *use,int **skip,double *x,double *y,double *z,double *charge,double chrgcut,double *eu,double **deqq);
void ehal(int natom,int *type,int *use,double *x,double *y,double *z,double vdwcut,int **skip,double **vrad,double **veps,double *evdw,double *e14);
void ehal1(int natom,int *type,int *use,double *x,double *y,double *z,double vdwcut,int **skip,double **vrad,double **veps,double *evdw,double *e14,
           double **devdw,double **de14);
void etorsion(int ntor,int **i14,int *use,int *type,double *x,double *y,double *z,int *phin1,int *phin2,int *phin3,int *phin4,int *phin5,int *phin6,float *vin1,float *vin2,
              float *vin3,float *vin4,float *vin5,float *vin6,double *etor);
void etorsion1(int natom,int ntor,int **i14,int *use,double *x,double *y,double *z,int *phin1,int *phin2,int *phin3,int *phin4,int *phin5,int *phin6,float *vin1,float *vin2,
             float *vin3,float *vin4,float *vin5,float *vin6,double *etor,double **detor);
void eopbend_wilson(int nang,int **,int *use,double *x,double *y,double *z,int *angin,float *copb,double *eopb);
void eopbend_wilson1(int natom,int nang,int **,int *use,double *x,double *y,double *z,int *angin,float *copb,double *eopb,double **deopb);
void estrbnd(int nstrbnd,int **,int **,int *use,double *x,double *y,double *z,float *anat,double *bl,double *bk,float *ksb1,
             float *ksb2,double *estrbnd);
void estrbnd1(int natom,int nstrbnd,int **,int **,int *use,double *x,double *y,double *z,float *anat,double *bl,double *bk,float *ksb1,
	      float *ksb2,double *estrbnd,double **destb);
void egeom(int natom,int *use,double *x,double *y,double *z,double *egeom);
void egeom1(int natom,int *use,double *x,double *y,double *z,double *egeom,double **degeom);
void esolv(int natom,double chrgcut,int **skip,int *use,double *charge,double *x,double *y,double *z,double *esolv);
void esolv1(int natom,double chrgcut,int **skip,int *use,double *charge,double *x,double *y,double *z,double *esolv,double **desolv,double *drb);

void hessian(int, double *,int *, int *, int *,double *);
void ebond2(int ia,int **,int **,int **,double *x,double *y,double *z,double *bl,double *bk,float **hessx,float **hessy,float **hessz);
void eangle2(int i,int nang,int **,int *angin, double *x,double *y,double *z,float *acon,float *anat,double **dea,float **hessx,float **hessy,float **hessz);
void etorsion2(int iatom,int ntor,int **i14,int *use,double *x,double *y,double *z,int *phin1,int *phin2,int *phin3,int *phin4,int *phin5,int *phin6,float *vin1,float *vin2,        
             float *vin3,float *vin4,float *vin5,float *vin6,float **hessx,float **hessy,float **hessz);
void ehal2(int iatom,int natom,int *type,int *use,double *x,double *y,double *z,double vdwcut,int **skip,double **vrad,double **veps,
           float **hessx,float **hessy,float **hessz);
void ebufcharge2(int i,int natom,int *use,int **skip,double *x,double *y,double *z,double *charge,double chrgcut,float **hessx,float **hessy,float **hessz);
void eopbend_wilson2(int i,int nang,int **,int *use,double *x,double *y,double *z,int *angin,float *copb,
             float **hessx,float **hessy,float **hessz);
void estrbnd2(int iatom,int nstrbnd,int **,int **,int *use,double *x,double *y,double *z,float *anat,double *bl,double *bk,float *ksb1,
	      float *ksb2,float **hessx,float **hessy,float **hessz);
void esolv2(int ia,int natom,double chrgcut,double *charge,double *x,double *y,double *z,float **hessx,float **hessy,float **hessz);
void egeom2(int ,int natom,int *use,double *x,double *y,double *z,float **hessx,float **hessy,float **hessz);
// gaff and amber
void eimptors(int nimptors,int **iiprop,int *use,int *type,double *x,double *y,double *z,float *v1in,float *v2in,float *v3in,int *ph1in,int *ph2in,int *ph3in, double *e_imp);
void eimptors1(int nimptors,int natom,int **iiprop,int *use,int *type,double *x,double *y,double *z,float *v1in,float *v2in,float *v3in,int *ph1in,int *ph2in,int *ph3in,
           double *e_imp,double **deimprop);
void eimptors2(int iatom,int nimptors,int **iiprop,int *use,int *type,double *x,double *y,double *z,float *v1in,float *v2in,float *v3in,int *ph1in,int *ph2in,int *ph3in,
           float **hessx,float **hessy,float **hessz);
void elj(int natom,int *type,int *use,double *x,double *y,double *z,double vdwcut,int **skip,double **vrad,double **veps,double *evdw,double *e14);
void elj1(int natom,int *type,int *use,double *x,double *y,double *z,double vdwcut,int **skip,double **vrad,double **veps,double *evdw,double *e14,
           double **devdw,double **de14);
void elj2(int iatom,int natom,int *type,int *use,double *x,double *y,double *z,double vdwcut,int **skip,double **vrad,double **veps,
           float **hessx,float **hessy,float **hessz);
void echarge(int natom,int *use, int **skip,double *x,double *y,double *z,double *charge,double chrgcut,double *eu);
void echarge1(int natom,int *use, int **skip,double *x,double *y,double *z,double *charge,double chrgcut,double *eu,double **deqq);
void echarge2(int i,int natom,int *use, int **skip,double *x,double *y,double *z,double *charge,double chrgcut,float **hessx,float **hessy,float **hessz);
//

int  kbond(void);
void kangle(void);
void kstrbnd(void);
void kcharge(int natom,int *type,int *atomnum,long int *flags,int **iat,int **bo,double *charge,double *sigma_charge,double *formal_charge);
void kopbend(int *type,int *tclass);
void ktorsion(int *type,int *tclass,int **iat,int **bo);
void kvdw(int natom,int *type,int *atomnum);
void ksolv(int natom,int *atomnum,int **iat,int **bo);

// gaff
void kimptors(void);

void get_memory(void);
void free_memory(void);
void gradient(void);
void read_datafiles(char *);
void get_added_const(void);
void hdel(int);
void set_atomtypes(int);
void type(void);
void zero_data(void);
void potoff(void);
int use_bond(void);
int use_angle(void);
int use_strbnd(void);
int use_opbend_wilson(void);
int use_tors(void);
int use_strtor(void);
int use_hal(void);
int use_charge(void);
int use_bufcharge(void);
int use_geom(void);
int use_solv(void);
int use_imptor(void);
int use_lj(void);

// ==================================================================
int setup_calculation()
{
  int i,j,nRet;
    char string[30];
    double etot;

    nRet = 0;
    potoff(); // turn off all force field terms before resetting
    if (minim_control.field == MMFF94)
    {
        set_field(MMFF94);
        set_atomtypes(MMFF94);
        if (minim_control.added_const)
           get_added_const();
    } else if (minim_control.field == GAFF)
      {
        set_field(GAFF);
        set_atomtypes(GAFF);
      } else
    {
        if (minim_control.added_const)
        {
           strcpy(string,"mmxconst.prm");
           zero_data();
           read_datafiles(string);
           get_added_const();
        }
        type();
        set_field(MMX);
    }

            
  // allocate memeory for derivatives
  get_memory();
 
  // setup bonds, angles, torsions, improper torsions and angles, allenes
  get_bonds();
  get_angles();
  get_torsions(natom,atom.type,atom.iat,atom.bo);

  attach();
  // need allene

  // set active atoms
  set_active(natom,atom.use,fx_atom.natom_fix,fx_atom.katom_fix);
  // setup non_bonded list of atoms to skip
  for (i=1; i <= natom; i++)
    {
         for (j=0; j < MAXIAT; j++)
             if (atom.iat[i][j] != 0 && atom.bo[i][j] != 9)
               {
                    skip[i][atom.iat[i][j]] = i;
                    skip[atom.iat[i][j]][i] = atom.iat[i][j];
               }
         for (j=0; j < attached.n13[i]; j++)
           {
               skip[i][attached.i13[j][i]] = i;
               skip[attached.i13[j][i]][i] = attached.i13[j][i];
           }
         for(j=0; j < attached.n14[i]; j++)
           {
            skip[i][attached.i14[j][i]] = -i;
            skip[attached.i14[j][i]][i] = -attached.i14[j][i];
           }
    }
/*  assign local geometry potential function parameters  */
    Missing_constants = FALSE;
    if (use_bond() || use_strbnd())  nRet = kbond();
     if (nRet == FALSE) // use end_calc to free memory
     {
         energies.total = 9999.;
         return FALSE;
     }
     if (use_angle() || use_strbnd()) kangle();
     if (use_angle() || use_opbend_wilson()) kopbend(atom.type,atom.tclass);
     if (use_tors()  || use_strtor()) ktorsion(atom.type,atom.tclass,atom.iat,atom.bo);
     if (use_strbnd()) kstrbnd();
     if (use_imptor()) kimptors();
      
     if (use_hal()) kvdw(natom,atom.type,atom.atomnum);
     if (use_charge() || use_bufcharge()) kcharge(natom,atom.type,atom.atomnum,atom.flags,atom.iat,atom.bo,
          atom.charge,atom.sigma_charge,atom.formal_charge);
     if (use_solv()) ksolv(natom,atom.atomnum,atom.iat,atom.bo);

       if (Missing_constants == TRUE)
       {
           energies.total = 9999.0;
           return (FALSE);
       }    
       if (minim_values.iprint == TRUE)
          etot = print_energy();
       else
          etot = energy();
       return TRUE;
}
// =================================================================
void end_calculation()
{
 
    free_memory();
}
// ==================================================================
void gradient()
{
    int i, j;

      energies.estr = 0.0;
      energies.ebend = 0.0;
      energies.estrbnd = 0.0;
      energies.e14 = 0.0;
      energies.evdw = 0.0;
      energies.etor = 0.0;
      energies.eu = 0.0;
      energies.eopb = 0.0;
      energies.eangang = 0.0;
      energies.estrtor = 0.0;
      energies.ehbond = 0.0;
      energies.efix = 0.0;
      energies.eimprop = 0.0;
      energies.eimptors = 0.0;
      energies.eurey = 0.0;
      energies.esolv = 0.0;
      energies.egeom =0.0;
      
      for (i=1; i <= natom; i++)
      {
        for (j=0; j < 3; j++)
        {
            deriv.deb[i][j] = 0.0;
            deriv.dea[i][j] = 0.0;
            deriv.destb[i][j] = 0.0;
            deriv.deopb[i][j] = 0.0;
            deriv.detor[i][j] = 0.0;
            deriv.de14[i][j] = 0.0;
            deriv.devdw[i][j]= 0.0;
            deriv.deqq[i][j] = 0.0;
            deriv.deaa[i][j] = 0.0;
            deriv.destor[i][j] = 0.0;
            deriv.deimprop[i][j] = 0.0;
            deriv.dehb[i][j] = 0.0;
            deriv.desolv[i][j] = 0.0;
            deriv.degeom[i][j] = 0.0;
        }
      }
   
      if (use_bond())ebond1(natom,bonds_ff.nbnd,bonds_ff.i12,atom.use,atom.x,atom.y,atom.z,bonds_ff.bl,bonds_ff.bk,&energies.estr,deriv.deb);
      if (use_angle())eangle1(angles.nang,natom, angles.i13,atom.use,atom.x,atom.y,atom.z,angles.angin,angles.angtype,angles.anat,angles.acon,&energies.ebend,deriv.dea);
      if (use_tors())etorsion1(natom,torsions.ntor,torsions.i14,atom.use,atom.x,atom.y,atom.z,torsions.ph1,torsions.ph2,torsions.ph3,torsions.ph4,torsions.ph5,torsions.ph6,
				 torsions.v1,torsions.v2,torsions.v3,torsions.v4,torsions.v5,torsions.v6,&energies.etor,deriv.detor);
      // mmff 
      if (use_strbnd())estrbnd1(natom,strbnd.nstrbnd,strbnd.isb,angles.i13,atom.use,atom.x,atom.y,atom.z,angles.anat,bonds_ff.bl,bonds_ff.bk,strbnd.ksb1,strbnd.ksb2,
				  &energies.estrbnd,deriv.destb);
      if (use_opbend_wilson()) eopbend_wilson1(natom,angles.nang,angles.i13,atom.use,atom.x,atom.y,atom.z,angles.angin,angles.copb,&energies.eopb,deriv.deopb);
      if (use_hal()) ehal1(natom,atom.type,atom.use,atom.x,atom.y,atom.z,cutoffs.vdwcut,skip,nonbond.vrad,nonbond.veps,&energies.evdw,&energies.e14,deriv.devdw,deriv.de14);
      if (use_bufcharge()) ebufcharge1(natom,atom.use,skip,atom.x,atom.y,atom.z,atom.charge,cutoffs.chrgcut,&energies.eu,deriv.deqq);
      // gaff
      if (use_imptor()) eimptors1(improp.nimptors,natom,improp.iiprop,atom.use,atom.type,atom.x,atom.y,atom.z,improp.v1,improp.v2,improp.v3,improp.ph1,improp.ph2,improp.ph3,
              &energies.eimptors,deriv.deimprop);
      if (use_lj()) elj1(natom,atom.type,atom.use,atom.x,atom.y,atom.z,cutoffs.vdwcut,skip,nonbond.vrad,nonbond.veps,&energies.evdw,&energies.e14,deriv.devdw,deriv.de14);
      if (use_charge()) echarge1(natom,atom.use,skip,atom.x,atom.y,atom.z,atom.charge,cutoffs.chrgcut,&energies.eu,deriv.deqq);

      if (use_geom()) egeom1(natom,atom.use,atom.x,atom.y,atom.z,&energies.egeom,deriv.degeom);
      if (use_solv()) esolv1(natom,cutoffs.chrgcut,skip,atom.use,atom.charge,atom.x,atom.y,atom.z,&energies.esolv,deriv.desolv,deriv.drb);
     
      energies.total =  energies.estr + energies.ebend + energies.etor + energies.estrbnd + energies.e14+
                        energies.evdw + energies.eu + energies.ehbond + energies.eangang + energies.estrtor +
                        energies.eimprop + energies.eimptors + energies.eopb + energies.eurey + energies.esolv + energies.egeom;
      for (i=1; i <= natom; i++)
      {
          for (j=0; j < 3; j++)
          {
              deriv.d1[i][j] = deriv.deb[i][j] + deriv.dea[i][j] + deriv.deaa[i][j] +
                               deriv.destb[i][j] + deriv.detor[i][j] + deriv.deopb[i][j] + deriv.dehb[i][j] +
                               deriv.destor[i][j] + deriv.deqq[i][j] + deriv.devdw[i][j] + deriv.de14[i][j] +
                               deriv.deimprop[i][j] + deriv.deub[i][j] + deriv.desolv[i][j] + deriv.degeom[i][j];

         }
      }
}
// ==================================================================
double get_total_deriv_x(int i)
{
  return (deriv.d1[i][0]);
}
double get_total_deriv_y(int i)
{
  return (deriv.d1[i][1]);
}
double get_total_deriv_z(int i)
{
  return (deriv.d1[i][2]);
}
// ============================================
void print_energy_totals()
{
     fprintf(pcmlogfile,"\nEnergy : %f\n",energies.total);
     fprintf(pcmlogfile,"Str: %f    Bnd: %f    Tor: %f\n", energies.estr, energies.ebend, energies.etor);
     fprintf(pcmlogfile,"StrBnd: %f VDW: %f     QQ: %f\n",energies.estrbnd,energies.evdw+energies.e14+energies.ehbond, energies.eu);
     fprintf(pcmlogfile,"OOP: %f AA: %f     Strtor: %f\n",energies.eopb,energies.eangang, energies.estrtor);
     if (use_solv()) fprintf(pcmlogfile,"Solv: %f \n",energies.esolv);
}
// ==================================================================
double get_total_energy()
{
  return energies.total;
}
// ==================================================================
double energy()
{
    double etot;
    
      energies.total = 0.0;
      energies.estr = 0.0;
      energies.ebend = 0.0;
      energies.etor = 0.0;
      energies.estrbnd = 0.0;
      energies.evdw = 0.0;
      energies.e14 = 0.0;
      energies.ehbond = 0.0;
      energies.eu = 0.0;
      energies.eangang = 0.0;
      energies.estrtor = 0.0;
      energies.eimprop = 0.0;
      energies.eimptors = 0.0;
      energies.eopb = 0.0;
      energies.eurey = 0.0;
      energies.esolv = 0.0;
      energies.egeom = 0.0;


      if (use_bond())ebond(bonds_ff.nbnd,bonds_ff.i12,atom.use,atom.x,atom.y,atom.z,bonds_ff.bl,bonds_ff.bk,&energies.estr);
      if (use_angle())eangle(angles.nang,angles.i13,atom.use,atom.x,atom.y,atom.z,angles.angin,angles.angtype,angles.anat,angles.acon,&energies.ebend);
      if (use_tors())etorsion(torsions.ntor,torsions.i14,atom.use,atom.type,atom.x,atom.y,atom.z,torsions.ph1,torsions.ph2,torsions.ph3,torsions.ph4,torsions.ph5,
           torsions.ph6,torsions.v1,torsions.v2,torsions.v3,torsions.v4,torsions.v5,torsions.v6,&energies.etor);
      // mmff 
      if (use_strbnd())estrbnd(strbnd.nstrbnd,strbnd.isb,angles.i13,atom.use,atom.x,atom.y,atom.z,angles.anat,bonds_ff.bl,bonds_ff.bk,strbnd.ksb1,strbnd.ksb2,
				  &energies.estrbnd);
      if (use_opbend_wilson()) eopbend_wilson(angles.nang,angles.i13,atom.use,atom.x,atom.y,atom.z,angles.angin,angles.copb,&energies.eopb);
      if (use_hal()) ehal(natom,atom.type,atom.use,atom.x,atom.y,atom.z,cutoffs.vdwcut,skip,nonbond.vrad,nonbond.veps,&energies.evdw,&energies.e14);
      if (use_bufcharge()) ebufcharge(natom,atom.use,skip,atom.x,atom.y,atom.z,atom.charge,cutoffs.chrgcut,&energies.eu);

      // gaff
      if (use_imptor()) eimptors(improp.nimptors,improp.iiprop,atom.use,atom.type,atom.x,atom.y,atom.z,improp.v1,improp.v2,improp.v3,improp.ph1,improp.ph2,improp.ph3,&energies.eimptors);
      if (use_lj()) elj(natom,atom.type,atom.use,atom.x,atom.y,atom.z,cutoffs.vdwcut,skip,nonbond.vrad,nonbond.veps,&energies.evdw,&energies.e14);
      if (use_charge()) echarge(natom,atom.use,skip,atom.x,atom.y,atom.z,atom.charge,cutoffs.chrgcut,&energies.eu);

      if (use_geom()) egeom(natom,atom.use,atom.x,atom.y,atom.z,&energies.egeom);
      if (use_solv()) esolv(natom,cutoffs.chrgcut,skip,atom.use,atom.charge,atom.x,atom.y,atom.z,&energies.esolv);
                              
      energies.total =  energies.estr + energies.ebend + energies.etor + energies.estrbnd
             + energies.evdw + energies.e14 + energies.ehbond + energies.eu + energies.eangang +
               energies.estrtor + energies.eimprop + energies.eimptors + energies.eopb + energies.eurey + energies.esolv + energies.egeom;
      etot = energies.total;
      return (etot);
}
/* =========================================================================== */ 
double print_energy()
{
    double etot;
    
      energies.total = 0.0;
      energies.estr = 0.0;
      energies.ebend = 0.0;
      energies.etor = 0.0;
      energies.estrbnd = 0.0;
      energies.evdw = 0.0;
      energies.e14 = 0.0;
      energies.ehbond = 0.0;
      energies.eu = 0.0;
      energies.eangang = 0.0;
      energies.estrtor = 0.0;
      energies.eimprop = 0.0;
      energies.eimptors = 0.0;
      energies.eopb = 0.0;
      energies.eurey = 0.0;
      energies.esolv = 0.0;
      energies.egeom = 0.0;


      if (use_bond())ebond(bonds_ff.nbnd,bonds_ff.i12,atom.use,atom.x,atom.y,atom.z,bonds_ff.bl,bonds_ff.bk,&energies.estr);
      if (use_angle())eangle(angles.nang,angles.i13,atom.use,atom.x,atom.y,atom.z,angles.angin,angles.angtype,angles.anat,angles.acon,&energies.ebend);
      if (use_tors())etorsion(torsions.ntor,torsions.i14,atom.use,atom.type,atom.x,atom.y,atom.z,torsions.ph1,torsions.ph2,torsions.ph3,torsions.ph4,torsions.ph5,
           torsions.ph6,torsions.v1,torsions.v2,torsions.v3,torsions.v4,torsions.v5,torsions.v6,&energies.etor);
      // mmff 
      if (use_opbend_wilson()) eopbend_wilson(angles.nang,angles.i13,atom.use,atom.x,atom.y,atom.z,angles.angin,angles.copb,&energies.eopb);
      if (use_strbnd())estrbnd(strbnd.nstrbnd,strbnd.isb,angles.i13,atom.use,atom.x,atom.y,atom.z,angles.anat,bonds_ff.bl,bonds_ff.bk,strbnd.ksb1,strbnd.ksb2,
				  &energies.estrbnd);
      if (use_hal()) ehal(natom,atom.type,atom.use,atom.x,atom.y,atom.z,cutoffs.vdwcut,skip,nonbond.vrad,nonbond.veps,&energies.evdw,&energies.e14);
      if (use_bufcharge()) ebufcharge(natom,atom.use,skip,atom.x,atom.y,atom.z,atom.charge,cutoffs.chrgcut,&energies.eu);

      // gaff
      if (use_imptor()) eimptors(improp.nimptors,improp.iiprop,atom.use,atom.type,atom.x,atom.y,atom.z,improp.v1,improp.v2,improp.v3,improp.ph1,improp.ph2,improp.ph3,&energies.eimptors);
      if (use_lj()) elj(natom,atom.type,atom.use,atom.x,atom.y,atom.z,cutoffs.vdwcut,skip,nonbond.vrad,nonbond.veps,&energies.evdw,&energies.e14);
      if (use_charge()) echarge(natom,atom.use,skip,atom.x,atom.y,atom.z,atom.charge,cutoffs.chrgcut,&energies.eu);
      
      if (use_geom()) egeom(natom,atom.use,atom.x,atom.y,atom.z,&energies.egeom);
      if (use_solv()) esolv(natom,cutoffs.chrgcut,skip,atom.use,atom.charge,atom.x,atom.y,atom.z,&energies.esolv);

      energies.total =  energies.estr + energies.ebend + energies.etor + energies.estrbnd
             + energies.evdw + energies.e14 + energies.ehbond + energies.eu + energies.eangang +
               energies.estrtor + energies.eimprop + energies.eimptors + energies.eopb + energies.eurey + energies.esolv + energies.egeom;
      etot = energies.total;
      
      return (etot);
}    
// ====================== hessian =======================
void hessian(int maxhess, double *h, int *hinit, int *hstop,int *hindex, double *hdiag)
{
   int i,j,k,nhess;
   double hmax;
   char hatext[60];
   
   int *keep;  //  keep[MAXATOM];

   keep = ivector(0,natom);   

   nhess = 0;
   for (i=1; i <= natom; i++)
   {

       hinit[nhess] = 0;
       hstop[nhess] = 0;
       hdiag[nhess] = 0.0;
       nhess++;
       hinit[nhess] = 0;
       hstop[nhess] = 0;
       hdiag[nhess] = 0.0;
       nhess++;
       hinit[nhess] = 0;
       hstop[nhess] = 0;
       hdiag[nhess] = 0.0;
       nhess++;       
   }
/*   periodic boundary conditions set here  */

 //   if (use_picalc) piseq(0,0);

    nhess= 0;
    for (i=1; i <= natom; i++)
    {
       for(k=1; k <= natom; k++)
       {
           for(j=0; j < 3; j++)
           {
              hess.hessx[k][j] = 0.0;
              hess.hessy[k][j] = 0.0;
              hess.hessz[k][j] = 0.0;
           }
       }
          
      if (atom.use[i])
      {
        if (use_bond())ebond2(i,bonds_ff.i12,atom.iat,atom.bo,atom.x,atom.y,atom.z,bonds_ff.bl,bonds_ff.bk,hess.hessx,hess.hessy,hess.hessz);
        if (use_angle())eangle2(i,angles.nang,angles.i13,angles.angin,atom.x,atom.y,atom.z,angles.acon,angles.anat,deriv.dea,hess.hessx,hess.hessy,hess.hessz); 
        if (use_tors())etorsion2(i,torsions.ntor,torsions.i14,atom.use,atom.x,atom.y,atom.z,torsions.ph1,torsions.ph2,torsions.ph3,torsions.ph4,torsions.ph5,
				   torsions.ph6,torsions.v1,torsions.v2,torsions.v3,torsions.v4,torsions.v5,torsions.v6,hess.hessx,hess.hessy,hess.hessz);
        
	// mmff
        if (use_strbnd())estrbnd2(i,strbnd.nstrbnd,strbnd.isb,angles.i13,atom.use,atom.x,atom.y,atom.z,angles.anat,bonds_ff.bl,bonds_ff.bk,strbnd.ksb1,strbnd.ksb2,
				    hess.hessx,hess.hessy,hess.hessz);
        if (use_opbend_wilson()) eopbend_wilson2(i,angles.nang,angles.i13,atom.use,atom.x,atom.y,atom.z,angles.angin,angles.copb,hess.hessx,hess.hessy,hess.hessz);
        if (use_hal())ehal2(i,natom,atom.type,atom.use,atom.x,atom.y,atom.z,cutoffs.vdwcut,skip,nonbond.vrad,nonbond.veps,hess.hessx,hess.hessy,hess.hessz);  
        if (use_bufcharge())ebufcharge2(i,natom,atom.use,skip,atom.x,atom.y,atom.z,atom.charge,cutoffs.chrgcut,hess.hessx,hess.hessy,hess.hessz);

	// gaff
        if (use_imptor()) eimptors2(i,improp.nimptors,improp.iiprop,atom.use,atom.type,atom.x,atom.y,atom.z,improp.v1,improp.v2,improp.v3,improp.ph1,improp.ph2,improp.ph3,hess.hessx,
				   hess.hessy,hess.hessz);
        if (use_lj()) elj2(i,natom,atom.type,atom.use,atom.x,atom.y,atom.z,cutoffs.vdwcut,skip,nonbond.vrad,nonbond.veps,hess.hessx,hess.hessy,hess.hessz);
	if (use_charge()) echarge2(i,natom,atom.use,skip,atom.x,atom.y,atom.z,atom.charge,cutoffs.chrgcut,hess.hessx,hess.hessy,hess.hessz);

        if (use_geom()) egeom2(i,natom,atom.use,atom.x,atom.y,atom.z,hess.hessx,hess.hessy,hess.hessz);
        if (use_solv()) esolv2(i,natom,cutoffs.chrgcut,atom.charge,atom.x,atom.y,atom.z,hess.hessx,hess.hessy,hess.hessz);

        hdiag[(i-1)*3]   += hess.hessx[i][0];
        hdiag[(i-1)*3+1] += hess.hessy[i][1];
        hdiag[(i-1)*3+2] += hess.hessz[i][2];

        for (k=i+1; k <= natom; k++)
        {
          keep[k] = FALSE;
          if (atom.use[k])
          {
              hmax = fabs(hess.hessx[k][0]);
              for (j=0; j < 3;j++)
              {
                  if (fabs(hess.hessx[k][j]) > hmax)
                     hmax = fabs(hess.hessx[k][j]);
              }
              for (j=0; j < 3;j++)
              {
                  if (fabs(hess.hessy[k][j]) > hmax)
                     hmax = fabs(hess.hessy[k][j]);
              }
              for (j=0; j < 3;j++)
              {
                  if (fabs(hess.hessz[k][j]) > hmax)
                     hmax = fabs(hess.hessz[k][j]);
              }
              if (hmax >= hesscut) keep[k] = TRUE;
          }
        }
      
        hinit[(i-1)*3] = nhess;  //  x component of new atom
      
        hindex[nhess] = 3*(i-1) + 1;      
        h[nhess] = hess.hessx[i][1];
        nhess++;
        hindex[nhess] = 3*(i-1) + 2;       
        h[nhess] = hess.hessx[i][2];
        nhess++;
      
        for (k=i+1; k <= natom; k++)
        {
          if (keep[k])
          {
              for (j=0; j < 3; j++)
              {
                  hindex[nhess] = 3*(k-1) +j;
                  h[nhess] = hess.hessx[k][j];
                  nhess++;
              }
          }
        }
        hstop[(i-1)*3] = nhess;
      
        hinit[(i-1)*3+1] = nhess;  // y component of atom i
        hindex[nhess] = 3*(i-1) + 2;
        h[nhess] = hess.hessy[i][2];
        nhess++;
        for (k=i+1; k <= natom; k++)
        {
          if (keep[k])
          {
              for (j=0; j < 3; j++)
              {
                  hindex[nhess] = 3*(k-1) +j;
                  h[nhess] = hess.hessy[k][j];
                  nhess++;
              }
          }
        }
        hstop[(i-1)*3+1] = nhess;
        hinit[(i-1)*3+2] = nhess;
        for (k=i+1; k <= natom; k++)
        {
          if (keep[k])
          {
              for (j=0; j < 3; j++)
              {
                  hindex[nhess] = 3*(k-1) +j;
                  h[nhess] = hess.hessz[k][j];
                  nhess++;
              }
          }
        }
        hstop[(i-1)*3+2] = nhess;

        if (nhess > maxhess)
        {
          sprintf(hatext,"Too many hessian elements: %d elements  %d max\n",nhess,maxhess);
          message_alert(hatext,"Error");
          fprintf(pcmlogfile,"Too many hessian elements !\n");
          strcpy(Openbox.fname,"Hesserr.pcm");
          exit(0);
        }
      }
    }
    free_ivector(keep, 0,natom);
}

