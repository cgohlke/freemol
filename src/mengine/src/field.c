#define EXTERN extern

#include "pcwin.h"
#include "pot.h"
#include "field.h"
#include "fix.h"
#include "job_control.h"

EXTERN struct t_minim_values {
        int iprint, ndc, nconst;
        float dielc;
        } minim_values;

void set_field(int type);
void potoff(void);
int get_field(void);
int use_solvation(void);
void set_solvation(int);
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

void set_field_name(char *);
void set_field_bondunit(float ftemp);
void set_field_bondcubic(float ftemp);
void set_field_bondquartic(float ftemp);
void set_field_angleunit(float ftemp);
void set_field_anglecubic(float ftemp);
void set_field_anglequartic(float ftemp);
void set_field_anglepentic(float ftemp);
void set_field_anglesextic(float ftemp);
void set_field_strbndunit(float ftemp);
void set_field_angangunit(float ftemp);
void set_field_strtorunit(float ftemp);
void set_field_torsionunit(float ftemp);
void set_field_vdwtype(char *name);
void set_field_radiustype(char *name);
void set_field_radiussize(char *name);
void set_field_radiusrule(char *name);
void set_field_epsrule(char *name);
void set_field_aterm(float ftemp);
void set_field_bterm(float ftemp);
void set_field_cterm(float ftemp);
void set_field_vdwscale(float ftemp);
void set_field_chrgscale(float ftemp);
void set_field_dielectric(float ftemp);

char * get_radiustype(void);
char * get_radiussize(void);
char * get_radiusrule(void);
char * get_epsrule(void);
// ======================
int get_field()
{
  return field.type;
}
// ====================
char * get_radiustype()
{
  return field.radiustype;
}
char * get_radiussize()
{
  return field.radiussize;
}
char * get_radiusrule()
{
  return field.radiusrule;
}
char * get_epsrule()
{
  return field.radiusrule;
}
// =================
void set_field_name(char *name)
{
  strcpy(field.name,name);
}
// =========== bonds =============
void set_field_bondunit(float ftemp)
{
  field.bondunit = ftemp;
}
void set_field_bondcubic(float ftemp)
{
  field.bond_cubic = ftemp;
}
void set_field_bondquartic(float ftemp)
{
  field.bond_quartic = ftemp;
}
// ============ angles ==============
void set_field_angleunit(float ftemp)
{
  field.angleunit = ftemp;
}
void set_field_anglecubic(float ftemp)
{
  field.angle_cubic = ftemp;
}
void set_field_anglequartic(float ftemp)
{
  field.angle_quartic = ftemp;
}
void set_field_anglepentic(float ftemp)
{
  field.angle_pentic = ftemp;
}
void set_field_anglesextic(float ftemp)
{
  field.angle_sextic = ftemp;
}
// =========== cross terms ===============
void set_field_strbndunit(float ftemp)
{
  field.str_bndunit = ftemp;
}
void set_field_angangunit(float ftemp)
{
  field.ang_angunit = ftemp;
}
void set_field_strtorunit(float ftemp)
{
  field.str_torunit = ftemp;
}
// =============== torsions ===========
void set_field_torsionunit(float ftemp)
{
  field.torsionunit = ftemp;
}
// ================ vdw ==============
void set_field_vdwtype(char *name)
{
  strcpy(field.vdwtype,name);
}
void set_field_radiustype(char *name)
{
  strcpy(field.radiustype,name);
}
void set_field_radiussize(char *name)
{
  strcpy(field.radiussize,name);
}
void set_field_radiusrule(char *name)
{
  strcpy(field.radiusrule,name);
}
void set_field_epsrule(char *name)
{
  strcpy(field.epsrule,name);
}
void set_field_aterm(float ftemp)
{
  field.a_expterm = ftemp;
}
void set_field_bterm(float ftemp)
{
  field.b_expterm = ftemp;
}
void set_field_cterm(float ftemp)
{
  field.c_expterm = ftemp;
}
void set_field_vdwscale(float ftemp)
{
  field.vdw_14scale = ftemp;
}
void set_field_dielectric(float ftemp)
{
  field.dielectric = ftemp;
}
void set_field_chrgscale(float ftemp)
{
  field.chg_14scale = ftemp;
}
// ================
void set_field(int type)
{
  field.type = type;
      pot.use_bond = FALSE;
      pot.use_angle = FALSE;
      pot.use_strbnd = FALSE;
      pot.use_urey = FALSE;
      pot.use_angang = FALSE;
      pot.use_opbend = FALSE;
      pot.use_improp = FALSE;
      pot.use_imptor = FALSE;
      pot.use_tors = FALSE;
      pot.use_strtor = FALSE;
      pot.use_tortor = FALSE;
      pot.use_vdw = FALSE;
      pot.use_lj = FALSE;
      pot.use_buck = FALSE;
      pot.use_hal = FALSE;
      pot.use_gauss = FALSE;
      pot.use_charge = FALSE;
      pot.use_bufcharge = FALSE;
      pot.use_chrgdpl = FALSE;
      pot.use_dipole = FALSE;
      pot.use_polar = FALSE;
      pot.use_geom = FALSE;
      pot.use_extra = FALSE;
      pot.use_picalc = FALSE;
      pot.use_hbond = FALSE;
      pot.use_coordb = FALSE;
      pot.use_opbend_wilson = FALSE;
      pot.use_highcoord = FALSE;

    if (type == MMX )
    {
        pot.use_bond = TRUE;
        pot.use_angle = TRUE;
        pot.use_strbnd = TRUE;
        pot.use_opbend = TRUE;
        pot.use_opbend_wilson = FALSE;
        pot.use_tors = TRUE;
        pot.use_vdw = TRUE;
        pot.use_buck = TRUE;
        pot.use_lj = FALSE;
        pot.use_hal = FALSE;
        pot.use_gauss = FALSE;
        pot.use_hbond = FALSE;
        
        pot.use_bufcharge= FALSE;
        if (minim_values.ndc == 4)
        {
             pot.use_charge = TRUE;
             pot.use_dipole = FALSE;
       }else
        {
             pot.use_charge = FALSE;
             pot.use_dipole = TRUE;
       }
        pot.use_urey = FALSE;
        
        pot.use_angang = FALSE;
        pot.use_improp = FALSE;
        pot.use_imptor = FALSE;
        pot.use_strtor = FALSE;
        pot.use_tortor = FALSE;
    } else if (type == GAFF)
    {
        pot.use_bond = TRUE;
        pot.use_angle = TRUE;
        pot.use_strbnd = FALSE;
        pot.use_opbend = FALSE;
        pot.use_opbend_wilson = FALSE;
        pot.use_tors = TRUE;
        pot.use_vdw = FALSE;
        pot.use_buck = FALSE;
        pot.use_lj = TRUE;
        pot.use_hal = FALSE;
        pot.use_hbond = FALSE;
        
        pot.use_gauss = FALSE;
        pot.use_charge = TRUE;
        pot.use_bufcharge= FALSE;
        pot.use_dipole = FALSE;
        pot.use_urey = FALSE;
        pot.use_improp = FALSE;
        pot.use_imptor = TRUE;
        pot.use_angang = FALSE;
        pot.use_strtor = FALSE;
    } else if (type == MMFF94)
    {
        pot.use_bond = TRUE;
        pot.use_angle = TRUE;
        pot.use_strbnd = TRUE;
        pot.use_opbend = FALSE;
        pot.use_opbend_wilson = TRUE;
        pot.use_tors = TRUE;
        
        pot.use_vdw = FALSE;
        pot.use_buck = FALSE;
        pot.use_lj = FALSE;
        pot.use_hal = TRUE;
        
        pot.use_hbond = FALSE;
        pot.use_picalc = FALSE;
        
        pot.use_gauss = FALSE;
        pot.use_charge = FALSE;
        pot.use_bufcharge= TRUE;
        pot.use_dipole = FALSE;
        pot.use_urey = FALSE;
        pot.use_improp = FALSE;
        pot.use_imptor = FALSE;
        pot.use_angang = FALSE;
        pot.use_strtor = FALSE;
    }
    if (fx_dist.ndfix > 0) pot.use_geom = TRUE;
    if (fx_angle.nafix > 0) pot.use_geom = TRUE;
    if (fx_torsion.ntfix > 0) pot.use_geom = TRUE;
    if (restrain_atom.natom_restrain > 0) pot.use_geom = TRUE;
    if (job_control.use_gbsa) pot.use_solv = TRUE;
}
// ============================
void potoff()
{
      pot.use_bond = FALSE;
      pot.use_angle = FALSE;
      pot.use_strbnd = FALSE;
      pot.use_urey = FALSE;
      pot.use_angang = FALSE;
      pot.use_opbend = FALSE;
      pot.use_improp = FALSE;
      pot.use_imptor = FALSE;
      pot.use_tors = FALSE;
      pot.use_strtor = FALSE;
      pot.use_tortor = FALSE;
      pot.use_vdw = FALSE;
      pot.use_lj = FALSE;
      pot.use_buck = FALSE;
      pot.use_hal = FALSE;
      pot.use_gauss = FALSE;
      pot.use_charge = FALSE;
      pot.use_bufcharge = FALSE;
      pot.use_chrgdpl = FALSE;
      pot.use_dipole = FALSE;
      pot.use_polar = FALSE;
      pot.use_solv = FALSE;
      pot.use_geom = FALSE;
      pot.use_extra = FALSE;
      pot.use_picalc = FALSE;
      pot.use_hbond = FALSE;
      pot.use_coordb = FALSE;
      pot.use_opbend_wilson = FALSE;
}
// ================ Solvation information  ===========
static int SOLVATION = FALSE;

int use_solvation()
{
  if (SOLVATION)
    return TRUE;
  else
    return FALSE;
}
void set_solvation(int mode)
{
  if (mode)
    SOLVATION = TRUE;
  else
    SOLVATION = FALSE;
}
// ==================
int use_bond(void)
{
  if (pot.use_bond)
    return TRUE;
  else 
    return FALSE;
}
int use_angle(void)
{
  if (pot.use_angle)
    return TRUE;
  else 
    return FALSE;
}
int use_strbnd(void)
{
  if (pot.use_strbnd)
    return TRUE;
  else 
    return FALSE;
}
int use_opbend_wilson(void)
{
  if (pot.use_opbend_wilson)
    return TRUE;
  else 
    return FALSE;
}
int use_tors(void)
{
  if (pot.use_tors)
    return TRUE;
  else 
    return FALSE;
}
int use_strtor(void)
{
  if (pot.use_strtor)
    return TRUE;
  else 
    return FALSE;
}
int use_imptor(void)
{
  if (pot.use_imptor)
    return TRUE;
  else 
    return FALSE;
}
int use_hal(void)
{
  if (pot.use_hal)
    return TRUE;
  else 
    return FALSE;
}
int use_lj(void)
{
  if (pot.use_lj)
    return TRUE;
  else 
    return FALSE;
}
int use_charge(void)
{
  if (pot.use_charge)
    return TRUE;
  else 
    return FALSE;
}
int use_bufcharge(void)
{
  if (pot.use_bufcharge)
    return TRUE;
  else 
    return FALSE;
}
int use_geom(void)
{
  if (pot.use_geom)
    return TRUE;
  else 
    return FALSE;
}
int use_solv(void)
{
  if (pot.use_solv)
    return TRUE;
  else 
    return FALSE;
}
