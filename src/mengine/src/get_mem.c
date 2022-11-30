#define EXTERN extern

#include "pcwin.h"
#include "pcmod.h"
#include "angles.h"
#include "attached.h"
#include "bonds_ff.h"
#include "derivs.h"
#include "hess.h"
#include "utility.h"
#include "nonbond.h"
#include "torsions.h"

void get_molecule_memory(int);
int use_solv(void);
void allocate_rings(int);
void max_torsions(int natom,int *type,int **iat,int **bo);
void get_memory(void);
void free_memory(void);
void free_rings(int);
void free_molecule_memory(void);

#define STILL 1
#define HCT   2

EXTERN struct t_solvent {
    int type;
    double EPSin, EPSsolv;
    double doffset, p1,p2,p3,p4,p5;
    double *shct,*asolv,*rsolv,*vsolv,*gpol,*rborn;
    } solvent;

EXTERN struct t_strbnd {
        int nstrbnd, **isb;
        float *ksb1, *ksb2;
        } strbnd;

// we can not tell if the input structure has hydrogens - coming from pymol it should
// but coming from somewhere else it might not. So allocate enough memory to add hydrogens
// =======================================
void get_molecule_memory(int niatom)
{
  int i,j;

  if (niatom == 0)
    return;
  MAXATOM = niatom + 2*niatom +10;
  MAXBND = 6*MAXATOM/4;
  MAXANG = 12*MAXATOM/4;
  MAXTOR = 6*MAXATOM;

  atom.type = ivector(0,MAXATOM);
  atom.tclass = ivector(0,MAXATOM);
  atom.mmx_type = ivector(0,MAXATOM);
  atom.gaff_type = ivector(0,MAXATOM);
  atom.mmff_type = ivector(0,MAXATOM);
  atom.atomnum = ivector(0,MAXATOM);
  atom.use = ivector(0,MAXATOM);
  atom.input_charge = ivector(0,MAXATOM);
  atom.flags = ilvector(0,MAXATOM);
  atom.nconnect = ivector(0,MAXATOM);
  atom.tbo = ivector(0,MAXATOM);
  atom.iat = imatrix(0, MAXATOM,0,MAXIAT);
  atom.bo = imatrix(0,MAXATOM, 0, MAXIAT);
  atom.x = dvector(0,MAXATOM);
  atom.y = dvector(0,MAXATOM);
  atom.z = dvector(0,MAXATOM);
  atom.charge = dvector(0,MAXATOM);
  atom.formal_charge = dvector(0,MAXATOM);
  atom.sigma_charge = dvector(0,MAXATOM);
  atom.atomwt = dvector(0,MAXATOM);
  atom.radius = dvector(0,MAXATOM);
  atom.name = malloc( (MAXATOM)*sizeof(LABEL));
 
  for (i=0; i < MAXATOM; i++)
    {
      atom.type[i] = 0;
      atom.tclass[i] = 0;
      atom.mmx_type[i] = 0;
      atom.gaff_type[i] = 0;
      atom.mmff_type[i] = 0;
      atom.atomnum[i] = 0;
      atom.use[i] = 0;
      atom.input_charge[i] = 0;
      atom.nconnect[i] = 0;
      atom.tbo[i] = 0;
      atom.flags[i] = 0;
      atom.x[i] = 0.0;
      atom.y[i] = 0.0;
      atom.z[i] = 0.0;
      atom.charge[i] = 0.0;
      atom.formal_charge[i] = 0.0;
      atom.sigma_charge[i] = 0.0;
      atom.atomwt[i] = 0.0;
      for (j=0; j < MAXIAT; j++)
	{
	  atom.iat[i][j] = 0;
	  atom.bo[i][j] = 0;
	}
      strcpy(atom.name[i],"");
    }
  allocate_rings(niatom); // allocate space for ring data
}
// ====================================
void free_molecule_memory()
{
  if (natom == 0)
    return;
  free_rings(natom);
  free_ivector(atom.type,0,MAXATOM);
  free_ivector(atom.tclass ,0,MAXATOM);
  free_ivector(atom.mmx_type ,0,MAXATOM);
  free_ivector(atom.gaff_type ,0,MAXATOM);
  free_ivector(atom.mmff_type ,0,MAXATOM);
  free_ivector(atom.atomnum ,0,MAXATOM);
  free_ivector(atom.use ,0,MAXATOM);
  free_ivector(atom.input_charge ,0,MAXATOM);
  free_ilvector(atom.flags ,0,MAXATOM);
  free_ivector(atom.nconnect ,0,MAXATOM);
  free_ivector(atom.tbo ,0,MAXATOM);
  free_imatrix(atom.iat ,0, MAXATOM,0,MAXIAT);
  free_imatrix(atom.bo ,0,MAXATOM, 0, MAXIAT);
  free_dvector(atom.x ,0,MAXATOM);
  free_dvector(atom.y ,0,MAXATOM);
  free_dvector(atom.z ,0,MAXATOM);
  free_dvector(atom.charge ,0,MAXATOM);
  free_dvector(atom.formal_charge ,0,MAXATOM);
  free_dvector(atom.sigma_charge ,0,MAXATOM);
  free_dvector(atom.atomwt ,0,MAXATOM);
  free_dvector(atom.radius ,0,MAXATOM);
  free(atom.name);

}
/* ================================================================ */
void get_memory()
{
    int i, j;
    int ntor;
    int ntypes, found, itype[MAXATOMTYPE];
    
    skip = imatrix(0,natom+1,0,natom+1);
    for (i=1; i <= natom; i++)
    {
      for (j=1; j <=natom; j++)
	skip[i][j] = 0;
    }

    deriv.d1     = dmatrix(0,natom+1, 0,3);
    deriv.deb    = dmatrix(0,natom+1, 0,3);
    deriv.dea    = dmatrix(0,natom+1, 0,3);
    deriv.destb  = dmatrix(0,natom+1, 0,3);
    deriv.deopb  = dmatrix(0,natom+1, 0,3);
    deriv.detor  = dmatrix(0,natom+1, 0,3);
    deriv.de14  = dmatrix(0,natom+1, 0,3);
    deriv.devdw  = dmatrix(0,natom+1, 0,3);
    deriv.deqq   = dmatrix(0,natom+1, 0,3);
    deriv.deaa   = dmatrix(0,natom+1, 0,3);
    deriv.destor = dmatrix(0,natom+1, 0,3);   
    deriv.dehb = dmatrix(0,natom+1, 0,3);   
    deriv.deimprop = dmatrix(0,natom+1, 0,3);   
    deriv.deub = dmatrix(0,natom+1, 0,3);   
    deriv.desolv = dmatrix(0,natom+1, 0,3);
    deriv.degeom = dmatrix(0,natom+1, 0,3);
    deriv.drb = dvector(0,natom+1);  
    
    ntypes = 0;
    nonbond.maxnbtype = 0;
    
    for (i=1; i <= natom; i++)
    {
        found = FALSE;
        if (atom.type[i] > nonbond.maxnbtype)
           nonbond.maxnbtype = atom.type[i];
        for (j=0; j < ntypes; j++)
        {
            if (atom.type[i] == itype[j])
            {
                found = TRUE;
                break;
            }
        }
        if (found == FALSE)
        {
            itype[ntypes] = atom.type[i];
            ntypes++;
        }
    }

    j=0;
    for (i=0; i <= ntypes; i++)
       j+= i;
       
    nonbond.nonbond = j;    
    nonbond.iNBtype = imatrix(0,nonbond.maxnbtype+1, 0,nonbond.maxnbtype+1);
    nonbond.ipif =    imatrix(0,nonbond.maxnbtype+1, 0,nonbond.maxnbtype+1);
    nonbond.vrad = dmatrix(0,nonbond.maxnbtype+1, 0,nonbond.maxnbtype+1);
    nonbond.veps = dmatrix(0,nonbond.maxnbtype+1, 0,nonbond.maxnbtype+1);
    nonbond.vrad14 = dmatrix(0,nonbond.maxnbtype+1, 0,nonbond.maxnbtype+1);
    nonbond.veps14 = dmatrix(0,nonbond.maxnbtype+1, 0,nonbond.maxnbtype+1);

    // bonds
    bonds_ff.i12 = imatrix(0,MAXBND, 0,2);
    bonds_ff.idpl = ivector(0,MAXBND);
    bonds_ff.index = ivector(0,MAXBND);
    bonds_ff.bk = dvector(0,MAXBND);
    bonds_ff.bl = dvector(0,MAXBND);
    bonds_ff.bmom = dvector(0,MAXBND);
    // angles
    angles.i13 = imatrix(0,MAXANG,0,4);
    angles.index = ivector(0,MAXANG);
    angles.angin = ivector(0,MAXANG);
    angles.angtype = ivector(0,MAXANG);
    angles.acon = vector(0,MAXANG);
    angles.anat = vector(0,MAXANG);
    angles.copb = vector(0,MAXANG);
    // strbnd
   strbnd.isb = imatrix(0,MAXANG, 0,3);
   strbnd.ksb1 = vector(0,MAXANG);
   strbnd.ksb2 = vector(0,MAXANG);
    // torsions
    max_torsions(natom,atom.type,atom.iat,atom.bo);
    ntor = torsions.ntor;
    if (ntor > 0)
    {
       torsions.i14 = imatrix(0,ntor,0,4);
       torsions.v1 = vector(0,ntor);
       torsions.v2 = vector(0,ntor);
       torsions.v3 = vector(0,ntor);
       torsions.v4 = vector(0,ntor);
       torsions.v5 = vector(0,ntor);
       torsions.v6 = vector(0,ntor);
       torsions.ph1 = ivector(0,ntor);
       torsions.ph2 = ivector(0,ntor);
       torsions.ph3 = ivector(0,ntor);
       torsions.ph4 = ivector(0,ntor);
       torsions.ph5 = ivector(0,ntor);
       torsions.ph6 = ivector(0,ntor);
       for (i=0; i < ntor; i++)
       {
           torsions.v1[i] = 0.0;
           torsions.v2[i] = 0.0;
           torsions.v3[i] = 0.0;
           torsions.v4[i] = 0.0;
           torsions.v5[i] = 0.0;
           torsions.v6[i] = 0.0;
           torsions.ph1[i] = 0;
           torsions.ph2[i] = 0;
           torsions.ph3[i] = 0;
           torsions.ph4[i] = 0;
           torsions.ph5[i] = 0;
           torsions.ph6[i] = 0;
       }
    }    
    hess.hessx = matrix(0,natom+1, 0,3);
    hess.hessy = matrix(0,natom+1, 0,3);
    hess.hessz = matrix(0,natom+1, 0,3);

    attached.n13 = ivector(0,natom+1);
    attached.n14 = ivector(0,natom+1);
    attached.i13 = imatrix(0, 20, 0,natom+1);
    attached.i14 = imatrix(0,144, 0,natom+1);

    for(i=0; i <= natom; i++)
    {
        for (j=0; j < 3; j++)
        {
            deriv.d1[i][j] = 0.0;
            deriv.deb[i][j] = 0.0;
            deriv.dea[i][j] = 0.0;
            deriv.destb[i][j] = 0.0;
            deriv.deopb[i][j] = 0.0;
            deriv.detor[i][j] = 0.0;
            deriv.de14[i][j] = 0.0;
            deriv.devdw[i][j] = 0.0;
            deriv.deqq[i][j] = 0.0;
            deriv.deaa[i][j] = 0.0;
            deriv.destor[i][j] = 0.0;
            deriv.dehb[i][j] = 0.0;
            deriv.deimprop[i][j] = 0.0;
            deriv.deub[i][j] = 0.0;

            hess.hessx[i][j] = 0.0;
            hess.hessy[i][j] = 0.0;
            hess.hessz[i][j] = 0.0;
        }
    }
    if (use_solv())
    {
        solvent.asolv = dvector(0,natom+1);
        solvent.rsolv = dvector(0,natom+1);
        solvent.rborn = dvector(0,natom+1);
        if (solvent.type == STILL)
        {
             solvent.vsolv = dvector(0,natom+1);
             solvent.gpol = dvector(0,natom+1);
        } else if (solvent.type == HCT)
        {
             solvent.shct = dvector(0,natom+1);
        }
    }  
}

void free_memory()
{
    int ntor;
    free_imatrix(skip, 0, natom+1,0,natom+1);      
    free_dmatrix(deriv.d1,      0,natom+1, 0,3);
    free_dmatrix(deriv.deb,     0,natom+1, 0,3);
    free_dmatrix(deriv.dea,     0,natom+1, 0,3);
    free_dmatrix(deriv.destb,   0,natom+1, 0,3);
    free_dmatrix(deriv.deopb,   0,natom+1, 0,3);
    free_dmatrix(deriv.detor,   0,natom+1, 0,3);
    free_dmatrix(deriv.de14,    0,natom+1, 0,3);
    free_dmatrix(deriv.devdw,   0,natom+1, 0,3);
    free_dmatrix(deriv.deqq,    0,natom+1, 0,3);
    free_dmatrix(deriv.deaa,    0,natom+1, 0,3);
    free_dmatrix(deriv.destor,  0,natom+1, 0,3);
    free_dmatrix(deriv.dehb,    0,natom+1, 0,3);
    free_dmatrix(deriv.deimprop,0,natom+1, 0,3);
    free_dmatrix(deriv.deub,    0,natom+1, 0,3);
    free_dmatrix(deriv.desolv ,0,natom+1, 0,3);   
    free_dmatrix(deriv.degeom ,0,natom+1, 0,3);   
    free_dvector(deriv.drb ,0,natom+1);  
    
    free_imatrix(nonbond.iNBtype ,0,nonbond.maxnbtype+1, 0,nonbond.maxnbtype+1);
    free_imatrix(nonbond.ipif ,0,nonbond.maxnbtype+1, 0,nonbond.maxnbtype+1);
    free_dmatrix(nonbond.vrad ,0,nonbond.maxnbtype+1, 0,nonbond.maxnbtype+1);
    free_dmatrix(nonbond.veps ,0,nonbond.maxnbtype+1, 0,nonbond.maxnbtype+1);
    free_dmatrix(nonbond.vrad14 ,0,nonbond.maxnbtype+1, 0,nonbond.maxnbtype+1);
    free_dmatrix(nonbond.veps14 ,0,nonbond.maxnbtype+1, 0,nonbond.maxnbtype+1);

    //bonds
    free_imatrix(bonds_ff.i12 ,0,MAXBND, 0,2);
    free_ivector(bonds_ff.idpl ,0,MAXBND);
    free_ivector(bonds_ff.index ,0,MAXBND);
    free_dvector(bonds_ff.bk ,0,MAXBND);
    free_dvector(bonds_ff.bl ,0,MAXBND);
    free_dvector(bonds_ff.bmom ,0,MAXBND);
    // angles
    free_imatrix(angles.i13 ,0,MAXANG,0,4);
    free_ivector(angles.index ,0,MAXANG);
    free_ivector(angles.angin ,0,MAXANG);
    free_ivector(angles.angtype ,0,MAXANG);
    free_vector(angles.acon ,0,MAXANG);
    free_vector(angles.anat ,0,MAXANG);
    free_vector(angles.copb ,0,MAXANG);
    // strbnd
    free_imatrix(strbnd.isb ,0,MAXANG, 0,3);
    free_vector(strbnd.ksb1 ,0,MAXANG);
    free_vector(strbnd.ksb2 ,0,MAXANG);

//    ntor = 9*natom;
    ntor = torsions.ntor;
    if (ntor > 0)
    {
       free_imatrix(torsions.i14 ,0,ntor,0,4);
       free_vector(torsions.v1 ,0,ntor);
       free_vector(torsions.v2 ,0,ntor);
       free_vector(torsions.v3 ,0,ntor);
       free_vector(torsions.v4 ,0,ntor);
       free_vector(torsions.v5 ,0,ntor);
       free_vector(torsions.v6 ,0,ntor);
       free_ivector(torsions.ph1 ,0,ntor);
       free_ivector(torsions.ph2 ,0,ntor);
       free_ivector(torsions.ph3 ,0,ntor);
       free_ivector(torsions.ph4 ,0,ntor);
       free_ivector(torsions.ph5 ,0,ntor);
       free_ivector(torsions.ph6 ,0,ntor);
    }

    free_matrix(hess.hessx, 0,natom+1, 0,3);
    free_matrix(hess.hessy, 0,natom+1, 0,3);
    free_matrix(hess.hessz, 0,natom+1, 0,3);

    free_ivector(attached.n13 ,0,natom+1);
    free_ivector(attached.n14 ,0,natom+1);
    free_imatrix(attached.i13 ,0, 20, 0,natom+1);
    free_imatrix(attached.i14 ,0,144, 0,natom+1);

    if (use_solv())
    {
        free_dvector(solvent.asolv ,0,natom+1);
        free_dvector(solvent.rsolv ,0,natom+1);
        free_dvector(solvent.rborn ,0,natom+1);
        if (solvent.type == STILL)
        {
             free_dvector(solvent.vsolv ,0,natom+1);
             free_dvector(solvent.gpol ,0,natom+1);
        } else if (solvent.type == HCT)
        {
             free_dvector(solvent.shct ,0,natom+1);
        }
    }
}

