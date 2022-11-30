#define EXTERN extern
#include "pcwin.h"
#include "pcmod.h"
#include "atom_k.h"
#include "draw.h"

EXTERN struct t_minim_values {
        int iprint, ndc, nconst;
        float dielc;
        } minim_values;
EXTERN struct t_minim_control {
        int type, method, field, added_const;
        char added_path[256],added_name[256];
        } minim_control;

EXTERN struct {
  int gaff, mmff, amber, opls;
        } AtomTypes[];

EXTERN struct ElementType { char symbol[3];
                             int atomnum;
                             float weight, covradius, vdwradius;
                              int s,p,d,f, type;
                            } Elements[];

EXTERN int gaff_mmx[], mmff_mmx[];

int get_field(void);
int get_hydrid(int);
int is_bond(int,int);
void set_connectivity(void);
int gaff_mmxtype(int);

// ===================
void set_connectivity()
{
  int i,j;
  for (i=1; i <= natom; i++)
    {
      atom.nconnect[i] = 0;
      atom.tbo[i] = 0;
      for (j=0; j < MAXIAT; j++)
	{
	  if (atom.iat[i][j] != 0) // old test included metals and lone pairs
	    {
	      atom.nconnect[i] += 1;
	      atom.tbo[i] += atom.bo[i][j];
	    }
	}
    }
}
// ===============================================
void set_atomdata(int ia, int mmxtype, int gafftype, int mmfftype)
{
  int type,field;
    
  type = 0;
        if (mmxtype > 0)
           atom.mmx_type[ia] = mmxtype;
        if (gafftype > 0)
           atom.gaff_type[ia] = gafftype;
        if (mmfftype > 0)
           atom.mmff_type[ia] = mmfftype;
        field = get_field();

        if (field == MMX || field == MM2)
        {
           type = mmxtype;
        } else if (field == GAFF)
        {
           type = gafftype;
	} else if (field == MMFF94)
        {
            type = mmfftype;
        }
        if ( type < 300) 
        {
           atom.tclass[ia] = atom_k.tclass1[type];
           atom.atomnum[ia] = atom_k.number[type];
           atom.atomwt[ia] = atom_k.weight[type];
          strcpy(atom.name[ia], atom_k.symbol[type]);
        }        
}
// ===============================================
void set_atomtype(int ia, int mmxtype, int gafftype, int mmfftype)
{
        if (mmxtype > 0)
           atom.mmx_type[ia] = mmxtype;
        if (gafftype > 0)
           atom.gaff_type[ia] = gafftype;
        if (mmfftype > 0)
           atom.mmff_type[ia] = mmfftype;           
}
// ======================================================
int make_atom(int type, float x, float y, float z,char *name)
{
  int i,iz,field;
     int mmxtype, gafftype, mmfftype;

     natom++;
     if (natom >= MAXATOM)
     {
         natom--;
         message_alert("Max atoms exceeded in makeatom","PCMODEL Error");
         return -1;
     }
     field = get_field();
        iz = strlen(name);
        if (iz > 0) // got an atom name use it 
        {
            for (i=1; i <= atom_k.natomtype; i++)
            {
                if (strcmp(name,atom_k.symbol[i]) == 0) // names match
                {
                    strcpy(atom.name[natom],name);
                    atom.tclass[natom] = atom_k.tclass1[i];
                    atom.atomwt[natom] = atom_k.weight[i];
                    atom.atomnum[natom] = atom_k.number[i];
                    atom.type[natom] = atom_k.type[i];
                    if (type == 0)
                    {
                        type = atom.type[natom];  // this is an mmff type
                        mmxtype = mmff_mmxtype(type);
                        gafftype = 0;
                        mmfftype = type;
                        set_atomtype(natom,mmxtype,gafftype,mmfftype);
                        if (mmfftype == 0) return -1;
                    } else if (field == MMX || field == MM2)
                    {
                          atom.mmx_type[natom] = type;
                          mmxtype = type;
                          gafftype = AtomTypes[type-1].gaff;
                          mmfftype = AtomTypes[type-1].mmff;
                          set_atomtype(natom,mmxtype,gafftype,mmfftype);
                    } else if (field == GAFF)
                    {
                          atom.gaff_type[natom] = type;
                          mmxtype = gaff_mmxtype(type);
                          gafftype = type;
                          mmfftype = AtomTypes[mmxtype-1].mmff;
                          set_atomtype(natom,mmxtype,gafftype,mmfftype);
                    } else if (field == MMFF94)
                    {
                          atom.mmff_type[natom] = type;
                          mmxtype = mmff_mmxtype(type);
                          gafftype = AtomTypes[mmxtype-1].gaff;
                          mmfftype = type;
                          set_atomtype(natom,mmxtype,gafftype,mmfftype);
                    }

                    atom.x[natom] = x; atom.y[natom] = y; atom.z[natom] = z;
                    atom.flags[natom] = 0;
                    return natom;
                } 
            }
           for( i = 0; i < 103; i++ )
           {
              if( strcasecmp(Elements[i].symbol,name) == 0 )
              {
                  strcpy(atom.name[natom],Elements[i].symbol);
                  atom.atomwt[natom] = Elements[i].weight;
                  atom.atomnum[natom] = Elements[i].atomnum;
                  atom.type[natom] = Elements[i].type;    // set generic type
                  atom.tclass[natom] = Elements[i].type;    // set generic type
                  atom.mmx_type[natom] = Elements[i].type;    // set generic type
                  if (type > 0) 
                  {
                       atom.type[natom] = type; // overwrite type if type is passed in
                       if (field == MMX || field == MM2)
                          atom.mmx_type[natom] = type;
                       else if (field == GAFF)
                          atom.gaff_type[natom] = type;
                       else if (field == MMFF94)
                          atom.mmff_type[natom] = type;
                  }
                  if (type != 0) atom.type[natom] = type;
                  atom.x[natom] = x; atom.y[natom] = y; atom.z[natom] = z;
                  if (type > 299)
                    return -1;
                  else if (atom.mmff_type[natom] == 0)
                    return -1;
                  else
                   return natom;
              }
            }
        } else
        {
           message_alert("Input atom does not have valid atomic symbol or atom type","Error");
           return -1;
        }
        return -1;
}
// ====================================================
void make_bond(int ia1, int ia2, int bo)
{
    int loop,loop1;

   if (ia1 < 1 || ia1 > MAXATOM) return;
   if (ia2 < 1 || ia2 > MAXATOM) return;
   /*   find if bond exists        */
   for (loop = 0; loop < atom.nconnect[ia1]; loop++)
   {
     if (atom.iat[ia1][loop] == ia2)
     {
        if (bo == 1)
           atom.bo[ia1][loop] += bo;
        else
           atom.bo[ia1][loop] = bo;
        
        for (loop1 = 0; loop1 < atom.nconnect[ia2]; loop1++)
        {
         if (atom.iat[ia2][loop1] == ia1)
         {
             if ( bo == 1)
             {
                atom.bo[ia2][loop1] += bo; /* bond exists return */
                return;
             } else
             {
                atom.bo[ia2][loop1] = bo; /* bond exists return */
                return;
             }
         }
        }
      }
   }
/*  bond does not exist create it               */     
   for (loop = 0; loop < MAXIAT; loop++) 
   {
      if (atom.iat[ia1][loop] == 0) 
      {
          atom.iat[ia1][loop] = ia2;
          atom.bo[ia1][loop] = bo;
	  atom.nconnect[ia1]++;
          break;
      }
   }
   for (loop = 0; loop < MAXIAT; loop++) 
   {
      if (atom.iat[ia2][loop] == 0) 
      {
          atom.iat[ia2][loop] = ia1;
          atom.bo[ia2][loop] = bo;
	  atom.nconnect[ia2]++;
          break;
      }
   }
}
// ===============================================
void deleteatom(int i)
{
   int j,k, ibo;
   atom.type[i] = 0;
   atom.x[i] = 0.0F;
   atom.y[i] = 0.0F;
   atom.z[i] = 0.0F;
   atom.flags[i] = 0;
   atom.nconnect[i] = 0;
   for (j=0; j<MAXIAT; j++)
   {
      if (atom.iat[i][j] != 0 )
      {
         if (atom.atomnum[atom.iat[i][j]] == 1 || atom.atomnum[atom.iat[i][j]] == 0 )
         {
             atom.type[atom.iat[i][j]] = 0;
             atom.x[atom.iat[i][j]] = 0.0F;
             atom.y[atom.iat[i][j]] = 0.0F;
             atom.z[atom.iat[i][j]] = 0.0F;
             atom.flags[atom.iat[i][j]] = 0;
	     atom.nconnect[atom.iat[i][j]] = 0;
         }
         ibo = atom.bo[i][j];
         for (k=1; k <= ibo; k++)
            deletebond(i,atom.iat[i][j]);
         atom.iat[i][j] = 0;
     }
   }
   natom--;
}
// ==================================================
void deletebond(int i, int j) 
{
   int i1,j1;
   for (i1=0; i1 < atom.nconnect[i]; i1++) 
   {
      if (atom.iat[i][i1] == j) 
      {
        if (atom.bo[i][i1] == 9)  /* if deleting a coordinated bond - delete it completely */
        {
            atom.bo[i][i1] = 0;
            atom.iat[i][i1] = 0;
        }
        atom.bo[i][i1] -= 1;
        if (atom.bo[i][i1] <= 0 )  // completely delete bond
	  {
            atom.iat[i][i1] = 0;
	    atom.nconnect[i] = 0;
	  }
        break;
      }
   }
   for (j1=0; j1 < atom.nconnect[j]; j1++) 
   {
      if (atom.iat[j][j1] == i) 
      {
        if (atom.bo[j][j1] == 9)  /* if deleting a coordinated bond - delete it completely */
        {
            atom.bo[j][j1] = 0;
            atom.iat[j][j1] = 0;
        }
        atom.bo[j][j1] -= 1;
        if (atom.bo[j][j1] <= 0)
	  {
            atom.iat[j][j1] = 0;
	    atom.nconnect[j] = 0;
	  }
        break;
      }
   }
}
// ========================================================
int gaff_mmxtype(int gafftype)
{
    int i;

    if (gafftype < 58)
       i = gaff_mmx[gafftype-1];
    else
       i = gafftype;
    return(i);
}
// ==================================================
int mmff_mmxtype(int mmfftype)
{
    int i;
    if (mmfftype < 100)
       i = mmff_mmx[mmfftype-1];
    else
       i = mmfftype;
    return(i);
}
/* =============================================== */
// get hybridization of atom - tetrahedral = 1, planar = 2, linear = 3
//
int get_hybrid(int ia)
{
    int itype;

    itype = atom.mmx_type[ia];
    if (itype == 2 || itype == 3 || itype == 7 || itype == 9 || itype == 25 ||
        itype == 26 || itype == 29 || itype == 30 || itype == 37 || itype == 38 ||
        itype == 40 || itype == 57 )
        return 2;
    else if (itype == 4 || itype == 10)
        return 3;
    else
        return 1;
}
/* ------------------------------------*/    
int is_bond(int ia, int ib)
{
    int i;
    for (i=0; i < atom.nconnect[ia]; i++)
    {
        if (atom.iat[ia][i] == ib)
          return TRUE;
    }
    return FALSE;
}
