#define EXTERN extern

#include "pcwin.h"
#include "pcmod.h"
#include "bonds_ff.h"
#include "atom_k.h"

int is_ring31(int);
int is_ring41(int);
int is_ring51(int);
int is_ring61(int);
int is_ring42(int, int);
int is_ring52(int, int);
int is_ring62(int,int);
void get_ring62(int,int,int *);
void numeral(int, char *, int);
int is_delocalbond(int ia, int ib);
void message_alert(char *, char *);
double SIGN(double ,double  ); 
int is_bond(int,int);
float get_bond(int, int);
int have_ring3(void);
int have_ring4(void);
int have_ring5(void);
int have_ring6(void);
int get_field(void);
char *get_structure_title(void);
int kbond(void);
int is_ring_aromatic6(int *);

EXTERN struct  t_bondk1 {
        int use_ring3, use_ring4, use_ring5;
        int nbnd, nbnd3, nbnd4, nbnd5, ndeloc;
        char kb[MAXBONDCONST][7], kb3[MAXBOND3CONST][7],kb4[MAXBOND4CONST][7],kb5[MAXBOND5CONST][7];
        char kbdel[MAXBONDDELOC][7];
        float s[MAXBONDCONST], t[MAXBONDCONST];
        float s3[MAXBOND3CONST], t3[MAXBOND3CONST];
        float s4[MAXBOND4CONST], t4[MAXBOND4CONST];
        float s5[MAXBOND5CONST], t5[MAXBOND5CONST];
        float sdel[MAXBONDDELOC], tdel[MAXBONDDELOC];
        }  bondk1;
                        
EXTERN int Missing_constants;

int kbond()
{
    long int pi_mask;
    int field;
    int i, j,k, iit, kit, ierr,nbpi;
    int igeni, igenk;
    int ia, ib;
    char hatext[80];
    char pa[4],pb[4], pt[7], pt_gen[7];

    pi_mask = 1L << PI_MASK;   /*  this mask is for pi atoms - now uses bit 0 */
    nbpi = 0;
    if( natom != 0 )
    {
      field = get_field();
      for( i = 0; i < bonds_ff.nbnd; i++ )
      {
         ia = bonds_ff.i12[i][0];
         ib = bonds_ff.i12[i][1];
         iit = atom.type[bonds_ff.i12[i][0]];
         kit = atom.type[bonds_ff.i12[i][1]];

         igeni = atom_k.tclass1[iit];
         igenk = atom_k.tclass1[kit];

// metal general class
         if (iit > 300)
           igeni = 300;
         if (kit > 300)
           igenk = 300;

         if( field == MMX)
         {
             if (iit == 40 && ( kit != 2 && kit != 3 && kit != 4 && kit != 40) )
                iit = 2;
             if (kit == 40 && ( iit != 2 && iit != 3 && iit != 4 && iit != 40) )
                kit = 2;
         }

         numeral(iit, pa, 3);
         numeral(kit, pb, 3);
         strcpy(pt,pa);
         strcat(pt,pb);
         if( iit <= kit )
         {
            strcpy(pt,pa);
            strcat(pt,pb);
         }else
         {
            strcpy(pt,pb);
            strcat(pt,pa);
         }
         numeral(igeni, pa, 3);
         numeral(igenk, pb, 3);
         if( igeni <= igenk )
         {
            strcpy(pt_gen,pa);
            strcat(pt_gen,pb);
         }else
         {
            strcpy(pt_gen,pb);
            strcat(pt_gen,pa);
         }

       ierr = FALSE;
//  check 3 membered ring
       if ((have_ring3()) && bondk1.nbnd3 > 0 && ( is_ring31(ia) && is_ring31(ib) ) )
      {
        for (j=0; j < bondk1.nbnd3; j++)
        {
          if ( (strcmp(bondk1.kb3[j], pt) == 0)  || (strcmp(bondk1.kb3[j],pt_gen) == 0) )
          {
              bonds_ff.bk[i] = bondk1.s3[j];
              bonds_ff.bl[i] = bondk1.t3[j];
              bonds_ff.index[i] = 0;
              ierr = TRUE;
              break;
          }
         }
         if (ierr != TRUE)
         {
             Missing_constants = TRUE;
             sprintf(hatext,"Bond constants missing for cyclopropane bond %d - %d of type %d-%d\n",bonds_ff.i12[i][0],bonds_ff.i12[i][1],
             atom.type[bonds_ff.i12[i][0]], atom.type[bonds_ff.i12[i][1]]);
             fprintf(pcmlogfile,"%s\n",hatext);
             //fprintf(stderr,"%s\n",hatext);
         }
      }
      if (ierr == TRUE)
       goto L_10;
//  check 4 membered ring
      if ( (have_ring4()) && bondk1.nbnd4 > 0 && ( is_ring41(ia) && is_ring41(ib) ) )
      {
         for (j=0; j < bondk1.nbnd4; j++)
         {
          if ( (strcmp(bondk1.kb4[j],pt) == 0)  || (strcmp(bondk1.kb4[j],pt_gen) == 0) )
          {
              bonds_ff.bk[i] = bondk1.s4[j];
              bonds_ff.bl[i] = bondk1.t4[j];
              bonds_ff.index[i] = 0;
              ierr = TRUE;
              break;
          }
         }
         if (ierr != TRUE)
         {
            Missing_constants = TRUE;
            sprintf(hatext,"Bond constants missing for cyclobutane bond %d - %d of type %d-%d\n",bonds_ff.i12[i][0],bonds_ff.i12[i][1],
             atom.type[bonds_ff.i12[i][0]], atom.type[bonds_ff.i12[i][1]]);
            fprintf(pcmlogfile,"%s\n",hatext);
            //            fprintf(stderr,"%s\n",hatext);
         }
      }
      if (ierr == TRUE)
       goto L_10;
//  check 5 membered ring
      if ( (have_ring5()) && bondk1.nbnd5 > 0 && (is_ring51(ia) != -1) && (is_ring51(ib) != -1))
      {
         for (j=0; j < bondk1.nbnd5; j++)
         {
          if (strcmp(bondk1.kb5[j],pt) == 0 )
          {
              bonds_ff.bk[i] = bondk1.s5[j];
              bonds_ff.bl[i] = bondk1.t5[j];
              bonds_ff.index[i] = 0;
              ierr = TRUE;
              break;
          }
         }
         // don't fail on missing param - check for regular param first
     }
      if (ierr == TRUE)
       goto L_10;
//  delocalized bonds in MMFF94
      if (field == MMFF94 && bondk1.ndeloc > 0)
      {
          if ( is_delocalbond(ia, ib) )
          {
             for (k = 0; k < MAXIAT; k++)
             {
                 if (atom.iat[ia][k] == ib)
                 {
                    if (atom.bo[ia][k] == 1) // single bond
                    {
                      for (j=0; j < bondk1.ndeloc; j++)
                      {
                          if (strcmp(bondk1.kbdel[j],pt) == 0)
                          {
                              bonds_ff.bk[i] = bondk1.sdel[j];
                              bonds_ff.bl[i] = bondk1.tdel[j];
                              bonds_ff.index[i] = 1;
                              ierr = TRUE;
                              break;
                          }
                      }
                    }
                 }
              }
          }
          if (ierr == TRUE)
            goto L_10;
      }
//  regular parameters
      for (j=0; j < bondk1.nbnd; j++)
      {
          if ( (strcmp(bondk1.kb[j], pt) == 0) )
          {
              bonds_ff.bk[i] = bondk1.s[j];
              bonds_ff.bl[i] = bondk1.t[j];
              bonds_ff.index[i] = 0;
              ierr = TRUE;
              if ( ((strcmp(pt,"  2  2") == 0) || (strcmp(pt,"  2  3") == 0))&& field == MMX)
              {
                  for (k=0; k < MAXIAT; k++)
                  {
                      if (atom.iat[ia][k] == ib)
                      {
                          if (atom.bo[ia][k] == 1)  // butadiene
                          {
                              bonds_ff.bl[i] = 1.48;
                              break;
                          }
                      }
                  }
              }
              break;
          }
      }
      if (ierr == TRUE)
         goto L_10;
//  
//  generalized parameters
      for (j=0; j < bondk1.nbnd; j++)
      {
          if ( (strcmp(bondk1.kb[j], pt_gen) == 0) )
          {
              bonds_ff.bk[i] = bondk1.s[j];
              bonds_ff.bl[i] = bondk1.t[j];
              bonds_ff.index[i] = 0;
              ierr = TRUE;
              break;
          }
      }
      if (ierr == TRUE)
         goto L_10;      
      if (ierr != TRUE)
      {
          Missing_constants = TRUE;
          sprintf(hatext,"Bond constants missing for bond %d - %d of type %d - %d %s %s\n",bonds_ff.i12[i][0],bonds_ff.i12[i][1],
		  atom.type[bonds_ff.i12[i][0]], atom.type[bonds_ff.i12[i][1]],pt,pt_gen);
          fprintf(pcmlogfile,"%s\n",hatext);
      }
L_10:
      continue;
          
    }
  }
  if (Missing_constants == TRUE)
  {
    fprintf(pcmlogfile,"Error assigning constants in: %s\n",get_structure_title());
      fprintf(pcmlogfile,"%s\n",hatext);
      return FALSE;
  }
          
  return TRUE;     
}
/* ------------------------------------------  */
int find_bond(int ia, int ib)
{
    int i;
    int iz;
    
    iz = -1;
    for (i=0; i < bonds_ff.nbnd; i++)
    {
        if (ia == bonds_ff.i12[i][0] && ib == bonds_ff.i12[i][1])
           return(i);
        if (ib == bonds_ff.i12[i][0] && ia == bonds_ff.i12[i][1])
           return(i);
    }
    return(iz);
}

/* ------------------------------------------  */
int is_delocalbond(int ia, int ib)
{
    // test for delocalized bond
    // if ia & ib are aromatic and part of 6 membered ring
    // if not aromatic is bond order 1

    long int aromatic_mask;
    int i, j;
    int jdbl, kdbl;
    int array[10];

    aromatic_mask = (1L << AROMATIC_MASK);

    if ( (atom.flags[ia] & aromatic_mask) && (atom.flags[ib] & aromatic_mask) )
      {
	if (is_ring62(ia,ib))
	  {
	    get_ring62(ia,ib,array);
	    if (is_ring_aromatic6(array))
	      return(FALSE);
	  }
      }

    if ( (atom.flags[ia] & aromatic_mask) && (atom.flags[ib] & aromatic_mask)
         && (is_ring52(ia,ib) == TRUE) )
      {
         return(FALSE);
      }
   if (atom.type[ia] == 37 && atom.type[ib] == 37 && (is_ring62(ia,ib) != FALSE))
      return FALSE;
   if (atom.type[ia] == 37 && atom.type[ib] == 37 && (is_ring52(ia,ib) != FALSE))
      return FALSE;

   for (i=0; i < MAXIAT; i++)
   {
       if (atom.iat[ia][i] == ib)
       {
           if (atom.bo[ia][i] == 1)
           {
             jdbl = FALSE;
             kdbl = FALSE;
             for (j=0; j < MAXIAT; j++)
             {
                 if (atom.bo[ia][j] >= 2 && atom.bo[ia][j] != 9)
                   jdbl = TRUE;
                 if (atom.bo[ib][j] >= 2 && atom.bo[ib][j] != 9)
                   kdbl = TRUE;
             }
             if (jdbl == TRUE && kdbl == TRUE)
                return (TRUE);
             else if (jdbl == TRUE && (atom.flags[ib] & aromatic_mask) )
                return (TRUE);
             else if (kdbl == TRUE && (atom.flags[ia] & aromatic_mask) )
                return (TRUE);            
             else
                return(FALSE);
           }
       }
   }
   return(FALSE);
}
// ==============================
int is_ring_aromatic6(int *array)
{
  int i;
  long int aromatic_mask;
  aromatic_mask = (1L << AROMATIC_MASK);
  for (i=0; i < 6; i++)
    {
      if (!(atom.flags[array[i]] & aromatic_mask))
	return FALSE;
    }
  return TRUE;  // all atoms of ring are aromatic
}
// ==================
#define IS_ODD(j)       ((j) & 1 )

double SIGN(double a,double b )         /* floating polong int SIGN transfer */
{
        return( b < 0.0 ? -fabs(a) : fabs(a) );
}

