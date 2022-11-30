#define EXTERN extern
#include "pcwin.h"
#include "pcmod.h"
#include "angles.h"
#include "bonds_ff.h"

int find_bond(int, int);
void numeral(int, char *, int);
int get_field(void);
void kstrbnd(void);

EXTERN struct t_crossterm_k {
        int nangang, nstrbnd, nstrtor;
        int ang_ang[MAXAA], stbnindex[MAXSTBN] ;
        char str_tor[MAXSTRTOR][7], stbn[MAXSTBN][10];
        float aacon[MAXAA][3], stbncon[MAXSTBN][3], str_torcon[MAXSTRTOR];
        } crossterm_k;

struct t_strbnd {
        int nstrbnd, **isb;
        float *ksb1, *ksb2;
        } strbnd;
struct stbndat {
        int irow, krow, lrow;
        float con1, con2;
      } stbn[] =
      {
       { 0,    1,    0,      0.15,      0.15 }, 
       { 0,    1,    1,      0.10,      0.30 }, 
       { 0,    1,    2,      0.05,      0.35 }, 
       { 0,    1,    3,      0.05,      0.35 }, 
       { 0,    1,    4,      0.05,      0.35 }, 
       { 0,    2,    0,      0.00,      0.00 }, 
       { 0,    2,    1,      0.00,      0.15 }, 
       { 0,    2,    2,      0.00,      0.15 }, 
       { 0,    2,    3,      0.00,      0.15 }, 
       { 0,    2,    4,      0.00,      0.15 }, 
       { 1,    1,    1,      0.30,      0.30 }, 
       { 1,    1,    2,      0.30,      0.50 }, 
       { 1,    1,    3,      0.30,      0.50 }, 
       { 1,    1,    4,      0.30,      0.50 }, 
       { 2,    1,    2,      0.50,      0.50 }, 
       { 2,    1,    3,      0.50,      0.50 }, 
       { 2,    1,    4,      0.50,      0.50 }, 
       { 3,    1,    3,      0.50,      0.50 }, 
       { 3,    1,    4,      0.50,      0.50 }, 
       { 4,    1,    4,      0.50,      0.50 }, 
       { 1,    2,    1,      0.30,      0.30 }, 
       { 1,    2,    2,      0.25,      0.25 }, 
       { 1,    2,    3,      0.25,      0.25 }, 
       { 1,    2,    4,      0.25,      0.25 }, 
       { 2,    2,    2,      0.25,      0.25 }, 
       { 2,    2,    3,      0.25,      0.25 }, 
       { 2,    2,    4,      0.25,      0.25 }, 
       { 3,    2,    3,      0.25,      0.25 }, 
       { 3,    2,    4,      0.25,      0.25 }, 
       { 4,    2,    4,      0.25,      0.25 }};

EXTERN int Missing_constants;
//EXTERN FILE *errfile;

void kstrbnd()
{
  int i, j, nb1, nb2, irow, krow, lrow,field;
   int ia, ib, ic, itb, ita, itc, ierr;
   int k, nstbn, index,nbo1,nbo2;
   char pa[4],pb[4],pc[4],pt[10],pt1[10];

   index = lrow = krow = irow = 0;
   strbnd.nstrbnd = 0;
   nb1 = -1;
   nb2 = -1;
   field = get_field();
   for (i=0; i < angles.nang; i++)
   {
       ierr = FALSE;
       ia = angles.i13[i][0];
       ib = angles.i13[i][1];
       ic = angles.i13[i][2];
       ita = atom.type[ia];
       itb = atom.type[ib];
       itc = atom.type[ic];
       if (field == MMX && (atom.type[ia] >= 300 || atom.type[ib] >= 300 || atom.type[ic] >= 300))
          goto L_10;
       if (field == MMX && (atom.type[ia] >= 20 || atom.type[ic] >= 20))
          goto L_10;
       nstbn = -1;
       numeral(ita,pa,3);
       numeral(itb,pb,3);
       numeral(itc,pc,3);
       if (ita <= itc)
       {
           strcpy(pt,pa);     strcat(pt,pb);  strcat(pt,pc);
       }else
       {
           strcpy(pt,pc);     strcat(pt,pb);  strcat(pt,pa);
       }
       itb = atom.tclass[ib];
       numeral(itb,pb,3);
       strcpy(pt1,"  0"); strcat(pt1,pb); strcat(pt1,"  0");
       for (j=0; j < crossterm_k.nstrbnd; j++)
       {
           if (strcmp(crossterm_k.stbn[j],pt) == 0 || strcmp(crossterm_k.stbn[j],pt1) == 0 )
           {
               nstbn = j;
               ierr = TRUE;
               break;
           }
       }
       if (ita < itc)
	 {
	   nb1 = find_bond(ia,ib);
	   nb2 = find_bond(ib,ic);
	 } else
	 {
	   nb2 = find_bond(ia,ib);
	   nb1 = find_bond(ib,ic);
	 }
       nbo1 = bonds_ff.index[nb1];
       nbo2 = bonds_ff.index[nb2];
       index = 0;

       if (angles.index[i] == 0)
          index = 0;
       else if (angles.index[i] == 1)
       {
          if (nbo1 == 1)
             index = 1;
	  else if (nbo2 == 1)
	     index = 2;
       } else if (angles.index[i] == 2)
       {
           if (nbo1 == 1 && nbo2 == 1)
             index = 3;
       } else if (angles.index[i] == 3)
          index = 5;
       else if (angles.index[i] == 4)
          index = 4;
       else if (angles.index[i] == 5)
       {
          if (nbo1 == 1)
             index = 6;
          else if (nbo2 == 1)
             index = 7;
       } else if (angles.index[i] == 6)
          index = 8;
       else if (angles.index[i] == 7)
       {
          if (nbo1 == 1)
             index = 9;
          else
             index = 10;
       } else if (angles.index[i] == 8)
          index = 11;

       if (crossterm_k.stbnindex[nstbn] != index)
	 {
          nstbn = -1;
	 }
       if (field == MMFF94)
       {
           ierr = TRUE;
           if (angles.anat[i] < 175.0)
           {
              if (nstbn == -1)  // constant not found  assign constants base on periodic table
              {
                 if (atom.atomnum[ia] <= 1)
                  irow = 0;
                 else if (atom.atomnum[ia] <= 10)
                  irow = 1;
                 else if (atom.atomnum[ia] <= 18)
                  irow = 2;
                 else if (atom.atomnum[ia] <= 36)
                  irow = 3;
                 else if (atom.atomnum[ia] <= 54)
                  irow = 4;
                  
                 if (atom.atomnum[ib] <= 1)
                  krow = 0;
                 else if (atom.atomnum[ib] <= 10)
                  krow = 1;
                 else if (atom.atomnum[ib] <= 18)
                  krow = 2;
                 else if (atom.atomnum[ib] <= 36)
                  krow = 3;
                 else if (atom.atomnum[ib] <= 54)
                  krow = 4;
                  
                 if (atom.atomnum[ic] <= 1)
                  lrow = 0;
                 else if (atom.atomnum[ic] <= 10)
                  lrow = 1;
                 else if (atom.atomnum[ic] <= 18)
                  lrow = 2;
                 else if (atom.atomnum[ic] <= 36)
                  lrow = 3;
                 else if (atom.atomnum[ic] <= 54)
                  lrow = 4;
                for (j=0; j < 30; j++)
                {
                  if ( irow == stbn[j].irow && krow == stbn[j].krow && lrow == stbn[j].lrow )
                  {
		    //		 printf("irow: %d %d %d - %d %d %d - %d %d - %d %d\n",ia,ib,ic,irow,krow,lrow,ita,itb,nb1,nb2);
                      strbnd.isb[strbnd.nstrbnd][0] = i;
		      if (ita < itc)
			{
			  strbnd.isb[strbnd.nstrbnd][1] = nb1;
			  strbnd.isb[strbnd.nstrbnd][2] = nb2;
			} else
			{
			  strbnd.isb[strbnd.nstrbnd][1] = nb2;
			  strbnd.isb[strbnd.nstrbnd][2] = nb1;
			}
		      strbnd.ksb1[strbnd.nstrbnd] = stbn[j].con1;
		      strbnd.ksb2[strbnd.nstrbnd] = stbn[j].con2;
                      strbnd.nstrbnd++;
                      break;
                  } else if (lrow == stbn[j].irow && krow == stbn[j].krow && irow == stbn[j].lrow)
                  {
                      strbnd.isb[strbnd.nstrbnd][0] = i;
		      //	 printf("lrow: %d %d %d - %d %d %d - %d %d - %d %d\n",ia,ib,ic,irow,krow,lrow,ita,itb,nb1,nb2);
		      if (ita < itc)
			{
			  strbnd.isb[strbnd.nstrbnd][1] = nb1;
			  strbnd.isb[strbnd.nstrbnd][2] = nb2;
			} else
			{
			  strbnd.isb[strbnd.nstrbnd][1] = nb2;
			  strbnd.isb[strbnd.nstrbnd][2] = nb1;
			}
		      strbnd.ksb1[strbnd.nstrbnd] = stbn[j].con2;
		      strbnd.ksb2[strbnd.nstrbnd] = stbn[j].con1;
                      strbnd.nstrbnd++;
                      break;
                  }
                }   
             }else
             {
                strbnd.isb[strbnd.nstrbnd][0] = i;
                if (angles.anat[i] > 175.0)
                {
                   strbnd.ksb1[strbnd.nstrbnd] = 0.0;
                   strbnd.ksb2[strbnd.nstrbnd] = 0.0;
                } else
                { 
                 if (ita < itc)
                 {
		   strbnd.isb[strbnd.nstrbnd][1] = nb1;
		   strbnd.isb[strbnd.nstrbnd][2] = nb2;
                   strbnd.ksb1[strbnd.nstrbnd] = crossterm_k.stbncon[nstbn][0];
                   strbnd.ksb2[strbnd.nstrbnd] = crossterm_k.stbncon[nstbn][1];
                 } else
                 {
		   strbnd.isb[strbnd.nstrbnd][1] = nb2;
		   strbnd.isb[strbnd.nstrbnd][2] = nb1;
                   strbnd.ksb1[strbnd.nstrbnd] = crossterm_k.stbncon[nstbn][1];
                   strbnd.ksb2[strbnd.nstrbnd] = crossterm_k.stbncon[nstbn][0];
                 }
                }
                strbnd.nstrbnd++;
             }
           }
       } else
       {         
             k = 0;
             if ( atom.atomnum[ia] <= 1)
             {
                 k++;
                 if (crossterm_k.stbncon[nstbn][2] == 0.0)
                   nb1 = -1;
             }
             if (atom.atomnum[ic] <= 1)
             {
                 k++;
                 if (crossterm_k.stbncon[nstbn][2] == 0.0)
                    nb2 = -1;
             }
             if (crossterm_k.stbncon[nstbn][k] != 0.0)
             {
                 strbnd.isb[strbnd.nstrbnd][0] = i;
                 strbnd.isb[strbnd.nstrbnd][1] = nb1;
                 strbnd.isb[strbnd.nstrbnd][2] = nb2;
                 strbnd.ksb1[strbnd.nstrbnd] = crossterm_k.stbncon[nstbn][k];           
                 strbnd.ksb2[strbnd.nstrbnd] = crossterm_k.stbncon[nstbn][k];           
                 strbnd.nstrbnd++;
             }
       }
//       if (ierr == FALSE)
//       {
//           Missing_constants = TRUE;
//           fprintf(errfile,"Missing strbnd constants for %d %d %d \n",ia,ib,ic);
//       }
L_10:
    continue;
   }
}
           
       
