#define EXTERN extern

#include "pcwin.h"
#include "pcmod.h"
#include "attached.h"
#include "bonds_ff.h"
 
int isbond(int, int);
void get_bonds(void);
int get_bondorder(int ia,int ib);
/* ================================================================== */
int isbond(int i, int j)
{
    int k;
    for (k=0; k < MAXIAT; k++)
    {
        if (atom.iat[i][k] == j && atom.bo[i][k] != 9)
            return TRUE;
    }
    return FALSE;
}
// ========================
int get_bondorder(int ia,int ib)
{
  int i;
  for (i=0; i < MAXIAT; i++)
    {
      if (atom.iat[ia][i] == ib)
	return (atom.bo[ia][i]);
    }
  return FALSE;
}
// ==========================
void get_bonds()
{
    int i, j;

    bonds_ff.nbnd = 0;    

    for (i=1; i <= natom; i++)
    {
        for (j=i+1; j <= natom; j++)
        {
           if (isbond(i,j))
           {
              bonds_ff.i12[bonds_ff.nbnd][0] = i;
              bonds_ff.i12[bonds_ff.nbnd][1] = j;
              bonds_ff.nbnd++;
              if (bonds_ff.nbnd > MAXBND)
              {
                 fprintf(pcmlogfile,"Error - Too many bonds!\nProgram will now quit.\n");
                 exit(0);
               }
           }
        }
    }
}
// ======================================
// make lists of 1,3 and 1,4 attachments             
void attach()
{
    int i,j,k,m,p;
    int jj,kk,mm, jji;

    for (i=1; i <= natom; i++)
    {
        attached.n13[i] = 0;
        for (j=0; j < MAXIAT; j++)
        {
            jj = atom.iat[i][j];
            jji = 0;
            for (mm = 0; mm < MAXIAT; mm++)
            {
                if (atom.iat[jj][mm] != 0 && atom.bo[jj][mm] != 9)
                   jji++;
            }
            
            if ( atom.type[jj] < 300 && jji < 7 )
            {
               if (jj != 0 && atom.bo[jj][j] != 9)
               {
                  for(k=0; k < MAXIAT; k++)
                  {
                      kk = atom.iat[jj][k];
                      if (kk != i && kk != 0 && atom.bo[jj][k] != 9)
                      {
                         for (m=0; m < MAXIAT; m++)
                         {
                           if (kk == atom.iat[i][m])
                            break;
                         }
                         attached.i13[attached.n13[i]][i] = kk;
                         attached.n13[i]++;
                      }
                  }
               }
            }
        }
    }
  // find 14 relations
    for (i=1; i <= natom; i++)
    {
        attached.n14[i] = 0;
        for (j=0; j < MAXIAT; j++)
        {
            jj = atom.iat[i][j];
            if (jj != 0 && atom.bo[i][j] != 9 && atom.type[jj] < 300)
            {
               for(k=0; k < MAXIAT; k++)
               {
                   kk = atom.iat[jj][k];
                   if (kk != 0 && atom.bo[jj][k] != 9 && atom.type[kk] < 300)
                   {
                      for (m = 0; m < MAXIAT; m++)
                      {
                        mm = atom.iat[kk][m];
                        if (mm != i && mm != 0 && atom.bo[kk][m] != 9)
                        {
                          for (p=0; p < MAXIAT; p++)
                          {
                            if (mm == atom.iat[i][p])
                              goto L_20;
                          }
                          for (p=0; p < attached.n13[i]; p++)
                          {
                             if (mm == attached.i13[p][i])
                               goto L_20;
                          }
                          attached.i14[attached.n14[i]][i] = mm;
                          attached.n14[i]++;
                        }
L_20:
                       continue;
                     }
                   }
               }
            }
        }
    }
}


