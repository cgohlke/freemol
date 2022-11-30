#define EXTERN extern
#include "pcwin.h"
#include "pcmod.h"
#include "atom_k.h"
#include "draw.h"

int get_field(void);
void hadd(void);
void hdel(int);
void type(void);
void nhadd(void);
void  hcoord(void);
static void revec(void);
static void  xyplan(float *, float *, float *);
static void xaxis(float *,float *,float *);
static void getvec(void);
static void reseq(void);
void set_atomtypes(int);

static double orient[4][3];        
EXTERN struct t_minim_control {
        int type, method, field, added_const;
        char added_path[256],added_name[256];
        } minim_control;
        
void hadd(void)
{
  int field;
  field = get_field();
    type();
    if (field == GAFF)
    {
        nhadd();
        hcoord();
    } else if (field == MMFF94)
    {
        nhadd();
        hcoord();
    } else
    {
        set_atomtypes(MMX);
        nhadd();
        hcoord();
    }
    type();
}
// ========================
void nhadd()
{
  int i, i1, ik, im, it, itads, j, jk, nhadds, newatom,field;
        int jji, bo,ncount;
        static int iadd[]={ 0,0,0,0,0,3,1,6,0,1,
                                  0,0,0,0,3,6,3,0,0,0,
                                  0,0,0,0,6,0,0,0,0,0,
                                  0,0,0,3,0,0,6,1,1,0,
                                  0,4,0,0,0,6,0,6,0,0,
                                  3,3,3,0,6,0,0,0,0,0,
                                  0,0,0,0,0,1,6,0,0,0,
                                  0,0,0,0,0};

        /*     this routine calculates the #'s of h to be added
         *     to carbon and packs it symbolically into iat(i,j) table
         *     as well as into the connection table via the appropriate
         *     manipulation instructions.  hydrogen coordinates
         *     are not added by this routine. lp are also added to oxygen. */
        /*     0 to disregard : 1 and 4 for lp's only
         *     2 and 5 for h's only and 3 and 6 for h and lp
         *     add 3 for trigonal atoms */
	field = get_field();
   ncount = natom;
   for( i = 1 ; i <= ncount; i++ )
   {
        if( atom.mmx_type[i] > 0 && (atom.mmx_type[i] < 300) && atom.atomnum[i] != 1 ) 
        {
           jji = 0;
           bo = 0;
           for (j=0; j < MAXIAT; j++)
           {
              if (atom.iat[i][j] != 0 && atom.bo[i][j] != 9)
              {
                 jji ++;
                 bo += atom.bo[i][j];
              }
           }
           nhadds = atom_def.ligands[atom.mmx_type[i]] - jji;

           if (atom.mmx_type[i] >= 11 && atom.mmx_type[i] <= 14)
               nhadds = 0;
           if (atom.mmx_type[i] == 41 || atom.mmx_type[i] == 46)
             nhadds = atom_def.ligands[atom.mmx_type[i]] - bo;
           if (atom.mmx_type[i] == 49 || atom.mmx_type[i] == 50 || atom.mmx_type[i] == 51)
             nhadds = atom_def.ligands[atom.mmx_type[i]] - bo;
           
/*  add atom and bond it to atom i      */
           if( nhadds > 0 )  
           {
              for( im = 1; im <= nhadds; im++ )
              {
                  newatom = make_atom(5,0.F,0.F,0.F,"H");
                  set_atomdata(newatom,5,5,5);
                  if (newatom == -1)
                  {
                      message_alert("Error in nhadd newatom = -1","Error");
                      hdel(0);
                      return;
                  }
                  make_bond(i,newatom,1);
              }
           }
        }
   }

/*       ***  check for lp's ***        */
   if (field == MMX)
   {     
      ncount = natom;
      for( i = 1 ; i <= ncount; i++ )
      {
           it = iadd[atom.mmx_type[i]-1];         
           if( it >= 1 )
           {
               itads = atom.mmx_type[i] ;
               if (it == 1 || it == 3 || it == 4 || it == 6 ) 
               {
                   nhadds = 4;
                   for( i1 = 0; i1 < MAXIAT; i1++ )
                   {
                       if( atom.bo[i][i1] != 0 && atom.bo[i][i1] != 9 )
                           nhadds = nhadds - atom.bo[i][i1];
                   }
                   if( itads == 6 || itads == 15 || itads == 34 || itads == 42 || itads == 53 || itads == 48 )
                   {
                        for( ik = 0; ik < MAXIAT; ik++ )
                        {
                            if( atom.iat[i][ik] != 0 && atom.bo[i][ik] != 9 )
                            {
                                jk = atom.mmx_type[atom.iat[i][ik]];
                                if( jk == 2 || jk == 3 || jk == 4 || jk == 9
                                 || jk == 10 || jk == 30 || jk == 37 || jk == 40 )
                                {
                                    nhadds = nhadds - 1;
                                    break;
                                }
                            }
                        }
                   }
                                
/*         correction for sulfoxide */
                   if( itads == 17 )
                       nhadds = nhadds + 1;

                   if (itads == 66)
                       nhadds = 2;

                   if( nhadds > 0 ) 
                   {
                      for( i1 = 1; i1 <= nhadds; i1++ )
                      {
                         newatom = make_atom(20, 0.F,0.F,0.F,"");
                         set_atomdata(newatom,20,0,0);
                         make_bond(i,newatom,1);
                      }
                   }
               }
           }
        }
   }
   return;
} 
// ==================================
void  hcoord()
{
        int i, i1, i3, i4, i5, i6, i7, i8, 
         i9, iadj[10], ii, ii1, ii2, ii3, iii, 
         ij, ijk[3], it, j, j1, j2, j3, j4, j5, j8,  
         j9, k1, k2, l1, na, nhadds, ni;
        float adj, dis, xpnt, ypnt, zpnt;

        getvec();
        for( i = 1; i <= natom; i++ )
        {
           it = atom.mmx_type[i];
           if( it < 300)
           {
              for( j = 0; j < 10; j++ )
                   iadj[j] = 0;

              for( j = 0; j < 3; j++ )
                   ijk[j] = 0;

              nhadds = 0;
              l1 = 0;
              for( j = 0; j < MAXIAT; j++ )
              {
                  if( atom.iat[i][j] != 0 && atom.bo[i][j] !=  9 )
                  {
                      if( atom.mmx_type[atom.iat[i][j]] == 5 || atom.mmx_type[atom.iat[i][j]] == 20 )
                           nhadds++ ;
                      else
                      {
                           iadj[l1] = atom.iat[i][j];
                           l1++;
                      }
                  }
              }
              if( nhadds != 0 )
              {
/*       **** calculate the coordinates of the center we
 *       are going to move to the x-axis */
                 xpnt = atom.x[i];
                 ypnt = atom.y[i];
                 zpnt = atom.z[i];
                 for( ii = 1; ii <= (natom + 4); ii++ )
                 {
                       atom.x[ii] = atom.x[ii] - xpnt;
                       atom.y[ii] = atom.y[ii] - ypnt;
                       atom.z[ii] = atom.z[ii] - zpnt;
                 }
                 xpnt = 0.F;
                 ypnt = 0.F;
                 zpnt = 0.F;
                 adj = 0.F;
                 for( i1 = 0; i1 < 10; i1++ )
                 {
                     if( iadj[i1] != 0 )
                     {
                         na = iadj[i1];
                         adj += 1.0001F ;
                         xpnt += atom.x[na];
                         ypnt += atom.y[na];
                         zpnt += atom.z[na];
                     }
                 }
                 if( adj > 0.01F )
                 {
                      xpnt /= adj;
                      ypnt /= adj;
                      zpnt /= adj;
                      xaxis( &xpnt, &ypnt, &zpnt );
                      if( xpnt > 0.F )
                      {
                         for( i3 = 1; i3 <= (natom + 4); i3++ )
                                atom.x[i3] = -atom.x[i3];
                      }
/*       ****  after moving, check the # of h added to the carbon i */
                      if( nhadds == 1 )
                      {
/*       **** there is only 1 h added to the carbon i */
                         for( i4 = 0; i4 < MAXIAT; i4++ )
                         {
                            if( atom.iat[i][i4] != 0 && atom.bo[i][i4] != 9 )
                            {
                                i5 = atom.iat[i][i4];
                                if( atom.mmx_type[i5] == 5 || atom.mmx_type[i5] == 20 )
                                {
/*           i5 is the atom to attach */
                                    if (adj < 2.0)  // fix for OH in MMFF
                                    {
                                        atom.x[i5] = 0.70;
                                        atom.y[i5] = 0.70;
                                        atom.z[i5] = 0.00;
                                    } else
                                    {
                                        atom.x[i5] = 1.10F;
                                        if( atom.mmx_type[i5] == 20 )
                                          atom.x[i5] = 0.6F;
                                        atom.y[i5] = 0.00F;
                                        atom.z[i5] = 0.00F;
                                    }
                                    if( atom.atomnum[i] == 16 )
                                    {
                                        atom.x[i5] = 0.F;
                                        atom.y[i5] = 0.3796F;
                                        atom.z[i5] = 1.043F;
                                    }
                                }
                            }
                         }
                      }
/*       **** there are 2 h added to the carbon i */
                      if( nhadds == 2 )
                      {
                          if( it == 2 || it == 3  || it == 29 || it == 30 || it == 9 || it == 41 || it == 46
                                      || it == 48 || it == 7 || it == 37 || it == 39
                                      || it == 42 || it == 38 || it == 40 || it == 44 || it == 66)
                          {
                               for( i8 = 0; i8 < 10; i8++ )
                               {
                                   if( iadj[i8] != 0 )
                                   {
                                       ii1 = iadj[i8];
                                       for( ii2 = 0; ii2 < MAXIAT; ii2++ )
                                       {
                                           ii3 = 0;
                                           if( atom.iat[ii1][ii2] != 0 && atom.iat[ii1][ii2] != i && atom.bo[ii1][ii2] != 9 )
                                           {
                                                  ii3 = atom.iat[ii1][ii2];
                                                  xpnt = atom.x[ii3];
                                                  ypnt = atom.y[ii3];
                                                  zpnt = atom.z[ii3];
                                                  xyplan( &xpnt, &ypnt, &zpnt );
                                                  ij = 0;
                                                  for( i9 = 0; i9 < MAXIAT; i9++ )
                                                  {
                                                      if( atom.iat[i][i9] != 0 && atom.bo[i][i9] != 9 )
                                                      {
                                                           iii = atom.iat[i][i9];
                                                           if( atom.mmx_type[iii] == 5 ||  atom.mmx_type[iii] == 20 )
                                                           {
                                                                ijk[ij] = iii;
                                                                ij++;
                                                                if( ij == 2 )
                                                                {
              /*           *** =ch2 *** */
                                                                      atom.x[ijk[0]] = 0.5582F;
                                                                      atom.x[ijk[1]] = 0.5582F;
                                                                      atom.y[ijk[0]] = 0.9478F;
                                                                      atom.y[ijk[1]] = -0.9478F;
                                                                      atom.z[ijk[0]] = 0.00F;
                                                                      atom.z[ijk[1]] = 0.00F;
                                                                      if( it == 7 || it == 38 || it ==  42 || it == 39 || it == 66)
                                                                      {
                                                                          atom.x[ijk[0]] = 0.30F;
                                                                          atom.x[ijk[1]] = 0.30F;
                                                                          atom.y[ijk[0]] = 0.52F;
                                                                          atom.y[ijk[1]] = -0.52F;
                                                                      }
                                                                      if( atom.mmx_type[i] == 37 )
                                                                      {
                                                                         atom.x[ijk[1]] = 0.35F;
                                                                         atom.y[ijk[1]] = -0.5F;
                                                                      }
                                                                }
                                                           }
                                                      }
                                                  }
                                                  break;
                                           }
                                           if (ii3 == 0) // C=O with no lone pairs
                                           {
                                                  ij = 0;
                                                  for( i9 = 0; i9 < MAXIAT; i9++ )
                                                  {
                                                      if( atom.iat[i][i9] != 0 && atom.bo[i][i9] != 9 )
                                                      {
                                                           iii = atom.iat[i][i9];
                                                           if( atom.mmx_type[iii] == 5 ||  atom.mmx_type[iii] == 20 )
                                                           {
                                                                ijk[ij] = iii;
                                                                ij++;
                                                                if( ij == 2 )
                                                                {
              /*           *** =ch2 *** */
                                                                      atom.x[ijk[0]] = 0.5582F;
                                                                      atom.x[ijk[1]] = 0.5582F;
                                                                      atom.y[ijk[0]] = 0.9478F;
                                                                      atom.y[ijk[1]] = -0.9478F;
                                                                      atom.z[ijk[0]] = 0.00F;
                                                                      atom.z[ijk[1]] = 0.00F;
                                                                      if( it == 7 || it == 38 || it ==  42 || it == 39 || it == 66)
                                                                      {
                                                                          atom.x[ijk[0]] = 0.30F;
                                                                          atom.x[ijk[1]] = 0.30F;
                                                                          atom.y[ijk[0]] = 0.52F;
                                                                          atom.y[ijk[1]] = -0.52F;
                                                                      }
                                                                      if( atom.mmx_type[i] == 37 )
                                                                      {
                                                                         atom.x[ijk[1]] = 0.35F;
                                                                         atom.y[ijk[1]] = -0.5F;
                                                                      }
                                                                }
                                                           }
                                                      }
                                                  }
                                                  break;
                                           }
                                       }
                                   }
                               }
                          } else
                          {
/*       ** there are 2 atoms attached to carbon i,
 *       one of which is "ni" and move these two atoms to xy-plane */
                              for( i6 = 0; i6 < 10; i6++ )
                              {
                                  if( iadj[i6] != 0 )
                                  {
                                      ni = iadj[i6];
                                      xpnt = atom.x[ni];
                                      ypnt = atom.y[ni];
                                      zpnt = atom.z[ni];
                                      xyplan( &xpnt, &ypnt, &zpnt );
                                  }
                              }
                              ij = 0;
                              for( i7 = 0; i7 < MAXIAT; i7++ )
                              {
                                  if( atom.iat[i][i7] != 0 && atom.bo[i][i7] != 9 )
                                  {
                                      iii = atom.iat[i][i7];
                                      if( atom.mmx_type[iii] == 5 || atom.mmx_type[iii] == 20 )
                                      {
                                          ijk[ij] = iii;
                                          ij++ ;
                                          if( ij == 2 )
                                          {
                  /*       *** -ch2- type 1 atom *** */
                                              atom.x[ijk[0]] = 0.6758F;
                                              atom.x[ijk[1]] = 0.6758F;
                                              atom.y[ijk[0]] = 0.00F;
                                              atom.y[ijk[1]] = 0.00F;
                                              atom.z[ijk[0]] = 0.8807F;
                                              atom.z[ijk[1]] = -0.8807F;
                                              if( atom.mmx_type[i] == 25 ||  atom.mmx_type[i] == 46 ||  atom.mmx_type[i] == 67)
                                              {
                                                  atom.x[ijk[0]] = 0.3796F;
                                                  atom.x[ijk[1]] = 0.3796F;
                                                  atom.y[ijk[0]] = 1.043F;
                                                  atom.y[ijk[1]] = -0.5215F;
                                                  atom.z[ijk[0]] = 0.0000F;
                                                  atom.z[ijk[1]] = 0.9032F;
                                              }
                                              if( atom.mmx_type[i] == 22 ||  atom.mmx_type[i] == 30 )
                                              {
                 /*         *** -ch2- type 22 and 3o atoms *** */
                                                   atom.x[ijk[0]] = 0.2F;
                                                   atom.x[ijk[1]] = 0.2F;
                                                   atom.y[ijk[0]] = 0.0F;
                                                   atom.y[ijk[1]] = 0.0F;
                                                   atom.z[ijk[0]] = 1.0F;
                                                   atom.z[ijk[1]] = -1.0F;
                                              }
                 /*       *** -o(lp)2- *** */
                                              if (atom.mmx_type[i] == 6 && (atom.mmx_type[iadj[0]] == 3 ||
                                                  atom.mmx_type[iadj[0]] == 2 || atom.mmx_type[iadj[0]] == 30 ||
                                                  atom.mmx_type[iadj[0]] == 37 || atom.mmx_type[iadj[0]] == 40) )
                                              {
                                                   atom.x[ijk[0]] = 0.450F;
                                                   atom.x[ijk[1]] = 0.2488F;
                                                   atom.z[ijk[0]] = 0.0F;
                                                   atom.z[ijk[1]] = 0.0F;
                                                   atom.y[ijk[0]] = 0.9500F;
                                                   atom.y[ijk[1]] = -0.5460F;
                                              }
                                              else if( atom.mmx_type[i] == 6 || atom.mmx_type[i] == 15 || atom.mmx_type[i] == 
                                                  42 || atom.mmx_type[i] == 34 || atom.mmx_type[i] == 35 || atom.mmx_type[i] == 53 )
                                              {
                                                   atom.x[ijk[0]] = 0.2488F;
                                                   atom.x[ijk[1]] = 0.2488F;
                                                   atom.y[ijk[0]] = 0.0F;
                                                   atom.y[ijk[1]] = 0.0F;
                                                   atom.z[ijk[0]] = 0.5460F;
                                                   atom.z[ijk[1]] = -0.5460F;
                                              }                                                 
                                          }
                                      }
                                  }
                              }
                          }
                      }
           /*  end of nhadds = 2
                                         *         **** there are 3 h added to the carbon i *** */
                      if( nhadds == 3 )
                      {
                         for( j1 = 0; j1 < 10; j1++ )
                         {
                             if( iadj[j1] != 0 )
                             {
                                j2 = iadj[j1];
                                if( atom.atomnum[j2] != 1 && atom.mmx_type[j2] != 20 )
                                     na = j2;
                                if( atom.mmx_type[j2] == 2 || atom.mmx_type[j2] == 3 )
                                {
                                    for( j3 = 0; j3 < MAXIAT; j3++ )
                                    {
                                         if( (atom.iat[j2][j3] != 0 && atom.iat[j2][j3] !=  i) &&
                                         atom.bo[j2][j3] !=  9 )
                                         {
                                             j4 = atom.iat[j2][j3];
                                             if( atom.mmx_type[j4] == 2 || (atom.mmx_type[j4] == 7 && atom.mmx_type[j4] ==  38) )
                                             {
                                                 xpnt = atom.x[j4];
                                                 ypnt = atom.y[j4];
                                                 zpnt = atom.z[j4];
                                                 xyplan( &xpnt, &ypnt, &zpnt );
                                                 if( atom.y[j4] < 0.F )
                                                 {
                                                     for( j5 = 1; j5 <= (natom + 4); j5++ )
                                                         atom.y[j5] = -atom.y[j5];
                                                 }
                                             }
                                         }
                                    }
                                }else
                                {
 /*  end of type 2 and type 3 adjacent* */
                                   for( j8 = 0; j8 < MAXIAT; j8++ )
                                   {
                                       if( (atom.iat[j2][j8] != 0 && atom.iat[j2][j8] != i)
                                        && atom.bo[j2][j8] != 9 )
                                        {
                                            j9 = atom.iat[j2][j8];
//                                            if( atom[j9].mmx_type != 5 &&  atom[j9].mmx_type != 20 )
                                            if(atom.mmx_type[j9] != 20 )
                                            {
                                                xpnt = atom.x[j9];
                                                ypnt = atom.y[j9];
                                                zpnt = atom.z[j9];
                                                xyplan( &xpnt, &ypnt, &zpnt );
                                                if( atom.y[j9] > 0.F )
                                                {
                                                    for( k1 = 1; k1 <= (natom + 4); k1++ )
                                                          atom.y[k1] = -atom.y[k1];
                                                }
                                            }
                                        }
                                   }
                                }
/* do type 4 here */
                                ij = 0;
                                for( k2 = 0; k2 < MAXIAT; k2++ )
                                {
                                    if( !(atom.iat[i][k2] == 0 || atom.bo[i][k2] == 9) )
                                    {
                                        iii = atom.iat[i][k2];
                                        if( !(atom.mmx_type[iii] != 5 &&  atom.mmx_type[iii] != 20) )
                                        {
                                            ijk[ij] = iii;
                                            ij++;
                                            if( ij == 3 )
                                            {
                                                dis = 0.F;
                                              /*       *** -ch3 ***
                                               *       *** add the hydrogens displaced along x enough
                                               *       to give a normal length c-c bond *** */
                                                atom.x[ijk[0]] = 0.3796F + dis;
                                                atom.x[ijk[1]] = 0.3796F + dis;
                                                atom.x[ijk[2]] = 0.3796F + dis;
                                                atom.y[ijk[0]] = 1.043F;
                                                atom.y[ijk[1]] = -0.5215F;
                                                atom.y[ijk[2]] = -0.5215F;
                                                atom.z[ijk[0]] = 0.00F;
                                                atom.z[ijk[1]] = 0.9032F;
                                                atom.z[ijk[2]] = -0.9032F;
                                            }
                                        }
                                    }
                                }
                             }
                         }
                      }
                 }
              }
           }
        }
        /*     *** reorient the structure in starting orientation and
         *     turn off atom markers *** */
        revec();
        return;
}
// ========================================================
void revec()
{
        int i;
        float xpnt, ypnt, zpnt;
        /*     *** this routine used the unit vectors and origin saved
         *     by 'subroutine getvec' to return the molecule to
         *     its original orientation *** */

            for( i = 1; i <= natom ; i++ )
            {
	      atom.x[i] = atom.x[i] - orient[3][0]; // atom.x[natom + 4];
	      atom.y[i] = atom.y[i] - orient[3][1]; //atom.y[natom + 4];
	      atom.z[i] = atom.z[i] - orient[3][2]; //atom.z[natom + 4];
            }
            xpnt = orient[0][0]; //atom.x[natom + 1];
            ypnt = orient[0][1]; //atom.y[natom + 1];
            zpnt = orient[0][2]; //atom.z[natom + 1];
            xaxis( &xpnt, &ypnt, &zpnt );
            xpnt = orient[1][0]; //atom.x[natom + 2];
            ypnt = orient[1][1]; //atom.y[natom + 2];
            zpnt = orient[1][2]; //atom.z[natom + 2];
            xyplan( &xpnt, &ypnt, &zpnt );
            if( orient[0][0] < 0. )
            {
               for( i = 1; i <= natom; i++ )
                     atom.x[i] = -atom.x[i];
            }
            if( orient[1][1] < 0. )
            {
               for( i = 1; i <= natom; i++ )
                     atom.y[i] = -atom.y[i];
            }
            if( orient[2][2] < 0. )
            {
               for( i = 1; i <= natom; i++ )
                     atom.z[i] = -atom.z[i];
            }
        return;
} 
/*     ------------------------ */
void  xyplan(float *xr, float *yr, float *zr)
{
        int i;
        float cos3, denom, savex, sin3;

        /*     this subroutine rotates the molecule to place atom (xr,yr,*zr)
         *     in the x,y-plane.  note0 y-coordinate will be +. */

        denom = sqrt( *yr**yr + *zr**zr );
        if( fabs( denom ) >= .00001F )
        {
           sin3 = *zr/denom;
           cos3 = *yr/denom;
           for( i = 1; i <= natom; i++ )
           {
               savex = atom.y[i];
               atom.y[i] = atom.y[i]*cos3 + atom.z[i]*sin3;
               atom.z[i] = atom.z[i]*cos3 - savex*sin3;
           }
	   for (i= 0; i < 4; i++)  // move unit vectors
	   {
	       savex = orient[i][1];
	       orient[i][1] = orient[i][1]*cos3 + orient[i][2]*sin3;
	       orient[i][2] = orient[i][2]*cos3 - savex*sin3;
	   }

           savex = *yr;
           *yr = *yr*cos3 + *zr*sin3;
           *zr = *zr*cos3 - savex*sin3;
        }
        return;
}
/*     ------------------------ */
void xaxis(float *xpnt,float *ypnt,float *zpnt)
{
        int i;
        float cos1, cos2, denom, save, sin1, sin2;

        /*     this subroutine rotates the molecule about z-axis and then y-axis to
         *     place point with coordinates (i*xpnt,iypnt,izpnt) along x-axis. */

        denom = sqrt( *xpnt**xpnt + *ypnt**ypnt );
        if( fabs( denom ) < .00001F )
        {
                sin1 = 0.F;
                cos1 = 1.F;
        }else
        {
                sin1 = *ypnt/denom;
                cos1 = *xpnt/denom;
                for( i = 1; i <= natom; i++ )
                {
                        save = atom.x[i];
                        atom.x[i] = atom.x[i]*cos1 + atom.y[i]*sin1;
                        atom.y[i] = atom.y[i]*cos1 - save*sin1;
                }
		for (i= 0; i < 4; i++) // rotate unit vectors
		  {
		    save = orient[i][0];
		    orient[i][0] = orient[i][0]*cos1 + orient[i][1]*sin1;
		    orient[i][1] = orient[i][1]*cos1 - save*sin1;
		  }
                save = *xpnt;
                *xpnt = *xpnt*cos1 + *ypnt*sin1;
                *ypnt = *ypnt*cos1 - save*sin1;
        }
        denom = sqrt( *xpnt**xpnt + *zpnt**zpnt );
        if( fabs( denom ) >= .00001F )
        {
                sin2 = *zpnt/denom;
                cos2 = *xpnt/denom;
                for( i = 1; i <= natom; i++ )
                {
                        save = atom.x[i];
                        atom.x[i] = atom.x[i]*cos2 + atom.z[i]*sin2;
                        atom.z[i] = atom.z[i]*cos2 - save*sin2;
                }
		for (i= 0; i < 4; i++)  // rotate unit vectors
		  {
		    save = orient[i][0];
		    orient[i][0] = orient[i][0]*cos2 + orient[i][2]*sin2;
		    orient[i][2] = orient[i][2]*cos2 - save*sin2;
		  }
                save = *xpnt;
                *xpnt = *xpnt*cos2 + *zpnt*sin2;
                *zpnt = *zpnt*cos2 - save*sin2;
        }
        return;
}
/* ==============================================  */
void getvec()
{

        /*     *** this routine saves unit vectors and origin in the
         *     n+1 - n+4 slots of x, y, and z *** */
  // new version uses array orient to save unit vectors

  orient[0][0] = 1.0F;
  orient[0][1] = 0.0F;
  orient[0][2] = 0.0F;
  orient[1][0] = 0.0F;
  orient[1][1] = 1.0F;
  orient[1][2] = 0.0F;
  orient[2][0] = 0.0F;
  orient[2][1] = 0.0F;
  orient[2][2] = 1.0F;
  orient[3][0] = 0.0F;
  orient[3][1] = 0.0F;
  orient[3][2] = 0.0F;
} 
/* ==============================================  */
void hdel(int lptest)
{
/*  lptest = 0 to remove H and lone pair  */
/*  lptest = 1 to remove lone pairs only  */
/*  lptest = 2 to remove lp & reset various groups for amber  */
   int i, iatt, iatt1, iatyp, it, it1,ntemp;
        /*     *** delete hydrogens *** */
   ntemp = natom;
   for( i = 1; i <= ntemp; i++ )
   {
       if( atom.mmx_type[i] != 0 )
       {
          iatyp = atom.mmx_type[i];
          if( (lptest == 0 && (iatyp == 5 || iatyp == 20 )) || (lptest >= 1 && iatyp == 20) ) 
          {
              iatt = atom.iat[i][0];
              iatt1 = atom.iat[i][1];
              if( iatt != 0 || iatt1 != 0)  
              {
                 it = atom.mmx_type[iatt];
                 it1 = atom.mmx_type[iatt1];
                 if (it >= 11 && it <= 14)   // halogens
                    goto L_10;
                 if( (it < 300) && (it1 < 300))   // attached metals
                    deleteatom(i);
               }
           }
         }
L_10:
         continue;
   }
     reseq();
     return;
} 
/* --------------------------  */
void reseq()
{
  int i, ia, ib, iplus, ixx, j, jincr, k,ncount;

        /*     ***this subroutine will repack the atom
         *     ***connection table, making sure to fill in
         *     **all "holes" caused by deleting atoms or bonds.
         *     ***the necessary common blocks involve arrays
         *     of the connection table (contb) */
    for( i = 1; i <= natom; i++ )
    {
        if( atom.type[i] != 0 )
        {
           for( ia = 0; ia < (MAXIAT - 1); ia++ )
           {
               if( atom.iat[i][ia] == 0 )
               {
                                        /*            *** for every zero entry in the array iat
                                         *             for atom i, want to move down all remaining
                                         *            entries of iat by one slot; after doing so
                                         *             want to store that zero entry in the last
                                         *            slot--iat *** */
                   for( ib = ia + 1; ib < MAXIAT; ib++ )
                   {
                       atom.iat[i][ib - 1] = atom.iat[i][ib];
                       atom.bo[i][ib - 1] = atom.bo[i][ib];
                   }
                   atom.iat[i][MAXIAT - 1] = 0;
               }
           }
        }
    }

/*     *** now want to renumber atoms to leave
 *     no "holes" caused by the deleted atoms *** */
     jincr = 0;
/*  i points to first zero entry
 *  kincr points to next real entry */
     for( i = 1; i <= (natom - 1); i++ )
     {
        iplus = 0;
        if( atom.type[i] == 0 )
        {
            for( jincr = i; jincr <= natom; jincr++ )
            {
                if( atom.type[jincr] != 0 )
                {
                    iplus = jincr;
                    break;
                }
            }
            if (iplus == 0)
                 break;
            for( j = 0; j < MAXIAT; j++ )
            {
                atom.iat[i][j] = atom.iat[iplus][j];
                atom.bo[i][j] = atom.bo[iplus][j];
            }

            atom.x[i] = atom.x[iplus];
            atom.y[i] = atom.y[iplus];
            atom.z[i] = atom.z[iplus];
            atom.type[i] = atom.type[iplus];
            atom.tclass[i] = atom.tclass[iplus];
            atom.mmx_type[i] = atom.mmx_type[iplus];
            atom.gaff_type[i] = atom.gaff_type[iplus];
            atom.mmff_type[i] = atom.mmff_type[iplus];
            atom.atomnum[i] = atom.atomnum[iplus];
            atom.formal_charge[i] = atom.formal_charge[iplus];
            atom.atomwt[i] = atom.atomwt[iplus];
            atom.use[i] = atom.use[iplus];
            strcpy(atom.name[i],atom.name[iplus]);
            atom.charge[i] = atom.charge[iplus];
            atom.flags[i] = atom.flags[iplus];
            atom.radius[i] = atom.radius[iplus];
        /*           *** delete the atom whose statistics were just
         *           stored in the i-th atom's connection table *** */
            for( ixx = 0; ixx < MAXIAT; ixx++ )
            {
                atom.iat[iplus][ixx] = 0;
                atom.bo[iplus][ixx] = 0;
            }
            atom.x[iplus] = 0.0F;
            atom.y[iplus] = 0.0F;
            atom.z[iplus] = 0.0F;
            atom.type[iplus] = 0;
            atom.tclass[iplus] = 0;
            atom.mmx_type[iplus] = 0;
            atom.gaff_type[iplus] = 0;
            atom.mmff_type[iplus] = 0;
            atom.atomnum[iplus] = 0;
            atom.formal_charge[iplus] = 0;
            atom.atomwt[iplus] = 0.0F;
            atom.use[iplus] = 0;             
            atom.charge[iplus] = 0.0F;
            atom.flags[iplus] = 0;
            atom.radius[iplus] = 0.0F;
            strcpy(atom.name[iplus],"");
 
            for( j = 0; j < MAXIAT; j++ )
            {
                if( atom.iat[i][j] != 0 )
                {
                    for( k = 0; k < MAXIAT; k++ )
                    {
                        if( atom.iat[atom.iat[i][j]][k] == iplus )
                             atom.iat[atom.iat[i][j]][k] = i;
                    }
                }
            }
        }
     }
/*     ***now want to do the same resequencing for iat array */

/* recalc number of atoms */
     ncount = natom;
     natom = 0;
     for (i = 1; i <= ncount; i++)
     {
        if (atom.type[i] != 0)
           natom++;
     }
     return;
}
