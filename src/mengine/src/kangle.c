#define EXTERN extern
#include "pcwin.h"
#include "pcmod.h"

#include "angles.h"
#include "bonds_ff.h"
#include "atom_k.h"

int is_ring43(int, int, int);
int is_ring53(int, int, int);
int is_ring63(int, int, int);
int is_delocalbond(int, int);
int isbond(int,int);
int find_bond(int,int);
long ipack(int ,int ,int ,int);
long kang(int,int,int);
void numeral(int,char *,int);
float get_general_angle(int);
int have_ring3(void);
int have_ring4(void);
int have_ring5(void);
int have_ring6(void);
int get_field(void);
void kangle(void);

EXTERN struct t_angk1 {
        int use_ang3, use_ang4, use_ang5;
        int nang, nang3, nang4, nang5;
        int ndel, ndel3, ndel4;
        char  ktype[MAXANGCONST][10],ktype3[MAXANG3CONST][10],ktype4[MAXANG4CONST][10],ktype5[MAXANG5CONST][10];
        char  kdel[MAXANGDEL][10],kdel3[MAXANG3DEL][10],kdel4[MAXANG4DEL][10];
        int index[MAXANGCONST],index3[MAXANG3CONST],index4[MAXANG4CONST],index5[MAXANG5CONST];
        int indexdel[MAXANGDEL],indexdel3[MAXANG3DEL],indexdel4[MAXANG4DEL];
        float con[MAXANGCONST], ang[MAXANGCONST][3];
        float con3[MAXANG3CONST], ang3[MAXANG3CONST][3];        
        float con4[MAXANG4CONST], ang4[MAXANG4CONST][3];        
        float con5[MAXANG5CONST], ang5[MAXANG5CONST][3];
        float condel[MAXANGDEL],  angdel[MAXANGDEL][3];      
        float condel3[MAXANG3DEL], angdel3[MAXANG3DEL][3];      
        float condel4[MAXANG4DEL], angdel4[MAXANG4DEL][3];      
         } angk1;

EXTERN int Missing_constants;

void kangle()
{
  int i, j, k, izero, ierr, nh, nRc,field;
    int ia, ib, ic, ie;
    int iaa, ibb, icc;
    int cl_a, cl_b, cl_c;
    int nbo1, nb1, nb2;
    char pa[4],pb[4],pc[4],pt[10], pt1[10], pt2[10],pt3[10], pt4[10], pt5[10], pt6[10], pt7[10];
    char zero[4];
    float bkan;

        if( angles.nang != 0 )
        {
	  field = get_field();
           for( i = 0; i < angles.nang; i++ )
           {
               ia = angles.i13[i][0];
               ib = angles.i13[i][1];
               ic = angles.i13[i][2];
               iaa = atom.type[ia];
               ibb = atom.type[ib];
               icc = atom.type[ic];
               cl_a = atom.tclass[ia];
               cl_b = atom.tclass[ib];
               cl_c = atom.tclass[ic];
               
/* get bond index each bond for MMFF94  */
               nb1 = find_bond(ia,ib);
               nb2 = find_bond(ib,ic);
                              
/* special metal cases */
               if( iaa >= 300)
                    cl_a = 300;
                       
               if( ibb >= 300)
                    cl_b = 300;
                        
               if( icc >= 300)
                    cl_c = 300;

/* special MMX cases  */
               if ( field == MMX)
               {                                               
                    if( iaa == 40 )
                         iaa = 2;
                       
                    if( ibb == 40 )
                         ibb = 2;
                        
                    if( icc == 40 )
                         icc = 2;
                         
                }
                numeral(iaa,pa,3);
                numeral(ibb,pb,3);
                numeral(icc,pc,3);
                if( iaa <= icc )
                {
                   strcpy(pt,pa);
                   strcat(pt,pb);
                   strcat(pt,pc);
                }
                if( iaa > icc )
                {
                   strcpy(pt,pc);
                   strcat(pt,pb);
                   strcat(pt,pa);
                }
                izero = 0;
                strcpy(zero,"  0");
                strcpy(pt1,pa); strcat(pt1,pb); strcat(pt1,zero);
                strcpy(pt2,pc); strcat(pt2,pb); strcat(pt2,zero);
                strcpy(pt3,zero); strcat(pt3,pb); strcat(pt3,pc);
                strcpy(pt4,zero); strcat(pt4,pb); strcat(pt4,pa);
                strcpy(pt5,zero); strcat(pt5,pb); strcat(pt5,zero);

                numeral(cl_a,pa,3);
                numeral(cl_b,pb,3);
                numeral(cl_c,pc,3);
                if (cl_a < cl_c)
                {
                   strcpy(pt6,pa); strcat(pt6,pb); strcat(pt6,pc);
                }else
                {
                   strcpy(pt6,pc); strcat(pt6,pb); strcat(pt6,pa);
                }

                if (field == MMFF94)
                {
                    cl_a = atom_k.tclass1[atom.type[ia]];
                    cl_b = atom_k.tclass1[atom.type[ib]];
                    cl_c = atom_k.tclass1[atom.type[ic]];
                    numeral(cl_a,pa,3);
                    numeral(cl_b,pb,3);
                    numeral(cl_c,pc,3);
                    if (cl_a < cl_c)
                    {
                       strcpy(pt7,pa); strcat(pt7,pb); strcat(pt7,pc);
                    }else
                    {
                       strcpy(pt7,pc); strcat(pt7,pb); strcat(pt7,pa);
                    }
                }
                angles.acon[i] = 0.0;
                angles.anat[i] = 0.0;
                angles.angtype[i] = HARMONIC;

                ierr = FALSE;
/*  check delocalized angles in MMFF94 cyclopropanes */
                if (field == MMFF94 && angk1.ndel3 > 0)
                {
                    if (isbond(ia,ic) )
                    {
                        nbo1 = bonds_ff.index[nb1] + bonds_ff.index[nb2];
                        if (nbo1 == 1 || nbo1 == 2)
                        { 
                           for(j=0; j < angk1.ndel3; j++)
                           {
                               if (strcmp(angk1.kdel3[j],pt) == 0)
                               {
                                   angles.acon[i] = angk1.condel3[j];
                                   angles.anat[i] = angk1.angdel3[j][0];
				   if (nbo1 == 1 && angk1.indexdel3[j] == 5)
				     {
					 angles.index[i] = 5;
					 ierr = TRUE;
					 break;
				     } else if (nbo1 == 2 && angk1.indexdel3[j] == 6)
				     {
					 angles.index[i] = 6;
					 ierr = TRUE;
					 break;
				     }
                               }
                           }
                        }
                    
                        if (ierr == TRUE)
                          goto L_10;
                    }
                } 
/*     check three membered rings */
                if ( (have_ring3()) && (angk1.nang3 > 0) )
                {
                    if (isbond(ia,ic))
                    {
/*                         search three membered ring parameters */
                      for (j = 0; j < angk1.nang3; j++)
                      {
                        if ( (strcmp(angk1.ktype3[j],pt) == 0) || (strcmp(angk1.ktype3[j],pt1) == 0) ||
                            (strcmp(angk1.ktype3[j],pt2) == 0) || (strcmp(angk1.ktype3[j],pt3) == 0) ||           
                            (strcmp(angk1.ktype3[j],pt4) == 0) || (strcmp(angk1.ktype3[j],pt5) == 0))
                        {
                           if ( angk1.ang3[j][1] == 0.0 && angk1.ang3[j][2] == 0)
                           {
                             angles.acon[i] = angk1.con3[j];
                             angles.anat[i] = angk1.ang3[j][0];
                             angles.index[i] = 3;
                             ierr = TRUE;
                             break;
                           } else
                           {
                              nh = 0;
                              for (k=0; k < MAXIAT; k++)
                              {
                                ie = atom.iat[ib][k];
                                if ( (atom.type[ie] == 5 || atom.type[ie] == 36) && (ie != ia) && (ie != ic))
                                  nh++;
                              }
                              angles.acon[i] = angk1.con3[j];
                              angles.anat[i] = angk1.ang3[j][nh];
                              angles.index[i] = 3;
                              ierr = TRUE;
                              break;
                           }
                        }
                     }
                     if (ierr == TRUE)
                       goto L_10;
                     else
                     {
                            Missing_constants = TRUE;
                     }
                    }
                }
/*  check delocalized angles in MMFF94 cyclobutanes */
                if (field == MMFF94 && angk1.ndel4 > 0)
                {
                    nRc = is_ring43(ia, ib, ic);
                    if (nRc == TRUE)
                    {
                        nbo1 = bonds_ff.index[nb1] + bonds_ff.index[nb2];
                        for(j=0; j < angk1.ndel4; j++)
                        {
                            if (strcmp(angk1.kdel4[j],pt) == 0)
                            {
                                angles.acon[i] = angk1.condel4[j];
                                angles.anat[i] = angk1.angdel4[j][0];
				if (nbo1 == 1 &&  angk1.indexdel4[j] == 7)
				  {
				    ierr = TRUE;
                                     angles.index[i] = 7;
				     break;
                                  } else if (nbo1 == 2 &&  angk1.indexdel4[j] == 8)
				  {
                                      angles.index[i] = 8;
				      ierr = TRUE;
				      break;
				  }
                               }
			   }
                        if (ierr == TRUE)
                          goto L_10;
                    }
                } 
/*     check four memebered rings */
                if ( (have_ring4()) && (angk1.nang4 > 0) )
                {
                    nRc = is_ring43(ia, ib, ic);
                    if (nRc == TRUE)
                    {
 /*                        search four membered ring parameters */
                      for (j = 0; j < angk1.nang4; j++)
                      {
                        if ( (strcmp(angk1.ktype4[j],pt) == 0) || (strcmp(angk1.ktype4[j],pt1) == 0) ||
                            (strcmp(angk1.ktype4[j],pt2) == 0) || (strcmp(angk1.ktype4[j],pt3) == 0) ||           
                            (strcmp(angk1.ktype4[j],pt4) == 0) || (strcmp(angk1.ktype4[j],pt5) == 0))
                        {
                           if ( angk1.ang4[j][1] == 0.0 && angk1.ang4[j][2] == 0)
                           {
                             angles.acon[i] = angk1.con4[j];
                             angles.anat[i] = angk1.ang4[j][0];
                             angles.index[i] = 4;
                             ierr = TRUE;
                             break;
                           } else
                           {
                              nh = 0;
                              for (k=0; k < MAXIAT; k++)
                              {
                                ie = atom.iat[ib][k];
                                if ( (atom.type[ie] == 5 || atom.type[ie] == 36) && (ie != ia) && (ie != ic))
                                  nh++;
                              }
                              angles.acon[i] = angk1.con4[j];
                              angles.anat[i] = angk1.ang4[j][nh];
                              angles.index[i] = 4;
                              ierr = TRUE;
                              break;
                           }
                        }
                     }
                     if (ierr == TRUE)
                        goto L_10;
                     else
                     {
                           Missing_constants = TRUE;                         
                     }
                    }
                }                    
/*     check five membered rings */
                if ( (have_ring5()) && (angk1.nang5 > 0) )
                {
                    nRc = is_ring53(ia,ib, ic);
                    if (nRc == TRUE)
                    {
 /*                        search five memebered ring parameters */
                      for (j = 0; j < angk1.nang5; j++)
                      {
                        if ( (strcmp(angk1.ktype5[j],pt) == 0) || (strcmp(angk1.ktype5[j],pt1) == 0) ||
                            (strcmp(angk1.ktype5[j],pt2) == 0) || (strcmp(angk1.ktype5[j],pt3) == 0) ||           
                            (strcmp(angk1.ktype5[j],pt4) == 0) || (strcmp(angk1.ktype5[j],pt5) == 0) )
                        {
                           if ( angk1.ang5[j][1] == 0.0 && angk1.ang5[j][2] == 0)
                           {
                             angles.acon[i] = angk1.con5[j];
                             angles.anat[i] = angk1.ang5[j][0];
                             angles.index[i] = 0;
                             ierr = TRUE;
                             break;
                           } else
                           {
                              nh = 0;
                              for (k=0; k < MAXIAT; k++)
                              {
                                ie = atom.iat[ib][k];
                                if ( (atom.type[ie] == 5 || atom.type[ie] == 36) && (ie != ia) && (ie != ic))
                                  nh++;
                              }
                              angles.acon[i] = angk1.con5[j];
                              angles.anat[i] = angk1.ang5[j][nh];
                              angles.index[i] = 0;
                              ierr = TRUE;
                              break;
                           }
                        }
                     }
                     if (ierr == TRUE)
                        goto L_10;
                     // don't fail on missing 5 ring parameters
                   }
                }
/*  check delocalized angles in MMFF94  */
                if (field == MMFF94 && angk1.ndel > 0)
                {
                      nbo1 = bonds_ff.index[nb1] + bonds_ff.index[nb2];
                      for(j=0; j < angk1.ndel; j++)
                      {
                         if (strcmp(angk1.kdel[j],pt) == 0)
                         {
                             angles.acon[i] = angk1.condel[j];
                             angles.anat[i] = angk1.angdel[j][0];
			     if (nbo1 == angk1.indexdel[j])
			       {
				 if (nbo1 == 1)
				   angles.index[i] = 1;
				 else if (nbo1 == 2)
				   angles.index[i] = 2;
				 ierr = TRUE;
				 break;
			       }
                         }
                      }
                      if (ierr == TRUE)
                        goto L_10;

                } 
/*     check regular parameters  for specific angles*/
                for (j = 0; j < angk1.nang; j++)
                {
                    if ( (strcmp(angk1.ktype[j],pt) == 0) )
                    {
                        if ( angk1.ang[j][1] == 0.0 && angk1.ang[j][2] == 0)
                        {
                          angles.acon[i] = angk1.con[j];
                          angles.anat[i] = angk1.ang[j][0];
                          angles.index[i] = 0;
                          ierr = TRUE;
                          break;
                        } else
                        {
                            nh = 0;
                            for (k=0; k < MAXIAT; k++)
                            {
                                ie = atom.iat[ib][k];
                                if ( (atom.atomnum[ie] == 1) && (ie != ia) && (ie != ic))
                                  nh++;
                            }
                            angles.acon[i] = angk1.con[j];
                            angles.anat[i] = angk1.ang[j][nh];
                            angles.index[i] = 0;
                            ierr = TRUE;
                            break;
                        }
                    }
                }
                if (ierr == TRUE)
                   goto L_10;
// did not find any specific look for generalized                         
                for (j = 0; j < angk1.nang; j++)
                {
                    if ( (strcmp(angk1.ktype[j],pt) == 0) || (strcmp(angk1.ktype[j],pt1) == 0)  ||
                         (strcmp(angk1.ktype[j],pt2) == 0) || (strcmp(angk1.ktype[j],pt3) == 0) ||           
                         (strcmp(angk1.ktype[j],pt4) == 0) || (strcmp(angk1.ktype[j],pt5) == 0) ||
                         (strcmp(angk1.ktype[j],pt6) == 0))
                    {
                        if ( angk1.ang[j][1] == 0.0 && angk1.ang[j][2] == 0)
                        {
                          angles.acon[i] = angk1.con[j];
                          angles.anat[i] = angk1.ang[j][0];
                          angles.index[i] = 0;
                          ierr = TRUE;
                          break;
                        } else
                        {
                            nh = 0;
                            for (k=0; k < MAXIAT; k++)
                            {
                                ie = atom.iat[ib][k];
                                if ( (atom.atomnum[ie] == 1) && (ie != ia) && (ie != ic))
                                  nh++;
                            }
                            angles.acon[i] = angk1.con[j];
                            angles.anat[i] = angk1.ang[j][nh];
                            angles.index[i] = 0;
                            ierr = TRUE;
                            break;
                        }
                    }
                }                         
                if (ierr == TRUE)
                   goto L_10;
/*     check MMFF for class1 parameters    */
             if (field == MMFF94)
             {
                for (j = 0; j < angk1.nang; j++)
                {
                    if ( (strcmp(angk1.ktype[j],pt7) == 0) )
                    {
                        if ( angk1.ang[j][1] == 0.0 && angk1.ang[j][2] == 0)
                        {
                          angles.acon[i] = angk1.con[j];
                          angles.anat[i] = angk1.ang[j][0];
                          angles.index[i] = 0;
                          ierr = TRUE;
                          break;
                        } else
                        {
                            nh = 0;
                            for (k=0; k < MAXIAT; k++)
                            {
                                ie = atom.iat[ib][k];
                                if ( (atom.atomnum[ie] == 1) && (ie != ia) && (ie != ic))
                                  nh++;
                            }
                            angles.acon[i] = angk1.con[j];
                            angles.anat[i] = angk1.ang[j][nh];
                            angles.index[i] = 0;
                            ierr = TRUE;
                            break;
                        }
                    }
                }
                if (ierr == TRUE)
                   goto L_10;
             }
             if (ierr != TRUE)
             {
                 bkan = get_general_angle(ibb);
                 angles.acon[i] = 0.80;
                 angles.anat[i] = bkan;
                 angles.index[i] = 0;
	     }
// ====  done with angle lookup =============
L_10:
	     continue;
	   }
	}
}
// ===================================
// generalized angles for MMFF with no hydrogen minimization
float get_general_angle(int ia)
{
    float angs[80] = {
       109.0, 120.0, 120.0, 180.0, 0.000, 109.0, 0.000, 109.0, 120.0, 118.0,
       0.000, 0.000, 0.000, 0.000, 109.0, 109.0, 120.0, 109.0, 109.0, 90.00,
       0.000, 60.00, 0.000, 0.000, 109.0, 109.0, 0.000, 0.000, 0.000, 90.00,
       0.000, 0.000, 0.000, 109.0, 0.000, 0.000, 120.0, 120.0, 120.0, 120.0,
       120.0, 180.0, 120.0, 120.0, 120.0, 120.0, 0.000, 0.000, 109.0, 0.000,
       109.0, 0.000, 180.0, 120.0, 120.0, 120.0, 120.0, 120.0, 120.0, 180.0,
       180.0, 120.0, 120.0, 120.0, 120.0, 120.0, 120.0, 109.0, 120.0, 109.0,
       0.000, 109.0, 109.0, 120.0, 120.0, 109.0, 0.000, 109.0, 109.0, 109.0 };

       if (ia < 80)
         return (angs[ia-1]);
       else
         return 109.0;
}   
/* --------------------------------- */
long kang(int i,int j,int k)
{
        long int kang_v;
        static long i10000 = 10000;
        static long i100 = 100;

        kang_v = i10000*i + i100*j + k;
        return( kang_v );
} 

/* -------------------- */
long ipack(int i,int j,int k,int l)
{
        long int ipack_v;
        static long itbig = 200000000;
        static long itmid = 2000000;
        static long ithou = 1000;

        ipack_v = itbig*i + itmid*j + ithou*k + l;
        return( ipack_v );
}
