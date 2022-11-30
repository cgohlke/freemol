#define EXTERN extern
#include "pcwin.h"
#include "pcmod.h"


void set_atomtype(int,int,int,int);
void set_atomtypes(int);
int is_ring31(int);
int is_ring41(int);
int is_ring51(int);
int is_ring61(int);
int is_cyclo5(int, int *);
int isbond(int,int);
int aromatic_5(int *,long int *,int *,int *,int *,int **,int **);
int find_rsize(int,int);
void get_rsize(int,int,int, int *);
int icompare(int, int *, int *);
void adjust_mmfftypes(void);
void deletebond(int,int);
int get_field(void);
int count_ewg(int);
void type_mmx(void);

// ==================
void type_mmx()
{
  int i, j, ij, jji, jjk, jjbo, jj_bo, iatype,nh,newg;
   int k, l, jatm, katm, latm, jbo,kbo,lbo, ismetal;
   int ktype,jtype;
   int  adjn, nplus, adjn1, noxide, icurr;
   int mmxtype,mm3type,mmfftype, gafftype, ia, ib, ia1, ib1;
   int ndouble, ntriple;
   int icycl3, icycl4, icycl5, icycl6;
   int full_ring=0, non_pi;
   int nox,nox_double,nnit,nsulf,nnh,nc;
   int array[7];
   int nit[6];
   long int aromatic_mask, mask6, type_mask;

      
   aromatic_mask = (1 << AROMATIC_MASK);
   mask6 = (1L << RING6);
   type_mask = (1L << NO_RETYPE);
   ia = 0;
   ib = 0;
   mm3type = 0;
   latm = katm = jatm = 0;
   k = 0;
   jji = 0;

   for (i=1; i <= natom; i++)
   {
      mmxtype = atom.mmx_type[i];
      gafftype = atom.gaff_type[i]; 
      mmfftype = atom.mmff_type[i];
      if (atom.atomnum[i] != 0 && !(atom.flags[i] & type_mask) )
      {
          if (mmxtype == 100 || mmxtype == 101 || mmxtype == 102 || mmxtype == 103 ) // user defined type 
          {
              goto L_20;
          }
          if (atom.atomnum[i] == 1) // hydrogens
          {
              if (atom.mmx_type[i] == 45)  // TS H 
              {
                  mmxtype = 45;
                  mm3type = 0;
                  mmfftype = 0;
                 goto L_10;
              }
              if (atom.mmx_type[i] == 36) // Deuterium
              {
                  mmxtype = 36;
                  mm3type = 36;
                  mmfftype = 5;
                   goto L_10;
              }                
              jji = atom.nconnect[i];
              if (jji == 2 ) // hydrogen bonded to two atoms or bonded to one and coordinated to another
              {
                  mmxtype = 70;
                  goto L_10;
              }
              if (atom.atomnum[atom.iat[i][0]] == 6)   // carbon
              {
                  mmxtype = 5;
                  mm3type = 5;
                  mmfftype = 5;
		  gafftype = 25; // default to H on aliphatic carbon
                  jatm = atom.iat[i][0];
                  for (j=0; j < 3; j++)   // Acetylene
                  {
                      if (atom.bo[jatm][j] == 3)
                      {
                        mm3type = 124;
                        goto L_10;
                      }
                  }
                  newg = count_ewg(jatm);
		  if (atom.flags[jatm]  & aromatic_mask)
		    {
		      gafftype = 25;
		      if (newg == 1) gafftype = 56;  // h4
		      if (newg == 2) gafftype = 57; // h5
		    } else
		    {
		      if (newg == 1) gafftype = 53;  //  h1
		      if (newg == 2) gafftype = 54;  //  h2
		      if (newg == 3) gafftype = 55;  //  h3
		    }

                  goto L_10;
              } else if (atom.atomnum[atom.iat[i][0]] == 7)  // nitrogen
              {
                  jatm = atom.iat[i][0];
                  iatype = atom.mmx_type[jatm];
                  noxide = FALSE;
                  ismetal = FALSE;
                  mmxtype = 23;
                  mm3type = 23;
                  mmfftype = 23;
		  gafftype = 28;
                  jj_bo = atom.tbo[jatm];
                  for (j=0; j < atom.nconnect[jatm]; j++)
                  {
                         if (atom.mmx_type[atom.iat[jatm][j]] >= 300)
                           ismetal = TRUE;
                         if (atom.atomnum[atom.iat[jatm][j]] == 8)
                         {
                             katm = atom.iat[jatm][j];
                             jjk = atom.nconnect[katm];
                             if (jjk == 1)
                                noxide = TRUE;
                         }
                  }
                  if (jj_bo == 4)   //  N+
                  {
                      mmxtype = 24;
                      mm3type = 48;
                      if (ismetal == TRUE)
                      {
                          mmxtype = 23;
                          mm3type = 23;
                      }
                      mmfftype = 36;
                      if (noxide == TRUE)
                         mmfftype = 23;
                      goto L_10;
                  }
                  if (atom.flags[jatm] & aromatic_mask)
                  {
                      mmxtype = 23;
                      mm3type = 23;
                      mmfftype = 23;
                      if (is_cyclo5(jatm,array) )
                      {
                          jjk = 0;
                          for (j=0; j < 5; j++)
                          {
                              if (atom.atomnum[array[j]] == 7)
                                 jjk++;
                          }
                          if (jjk == 2)
                          {
                              jj_bo = 0;
                              for (j=0; j < 5; j++)
                              {
                                  if (atom.atomnum[array[j]] == 7 && array[j] != jatm)
                                  {
                                      katm = array[j];
				      jj_bo = atom.tbo[katm];
                                      if (jj_bo == 4)
                                      {
                                          mmfftype = 36;
                                          goto L_10;
                                      }
                                  }
                              }
                          }
                      }
                      goto L_10;
                  }
                  if (atom.mmff_type[atom.iat[i][0]] == 56 || atom.mmff_type[atom.iat[i][0]] == 55)
                  {
                       mmfftype = 36;
                       goto L_10;
                  }
                  if (atom.mmff_type[atom.iat[i][0]] == 62)
                  {
                      mmfftype = 23;
                      goto L_10;
                  }
                  for (j=0; j < atom.nconnect[jatm]; j++)
                  {
                      if (atom.iat[jatm][j] != 0 && atom.mmx_type[atom.iat[jatm][j]] != 20)
                      {
                         if (atom.bo[jatm][j] == 2 && (atom.atomnum[atom.iat[jatm][j]] == 6 || atom.atomnum[atom.iat[jatm][j]] == 7))  // imine
                         {
                             mmxtype = 23;
                             mm3type = 23;
                             mmfftype = 27;
                             goto L_10;
                         }
                         if (atom.bo[jatm][j] == 2 && atom.atomnum[atom.iat[jatm][j]] == 16 )  // imine
                         {
                             mmxtype = 23;
                             mm3type = 23;
                             mmfftype = 28;
                             goto L_10;
                         }
                         if (atom.atomnum[atom.iat[jatm][j]] == 16 && atom.mmff_type[atom.iat[jatm][j]] == 18) // thioamide
                         {
                             mmxtype = 23;
                             mm3type = 23;
                             mmfftype = 28;
                             goto L_10;
                         }
                         if (atom.atomnum[atom.iat[jatm][j]] == 6 || atom.atomnum[atom.iat[jatm][j]] == 7)  // amide and enamine
                         {
                             katm = atom.iat[jatm][j];
                             for (k=0; k < atom.nconnect[katm]; k++)
                             {
                                 if (atom.iat[katm][k] != 0 && atom.iat[katm][k] != jatm)
                                 {
                                     if (atom.bo[katm][k] == 3 && atom.atomnum[atom.iat[katm][k]] == 7)
                                     {
                                         mmxtype = 23;
                                         mm3type = 23;
                                         if (atom.mmff_type[jatm] != 8)
                                            mmfftype = 28;
                                         else
                                            mmfftype = 23;
                                         goto L_10;
                                     }
                                     if (atom.bo[katm][k] == 2)
                                     {
                                         if (atom.atomnum[atom.iat[katm][k]] == 8 || atom.atomnum[atom.iat[katm][k]] == 16)  // amide
                                         {
                                             mmxtype = 23;
                                             mm3type = 28;
                                             mmfftype = 28;
                                             goto L_10;
                                         } else if (atom.atomnum[atom.iat[katm][k]] == 6)
                                         {
                                             mmxtype = 23;
                                             mm3type = 28;
                                             if (atom.mmff_type[jatm] != 8)
                                               mmfftype = 28;
                                             else
                                               mmfftype = 23;
                                             goto L_10;
                                         } else if ( atom.atomnum[atom.iat[katm][k]] == 7)  // HN-c=n and HN-c=n+
                                         {
                                           mmfftype = 28;
                                           noxide = FALSE;
                                           latm = atom.iat[katm][k];
					   jjk = atom.tbo[latm];
                                           for (l=0; l < atom.nconnect[latm]; l++)
                                             {
                                                 if (atom.atomnum[atom.iat[latm][l]] == 8)
                                                    noxide = TRUE;
                                             }
                                           if (jjk == 4 && noxide == FALSE && atom.mmff_type[jatm] != 40)
                                             mmfftype = 36;
                                           goto L_10;
                                         }
                                     }
                                 }
                             }
                         }
                      }
                  }
                  mmxtype = 23;  // amine 
                  mm3type = 23;
                  mmfftype = 23;
                  goto L_10;
              } else if (atom.atomnum[atom.iat[i][0]] == 8) // oxygen
              {
                  mmxtype = 21;
                  mm3type = 21;
                  mmfftype = 21;
		  gafftype = 27;
                  jatm = atom.iat[i][0];
                  if (atom.iat[jatm][0] != i)
                    katm = atom.iat[jatm][0];
                  else
                    katm = atom.iat[jatm][1];
                  jjk = 0;  // number of hydrogens attached to jatm
                  for (j=0; j < atom.nconnect[jatm]; j++)
                  {
                      if (atom.iat[jatm][j] != 0 && atom.atomnum[atom.iat[jatm][j]] == 1)
                        jjk++;
                  }
                  if (atom.atomnum[katm] == 1 && jjk == 2) // water
                  {
                      mmxtype = 21;
                      mm3type = 21;
                      mmfftype = 31;
                      goto L_10;
                  } else if (atom.atomnum[katm] == 1 && jjk == 3) // h3o+
                  {
                      mmxtype = 24;
                      mm3type = 21;
                      mmfftype = 50;
                      goto L_10;
                  }
                  if (atom.atomnum[katm] == 15) // h-o-p
                  {
                      mmxtype = 24;
                      mm3type = 24;
                      mmfftype = 24;
                      goto L_10;
                  }
                  if (atom.atomnum[katm] == 16) // h-o-s
                  {
                      mmxtype = 24;
                      mm3type = 24;
                      mmfftype = 33;
                      goto L_10;
                  }
                  jjk = atom.nconnect[jatm];  // total attachments to oxygen
                  if (jjk == 3 || jjk == 4)  // O+
                  {
                      mmxtype = 24;
                      mm3type = 21;
                      mmfftype = 50;
                      goto L_10;
                  }                      
                  for (j=0; j < atom.nconnect[katm]; j++)
                  {
                      if (atom.atomnum[atom.iat[katm][j]] == 8 && atom.iat[katm][j] != jatm)
                      {
                          if (atom.bo[katm][j] == 2)   // carboxyl
                          {
                              mmxtype = 24;
                              mm3type = 24;
                              mmfftype = 24;
                              goto L_10;
                          }
                      }
                      if (atom.atomnum[katm] == 6 && (atom.atomnum[atom.iat[katm][j]] == 6 || atom.atomnum[atom.iat[katm][j]] == 7)
                          && atom.iat[katm][j] != jatm)  // enol
                      {
                          if (atom.bo[katm][j] == 2)
                          {
                              mmxtype = 28;
                              mm3type = 73;
                              mmfftype = 29;
                              goto L_10;
                          }
                      }
                  }
                  if (jjk == 2)  // OH
                  {
                      mmxtype = 21;
                      mm3type = 21;
                      mmfftype = 21;
                      for (j=0; j < atom.nconnect[jatm]; j++)  // H-O=C
                      {
                          if (atom.iat[jatm][j] == katm && atom.bo[jatm][j] == 2)
                          {
                              mmxtype = 24;
                              mmfftype = 52;
                              goto L_10;
                          }
                      }
                      goto L_10;
                  }
                  goto L_10;
              } else if (atom.atomnum[atom.iat[i][0]] == 5) // boron
              {
                  mmxtype = 23;
                  mm3type = 5;
                  mmfftype = 71;
                  goto L_10;
              } else if (atom.atomnum[atom.iat[i][0]] == 9 || atom.atomnum[atom.iat[i][0]] == 17 ||
                          atom.atomnum[atom.iat[i][0]] == 35 || atom.atomnum[atom.iat[i][0]] == 53) // halogens
              {
                  mmxtype = 21;
                  mm3type = 5;
                  mmfftype = 71;
                  goto L_10;
              } else if (atom.atomnum[atom.iat[i][0]] == 16) // sulfur
              {
                  mmxtype = 21;
                  mm3type = 5;
                  mmfftype = 71;
		  gafftype = 30;
                  goto L_10;
              } else if (atom.atomnum[atom.iat[i][0]] == 15) // phosphorous
              {
                  mmxtype = 23;
                  mm3type = 5;
                  mmfftype = 71;
		  gafftype = 29;
                  goto L_10;
              } else  //  bridging hydrogens
              {
                  jji = atom.nconnect[i];
                  mm3type = 5;
                  if (jji == 2)
                    mmxtype = 70;
                  else
                  {
                    mmxtype = 5;
                    mm3type = 5;
                    mmfftype = 5;
                  }
                  goto L_10;
              }
          } else if (atom.atomnum[i] == 5) // boron
          {
              if (atom.mmx_type[i] == 27)  // four coordinate boron 
              {
                 mmxtype = 27;
                 mm3type = 27;
                 mmfftype = 0;
                 goto L_10;
              }
              if (atom.mmx_type[i] == 43)  // Transition state boron 
              {
                 mmxtype = 43;
                 mm3type = 0;
                 mmfftype = 0;
                 goto L_10;
              }
              jji = atom.nconnect[i];
              if (jji == 4)
              {
                 mmxtype = 27;
                 mm3type = 27;
                 mmfftype = 0;
              } else
              {
                 mmxtype = 26;
                 mm3type = 26;
                 mmfftype = 0;
              }
              goto L_10;
// =========================== Carbon ===============================
          } else if (atom.atomnum[i] == 6) // carbon
          {
              ndouble = 0;
              ntriple = 0;
              jji = atom.nconnect[i];
              icycl3 = find_rsize(3,i);
              icycl4 = find_rsize(4,i);
              icycl5 = find_rsize(5,i);
              icycl6 = find_rsize(6,i);
              
              for (j=0; j < jji; j++)
              {
                  if (atom.bo[i][j] == 3)
                     ntriple++;
                  if (atom.bo[i][j] == 2)
                     ndouble++;
              }
              //  check here for types that can not be done by rules but must be set by user
              //  and thus should not be changed
              if (mmxtype == 29)   // C radical
              {
                  mmxtype = 29;
                  mm3type = 29;
                  mmfftype = 0;
                  goto L_10;
              } else if (mmxtype == 30) // C cation
              {
                  mmxtype = 30;
                  mm3type = 30;
                  mmfftype = 0;
                  goto L_10;
              } else if (mmxtype == 48) // C anion
              {
                  mmxtype = 48;
                  mm3type = 0;
                  mmfftype = 0;
                  goto L_10;
              } else if (mmxtype == 40) // aromatic carbon
              {
                  mmxtype = 40;
                  mm3type = 50;
                  mmfftype = 37;
                  goto L_10;
              } else if (mmxtype == 49 || mmxtype == 50 || mmxtype == 51 || mmxtype == 52)  // TS atoms
              {
                  mm3type = 0;
                  mmfftype = 0;
                 goto L_10;
              }
//   check rings
              if (icycl5 >= 1 && icycl6 >= 1 && atom.flags[i] & aromatic_mask)  // CB and CN
              {
                  get_rsize(i,5,0,array);
                  for (j=0; j < 5; j++)
                  {
                      if ( !(atom.flags[array[j]] & aromatic_mask))
                         goto L_NOTPURINE;
                  }
                  nnh = 0;
                  nnit = 0;
                  for (j=0; j < 5; j++)
                  {
                      if (atom.atomnum[array[j]] == 7)
                      {
                          nnit++;
                          for (k=0; k < 4; k++)
                          {
                              if (atom.atomnum[atom.iat[array[j]][k]] == 1)
                                nnh++;
                          }
                      }
                  }
                  if (nnit == 2 && nnh == 2)
                  {
                      mmxtype = 2;
                      mm3type = 2;
                      mmfftype = 78;
                      goto L_10;
                  }
//
                  nnh = 0;
                  nnit = 0;
                  for (j=0; j < atom.nconnect[i]; j++)
                  {
                   if (atom.atomnum[atom.iat[i][j]] == 7)
                   {
                       nnit++;
                       for (k=0; k < atom.nconnect[atom.iat[i][j]]; k++)
                       {
                           if (atom.atomnum[atom.iat[atom.iat[i][j]][k]] == 1)
                           {
                               nnh++;
                           }
                       }
                   }
                  }
              }
L_NOTPURINE:
              if (icycl5 >= 1 && atom.flags[i] & aromatic_mask)
              {
                  mmxtype = 2;
                  mm3type = 2;
                  for (ij=0; ij < icycl5; ij++)
                  {
                      // do mmff types
                      get_rsize(i,5,ij,array);
		      if (aromatic_5(array,atom.flags,atom.atomnum,atom.nconnect,atom.input_charge,atom.iat,atom.bo))
                      {
                          icurr = -1;
                          nplus = FALSE;
                          for (j=0; j < 5; j++)
                          {
                              if (array[j] == i)
                                icurr = j;
                              if (atom.atomnum[array[j]] == 7)
                              {
                                  jjk = atom.nconnect[array[j]];
                                  jj_bo = atom.tbo[array[j]];
                                  if (jjk == 2 && jj_bo == 2) // divalent N anion
                                  {
                                      mmfftype = 78;
                                      goto L_10;
                                  }
                                  if (jj_bo == 4)
                                     nplus = TRUE;
                              }
                          }

                          // check alpha
                          ia = (icurr+4)%5;
                          ib = (icurr+6)%5;
                          if (atom.atomnum[array[ia]] == 7 && atom.atomnum[array[ib]] == 7)  // n=c-n
                          {
                              jatm = array[ia];
                              katm = array[ib];
                              jbo = atom.tbo[jatm];
                              kbo = atom.tbo[katm];
                              if ( jbo == 4 && kbo == 3)
                              {
                                  adjn = FALSE;
                                  for (j=0; j < atom.nconnect[katm]; j++)
                                  {
                                      if (atom.iat[katm][j] != 0 && atom.bo[katm][j] == 2)
                                         adjn = TRUE;
                                  }
                                  if (adjn == FALSE)
                                  {
                                     mmfftype = 80;
                                     goto L_10;
                                  }
                              }else if (jbo == 3 && kbo == 4)
                              {
                                  adjn = FALSE;
                                  for (j=0; j < atom.nconnect[jatm]; j++)
                                  {
                                      if (atom.iat[jatm][j] != 0 && atom.bo[jatm][j] == 2)
                                         adjn = TRUE;
                                  }
                                  if (adjn == FALSE)
                                  {
                                     mmfftype = 80;
                                     goto L_10;
                                  }
                              }
                          }
                          if (nplus == TRUE) // found N+ in ring
                          {
                              noxide = FALSE;
                              jjk = 0;
                              for (j=0; j < 5; j++)
                              {
                                  if (atom.atomnum[array[j]] == 7)
                                  {
                                      nit[jjk] = j;
                                      jjk++;
                                  }
                              }
                              if (jjk >= 2)
                              {
                                  jatm = array[nit[0]];
                                  katm = array[nit[1]];
                                  jbo = atom.tbo[jatm];
                                  kbo = atom.tbo[katm];
                                  adjn = FALSE;
                                  adjn1 = FALSE;
                                  for (k=0; k < atom.nconnect[jatm]; k++)
                                  {
                                       if (atom.iat[jatm][k] != 0 && atom.bo[jatm][k] == 2)
                                              adjn = TRUE;
				  }
                                  for (k=0; k < atom.nconnect[katm]; k++)
                                  {
                                       if (atom.iat[katm][k] != 0 && atom.bo[katm][k] == 2)
                                              adjn1 = TRUE;
                                  }
                                  if ( (jbo == 4 && kbo == 3 && adjn1 == FALSE) ||
                                       (jbo == 3 && kbo == 4 && adjn == FALSE)  )
                                  {
                                      noxide = FALSE;
                                      if (jbo == 4)
                                        ia1 = jatm;
                                      else
                                        ia1 = katm;
                                      for (k=0; k < atom.nconnect[ia1]; k++)
                                      {
                                          if (atom.iat[ia1][k] != 0 && atom.atomnum[atom.iat[ia1][k]] == 8)
                                          {
                                              ib1 = atom.iat[ia1][k];
                                              jjk = atom.nconnect[ib1];
                                              if (jjk == 1)
                                                 noxide = TRUE;
                                          }
                                      }
                                      if (noxide == FALSE)
                                      {             
                                         mmfftype = 78;
                                         goto L_10;
                                      }
                                  }
                              }
                          }

                          if (atom.atomnum[array[ia]] != 6 && atom.atomnum[array[ib]] != 6)  // x=c-x
                          {
                              jatm = array[ia];
                              katm = array[ib];
                              jbo = atom.tbo[jatm];
                              kbo = atom.tbo[katm];
                              if ( (jbo == 4 && kbo == 2) || (jbo == 2 && kbo == 4))
                              {
                                  for (j=0; j < atom.nconnect[i]; j++)
                                  {
                                      if (atom.iat[i][j] != 0 && atom.iat[i][j] != jatm && atom.iat[i][j] != katm)
                                      {
                                          if (atom.atomnum[atom.iat[i][j]] == 7)
                                          {
                                               mmfftype = 80;
                                               goto L_10;
                                          }
                                      }
                                  }
                              }
                          }
                          if (atom.atomnum[array[ia]] == 7)  // alpha n
                          {
                              jatm = array[ia];
                              adjn = FALSE;
                              for (j=0; j < atom.nconnect[jatm]; j++)
                              {
                                  if (atom.iat[jatm][j] != 0 && atom.bo[jatm][j] ==2)
                                     adjn = TRUE;
                              }
                              if (adjn == FALSE)
                              {
                                  mmfftype = 63;
                                  goto L_10;
                              }
                          }
                          if (atom.atomnum[array[ib]] == 7)  // alpha n
                          {
                              jatm = array[ib];
                              adjn = FALSE;
                              for (j=0; j < atom.nconnect[jatm]; j++)
                              {
                                  if (atom.iat[jatm][j] != 0 && atom.bo[jatm][j] == 2)
                                     adjn = TRUE;
                              }
                              if (adjn == FALSE)
                              {
                                  mmfftype = 63;
                                  goto L_10;
                              }
                          }
                          if (atom.atomnum[array[ia]] == 8 || atom.atomnum[array[ib]] == 8 ||  // alpha o or s
                              atom.atomnum[array[ia]] == 16 || atom.atomnum[array[ib]] == 16 )
                          {
                              mmfftype = 63;
                              goto L_10;
                          }
                          // check beta
                          ia = (icurr+3)%5;
                          ib = (icurr+7)%5;
                          if (atom.atomnum[array[ia]] == 7)  // c=x-n
                          {
                              jatm = array[ia];
                              adjn = FALSE;
                              for (j=0; j < atom.nconnect[jatm]; j++)
                              {
                                  if (atom.iat[jatm][j] != 0 && atom.bo[jatm][j] == 2)
                                     adjn = TRUE;
                              }
                              if (adjn == FALSE)
                              {
                                  mmfftype = 64;
                                  goto L_10;
                              }
                          }
                          if (atom.atomnum[array[ib]] == 7)  // c=x-n
                          {
                              jatm = array[ib];
                              adjn = FALSE;
                              for (j=0; j < atom.nconnect[jatm]; j++)
                              {
                                  if (atom.iat[jatm][j] != 0 && atom.bo[jatm][j] == 2)
                                     adjn = TRUE;
                              }
                              if (adjn == FALSE)
                              {
                                  mmfftype = 64;
                                  goto L_10;
                              }
                          }
                          if (atom.atomnum[array[ia]] == 8 || atom.atomnum[array[ib]] == 8 ||  // beta o or s
                              atom.atomnum[array[ia]] == 16 || atom.atomnum[array[ib]] == 16 )
                          {
                              mmfftype = 64;
                              goto L_10;
                          }
                      }
                  }
              }
              if (icycl6 >= 1 && atom.flags[i] & aromatic_mask)
              {
                  mmxtype = 2;
                  mm3type = 2;
                  mmfftype = 37;
                  for (k = 0; k < icycl6; k++)
                  {
                      get_rsize(i,6,k,array);
                      {
                          full_ring = TRUE;
                          non_pi = FALSE;
                          for (j=0; j < 6; j++)  // check all ring atoms are pi
                          {
                              if ( !(atom.flags[array[j]] & aromatic_mask) )
                                 non_pi = TRUE;
                          }
                          if (non_pi == TRUE)
                             full_ring = FALSE;
                      }
                      if (full_ring == TRUE)
                         break;
                  }
                  if (full_ring == FALSE) mmfftype = 2;
                  goto L_10;
              }
              if (icycl3 >= 1 )
              {
                  mmxtype = 22;
                  mm3type = 22;  // type 38 for cyclopropene, type 67 for cyclopropanone
                  mmfftype = 22;
		  if (ndouble >= 1) mmfftype = 2;
                  if (ndouble >= 1)
                  {
                      for (j=0; j < atom.nconnect[i]; j++)
                      {
                          if (atom.iat[i][j] != 0)
                          {
                              if (atom.bo[i][j] == 2 && atom.atomnum[atom.iat[i][j]] == 8)
                              {
                                  mmxtype = 3;
                                  mm3type = 67;
				  mmfftype = 3;
                                  goto L_10;
                              }
                          }
                      }
                     mm3type = 38;
                  }
                  goto L_10;
              }
//   start of check on bondorders
              if (ntriple == 1)
              {
                  mmxtype = 4;
                  mm3type = 4;
                  mmfftype = 4;
                  if (jji == 1 && atom.atomnum[atom.iat[i][0]] == 7)
                  {
                     mmfftype = 60;  // isonitrile carbon
                  }
                  if ( (atom.mmx_type[atom.iat[i][0]] >= 300) &&
                       atom.mmx_type[atom.iat[i][1]] == 46)
                  {
                         mmxtype = 63;
                         goto L_10;
                  }
                  if ( (atom.mmx_type[atom.iat[i][1]] >= 300) &&
                       atom.mmx_type[atom.iat[i][0]] == 46)
                  {
                         mmxtype = 63;
                         goto L_10;
                  }
                  for (j=0; j < atom.nconnect[i]; j++)
                  {
                      if ( (atom.mmx_type[atom.iat[i][j]] >= 300) && atom.bo[i][j] == 3)
                      {
                         mmxtype = 62;
                         goto L_10;
                      }
                  }
                  goto L_10;
              }
              if (ndouble == 2)   // allenes and ketenes
              {
                  mmxtype = 4;
                  mm3type = 68;
                  mmfftype = 4;
                  for (j=0; j < atom.nconnect[i]; j++)  // metal carbene
                  {
                      if ( atom.mmx_type[atom.iat[i][j]] >= 300)
                         mmxtype = 61;
                  }
                  goto L_10;
              }
              if (ndouble == 1)
              {
                  for (j=0; j < atom.nconnect[i]; j++)
                  {
                      if (atom.iat[i][j] != 0)
                      {
                          if (atom.bo[i][j] == 2)
                          {
                              jatm = atom.iat[i][j];
                              break;
                          }
                      }
                  }
                  if (atom.atomnum[jatm] == 15) // c=p
                  {
                      mmxtype = 2;
                      mm3type = 2;
                      mmfftype = 3;
                      goto L_10;
                  }
                  if (atom.atomnum[jatm] == 16) // c=s
                  {
                      mmxtype = 2;
                      mm3type = 2;
                      mmfftype = 3;
                      for (k=0; k < atom.nconnect[i]; k++)
                      {
                          if (atom.iat[i][k] != 0 && atom.bo[i][k] != 9)
                          {
                              if (atom.atomnum[atom.iat[i][k]] == 16 && atom.iat[i][k] != jatm)
                              {
                                  katm = atom.iat[i][k];
                                  jjk = atom.nconnect[katm];
                                  if (jjk == 1)  // thiocarboxylate
                                  {
                                      mmfftype = 41;
                                      goto L_10;
                                  }
                              }
                          }
                      }    
                      goto L_10;
                  }
                  if (atom.atomnum[jatm] == 8)  // C=O
                  {
                      mmxtype = 3;
                      mm3type = 3;
                      mmfftype = 3;
                      if (is_ring31(i))          // cyclopropanone
                          mm3type = 67;
                      else if (is_ring41(i))   // cyclobutanone
                          mm3type = 58;
                      else                     // carboxylate
                      {
                          for(k=0; k < atom.nconnect[i]; k++)
                          {
                              if (atom.iat[i][k] != 0)
                              {
                                  if (atom.atomnum[atom.iat[i][k]] == 8 && atom.iat[i][k] != jatm)
                                  {
                                      katm = atom.iat[i][k]; 
                                      jjk = atom.nconnect[katm];
                                      if (jjk == 1)
                                      {
                                           mmfftype = 41;
                                           break;
                                      }
                                  }
                              }
                          }
                      }
                      goto L_10;
                  } else if (atom.atomnum[jatm] == 7) // C=N
                  {
                      mmxtype = 2;
                      mm3type = 2;
                      mmfftype = 3;
                      if (jji == 3)
                      {
                          if (atom.iat[i][0] == jatm)
                          {
                              katm = atom.iat[i][1];
                              latm = atom.iat[i][2];
                          } else if (atom.iat[i][1] == jatm)
                          {
                              katm = atom.iat[i][0];
                              latm = atom.iat[i][2];
                          } else if (atom.iat[i][2] == jatm)
                          {
                              katm = atom.iat[i][0];
                              latm = atom.iat[i][1];
                          }
                          if (atom.atomnum[jatm] == 7 && atom.atomnum[katm] == 7 && atom.atomnum[latm] == 7)
                          {
                              jbo = atom.tbo[jatm];
                              kbo = atom.tbo[katm];
                              lbo = atom.tbo[latm];
                              adjn = FALSE;
                              for (j=0; j < atom.nconnect[katm]; j++)
                              {
                                       if (atom.bo[katm][j] == 2)
                                           adjn = TRUE;                                 
                                       if (atom.bo[latm][j] == 2)
                                           adjn = TRUE;                                 
                              }
                              if (jbo == 4 && kbo == 3 && lbo == 3 && adjn == FALSE)
                              {  
                                  mmfftype = 57;
                                  goto L_10;
                              }
                          }
                          if (atom.atomnum[jatm] == 7 && atom.atomnum[katm] == 7)
                          {
                              jbo = atom.tbo[jatm];
                              kbo = atom.tbo[katm];
                              adjn = FALSE;
                              for (j=0; j < atom.nconnect[katm]; j++)
                              {
                                       if (atom.bo[katm][j] == 2)
                                           adjn = TRUE;
                              }
                              if (jbo == 4 && kbo == 3 && adjn == FALSE)
                              {  
                                  mmfftype = 57;
                                  goto L_10;
                              }
                          }
                          if (atom.atomnum[jatm] == 7 && atom.atomnum[latm] == 7)
                          {
                              jbo = atom.tbo[jatm];
                              lbo = atom.tbo[latm];
                              adjn = FALSE;
                              for (j=0; j < atom.nconnect[latm]; j++)
                              {
                                       if (atom.bo[latm][j] == 2)
                                           adjn = TRUE;
                              }
                              if (jbo == 4 && lbo == 3 && adjn == FALSE)
                              {  
                                  mmfftype = 57;
                                  goto L_10;
                              }
                          }
                      }
                      goto L_10;
                  } else if (atom.atomnum[jatm] == 6) // C=C
                  {
                      mmxtype = 2;
                      mm3type = 2;
                      mmfftype = 2;
                      if (is_ring31(i))    // cyclopropene
                      {
                          mm3type = 38;
                          goto L_10;
                      }
                      if (is_ring41(i))  // cyclobutene
                      {
                          mmxtype = 57;
                          mm3type = 57;
                          mmfftype = 30;
                          goto L_10;
                      }

                      goto L_10;
                  } else   // default for cases not dealt with yet
                  {
                      mmxtype = 2;
                      mm3type = 2;
                      mmfftype = 2;
                      goto L_10;
                  }
                  goto L_10;
              }
      // get here with only single bonds to carbon
              mmxtype = 1;
              mm3type = 1;
              mmfftype = 1;
              if (is_ring31(i))
              {
                  mmxtype = 22;
                  mm3type = 22;
                  mmfftype = 22;
                  goto L_10;
              } else if (is_ring41(i))
              {
                  mmxtype = 56;
                  mm3type = 56;
                  mmfftype = 20;
                  goto L_10;
              } 
              goto L_10;
// ============   Nitrogen atom types =================
          } else if (atom.atomnum[i] == 7) // nitrogen
          {
              jji = atom.nconnect[i];
	      jjbo = atom.tbo[i];
              nc = 0;
              nh = 0;
	      nox = 0;
	      nox_double = 0;
	      nsulf = 0;
              ndouble = 0;
              ntriple = 0;
              nplus = FALSE;
              noxide = FALSE;
              for (j=0; j < jji; j++)
              {
                if (atom.atomnum[atom.iat[i][j]] == 1)
                  nh++;
                else if (atom.atomnum[atom.iat[i][j]] == 6)
                  nc++;
                else if (atom.atomnum[atom.iat[i][j]] == 7)
                  nnit++;
                else if (atom.atomnum[atom.iat[i][j]] == 8)
		  {
		    nox++;
		    if (atom.nconnect[atom.iat[i][j]] == 1)  // noxide
		      noxide = TRUE;
		  }
                else if (atom.atomnum[atom.iat[i][j]] == 16)
                  nsulf++;
		if (atom.bo[i][j] == 2)
		  ndouble++;
		if (atom.bo[i][j] == 3)
		  ntriple++;
              }
	      // triple bond atom types
	      if (ntriple  == 1)
		{
		  if (jji == 2)
		    {
		      mmfftype = 61;  // isonitrile
		      mmxtype = 68;
		    }
		  else 
		    {
		      mmxtype = 10;
		      mmfftype = 42;  // nitrile
		    }
		  goto L_10;
		}
        // two double bonds to N - azides
              if (ndouble == 2 && jji == 2)
              {
                  mmxtype = 37;
                  mmfftype = 53;
                  goto L_10;
              }
	      //  terminal azide and diazo and NSO
	      if (jji == 1 && ndouble == 1)
		{
		  mmxtype = 37;
		  mmfftype = 47;
		  goto L_10;
		}
	      // R-N=SO-R
	      if (jji == 2 && ndouble == 1 && nsulf == 1)
		{
		  for (j=0; j < jji; j++)
		    {
		      if (atom.atomnum[atom.iat[i][j]] == 16 && atom.bo[i][j] == 2)
			{
			  mmfftype = 48;
			  goto L_10;
			}
		    }
		}
	      // anionic divalent nitrogen - not aromatic and not R5
	      if (jji == 2 && jjbo == 2 && atom.input_charge[i] == -1 && !(atom.flags[i] & aromatic_mask) )
		{
		  mmfftype = 62;
		  goto L_10;
		}
	      // ammonium  N+ - sp3
	      if (jji == 4)
		{
		  mmxtype = 41;
		  if (noxide == FALSE)
		    mmfftype = 34;
		  else 
		    mmfftype = 68;
		  goto L_10;
		}
	      //
              icycl5 = find_rsize(5,i);
              icycl6 = find_rsize(6,i);
	      //
	      // looking for types 65, 66, 79, 39, 81,82
              if (icycl5 >= 1 && atom.flags[i] & aromatic_mask)
              {
		if (jji == 3 && jjbo == 3)  // default amine type in R5 & aromatic
		  mmfftype = 39;
		else if (jji == 2 && jjbo == 3) // default imine type in R5 & aromatic
		  mmfftype = 79;

		if (jji == 2 && jjbo == 2)  // divalent anion in R5 and aromatic
		  {
		    mmfftype = 76;
		    goto L_10;
		  }
                  if (jjbo == 4 && noxide == TRUE)
                  {
                      mmxtype = 41;
                      mmfftype = 82;
                      goto L_10;
                  }
                  if (jjbo == 4)
                  {
                      mmxtype = 41;
                      mmfftype = 81;
                      goto L_10;
                  }
                  // assume only one ring need to fix this
                  get_rsize(i,5,0,array);
                  if (ndouble == 1)  // check for alpha or beta N,O,S
                  {
                      mmxtype = 37;
                      icurr  = -1;
                      for (j=0; j < 5; j++)
                      {
                          if (array[j] == i)
                             icurr = j;
                      }
                      if (icurr == 0)
                      {
                          ia = array[1];
                          ib = array[4];
                      } else if (icurr == 4)
                      {
                          ia = array[0];
                          ib = array[3];
                      } else
                      {
                          ia = array[icurr-1];
                          ib = array[icurr+1];
                      }
                      
                      if ( (atom.atomnum[ia] == 7 && atom.nconnect[ia] == 3 && atom.tbo[ia] == 3) || atom.atomnum[ia] == 8 || atom.atomnum[ia] == 16)  // alpha
                      {
                             mmfftype = 65;
                             goto L_10;
                      }
                      if ( (atom.atomnum[ib] == 7 && atom.nconnect[ib] == 3 && atom.tbo[ib] == 3) || atom.atomnum[ib] == 8 || atom.atomnum[ib] == 16)
                      {
                              mmfftype = 65;
                              goto L_10;
                      }
                      // now check beta
                      if (icurr == 0)
                      {
                          ia = array[2]; ib = array[3];
                      } else if (icurr == 1)
                      {
                          ia = array[3]; ib = array[4];
                      } else if (icurr == 2)
                      {
                          ia = array[4]; ib = array[0];
                      } else if (icurr == 3)
                      {
                          ia = array[0]; ib = array[1];
                      } else if (icurr == 4)
                      {
                          ia = array[1]; ib = array[2];
                      }
                      if ((atom.atomnum[ia] == 7 && atom.nconnect[ia] == 3 && atom.tbo[ia] == 3) || atom.atomnum[ia] == 8 || atom.atomnum[ia] == 16)
                      {
                              mmfftype = 66;
                              goto L_10;                           
                      }
                      if ((atom.atomnum[ib] == 7 && atom.nconnect[ib] == 3 && atom.tbo[ib] == 3) || atom.atomnum[ib] == 8 || atom.atomnum[ib] == 16)
                      {
                              mmfftype = 66;
                              goto L_10;                            
                      }
                  }
		  //
                  // single bonds only to nitrogen in R5 and aromatic
                  icurr = -1;
                  nplus = FALSE;
                  for (j=0; j < 5; j++)
                  {
                      if (array[j] == i)
                         icurr = j;
                      if (atom.atomnum[array[j]] == 7)
                         nplus = TRUE;
                  }
                  if (nplus == TRUE) // looking for N-c=N+
                  {
                      if (icurr == 0)
                      {
                          ia = array[2]; ib = array[3];
                      } else if (icurr == 1)
                      {
                          ia = array[3]; ib = array[4];
                      } else if (icurr == 2)
                      {
                          ia = array[4]; ib = array[0];
                      } else if (icurr == 3)
                      {
                          ia = array[0]; ib = array[1];
                      } else if (icurr == 4)
                      {
                          ia = array[1]; ib = array[2];
                      }
                      if (atom.atomnum[ia] == 7)
                      {
                          jj_bo = atom.tbo[ia];
                          noxide = FALSE;
                          for (k=0; k < atom.nconnect[ia]; k++)
                          {
                              if (atom.atomnum[atom.iat[ia][k]] == 8)
                              {
                                  katm = atom.iat[ia][k];
                                  jjk = atom.nconnect[katm];
                                  if (jjk == 1)
                                     noxide = TRUE;
                              }
                          }
                          if (jj_bo == 4 && noxide == FALSE)
                          {
                              mmfftype = 81;  // histidine
                              goto L_10;
                          }
                      }
                      if (atom.atomnum[ib] == 7)
                      {
                          jj_bo = atom.tbo[ib];
                          noxide = FALSE;
                          for (k=0; k < atom.nconnect[ib]; k++)
                          {
                              if (atom.atomnum[atom.iat[ib][k]] == 8)
                              {
                                  katm = atom.iat[ib][k];
                                  jjk = atom.nconnect[katm];
                                  if (jjk == 1)
                                     noxide = TRUE;
                              }
                          }
                          if (jj_bo == 4 && noxide == FALSE)
                          {
                              mmfftype = 81;  // histidine
                              goto L_10;
                          }
                      }
                  }
                  // failed tests use general types
                  goto L_10;
              }
	      // R6 and aromatic = types 38, 58, 69
              if (icycl6 >= 1 && atom.flags[i] & aromatic_mask)
              {
                  mmxtype = 37;
                  mmfftype = 38;
                  if (jjbo == 4 && noxide == TRUE)
                  {
                      mmxtype = 41;
                      mmfftype = 69;
                  } else if (jjbo == 4)
                  {
                      mmxtype = 41;
                      mmfftype = 58;
                  }
                  goto L_10;
              }
// non cyclic systems
          // single double bond to N
              if (ndouble == 1)
              {
                  mmxtype = 37;
                  mmfftype = 9;
                  if (jji == 3 && jjbo == 4) //=N+
                     mmxtype = 41;
                  
                  for (j=0; j < jji; j++)
                  {
                      if (atom.iat[i][j] != 0 && atom.bo[i][j] == 2)
                      {
                          jatm = atom.iat[i][j];
                          break;
                      }
                  }
                  if (atom.atomnum[jatm] == 8)  // N=O
                  {
                      if (jjbo == 4 && noxide == TRUE)   // nitro
                      {
                          mmxtype = 41;
                          mm3type = 46;
                          mmfftype = 45;
                          goto L_10;
                      } else if (jjbo == 3)  // nitroso
                      {
                          mmxtype = 37;
                          mm3type = 0;
                          mmfftype = 46;
                          goto L_10;
                      }
                  } else if (atom.atomnum[jatm] == 7) // N=N
                  {
                      mmfftype = 9;
		      if (jjbo == 4 && noxide == TRUE) // o-n=n
			mmfftype = 67;
                      goto L_10;                     
                  } else if (atom.atomnum[jatm] == 6) // N=C
                  {
		      mmfftype = 9;
                      jjk = atom.nconnect[jatm];
                      if (jjk == 3)
                      {
                          if (atom.iat[jatm][0] == i)
                          {
                              katm = atom.iat[jatm][1];
                              latm = atom.iat[jatm][2];
                          } else if (atom.iat[jatm][1] == i)
                          {
                              katm = atom.iat[jatm][0];
                              latm = atom.iat[jatm][2];
                          } else if (atom.iat[jatm][2] == i)
                          {
                              katm = atom.iat[jatm][0];
                              latm = atom.iat[jatm][1];
                          }
                          jbo = atom.tbo[jatm];
                          kbo = atom.tbo[katm];
                          lbo = atom.tbo[latm];
                          adjn = FALSE;
                          if (atom.atomnum[katm] == 7 && atom.atomnum[latm] == 7)
                          {
                             for (k=0; k < atom.nconnect[katm]; k++)
                             {
                                    if (atom.bo[katm][k] == 2)
			               adjn = TRUE;
			     }
                             for (k=0; k < atom.nconnect[latm]; k++)
			       {
                                    if (atom.bo[latm][k] == 2)
				      adjn = TRUE;
                             }
                             if (atom.tbo[i] == 4 && jbo == 4 && (kbo == 3 && lbo == 3) && adjn == FALSE)  ///  N+=C(N)N
                             {
                                mmfftype = 56;
                                goto L_10;
                             } else if (atom.tbo[i] == 3 && (kbo == 3 && lbo == 3) && adjn == FALSE)  ///  N=C(N)N
			       {
				 mmfftype = 9;
				 goto L_10;
			       }
                          }
                          if (atom.atomnum[katm] == 7)
                          {
                              jbo = atom.tbo[i];
                              kbo = atom.tbo[katm];
                              adjn = FALSE;
                              for (k=0; k < atom.nconnect[katm]; k++)
                              {
                                    if (atom.bo[katm][k] == 2)
                                       adjn = TRUE;
                              }
                              if (jbo == 4 && kbo == 3 && adjn == FALSE)
                              {
                                  mmfftype = 55;
                                  goto L_10;
                              }
                          }
                          if (atom.atomnum[latm] == 7)
                          {
                              jbo = atom.tbo[i];
                              lbo = atom.tbo[latm];
                              adjn = FALSE;
                              for (k=0; k < atom.nconnect[latm]; k++)
                              {
                                    if (atom.bo[latm][k] == 2)
                                       adjn = TRUE;
                              }
                              if (jbo == 4 && lbo == 3 && adjn == FALSE)
                              {
                                  mmfftype = 55;
                                  goto L_10;
                              }
                          }         
		      }
                      if (jjbo == 4) // c=n+
                      {
                          mmxtype = 41;
                          for (j=0; j < atom.nconnect[i]; j++)
                          {
                                  if (atom.atomnum[atom.iat[i][j]] == 8) // c=n+-o
                                  {
                                      katm = atom.iat[i][j];
                                      jjk = atom.nconnect[katm];
                                      if (jjk == 1)
                                        mmfftype = 67;
                                      else if (jjk == 2)
                                        mmfftype = 54;
                                      goto L_10;
                                  }
                          }
                          mmfftype = 54;
                          goto L_10;
                      }
                      mmxtype = 37;
                      mmfftype = 9;
                      goto L_10;                 
	          } // got here with not assignment check other side for mmff
		  // goto L_10;
	      }
        // goto here with only single bonds to nitrogen
               mmxtype = 8;
               mmfftype = 8;
	       // look for adjacent double bonds and triple bonds
               ndouble = 0;
	       ntriple = 0;
               for (j=0; j < jji; j++)
               {
                  katm = atom.iat[i][j];
                  for (k=0; k < atom.nconnect[katm]; k++)
                  {
                     if (atom.bo[katm][k] == 2)
                       ndouble++;
                     if (atom.bo[katm][k] == 3)
                        ntriple++;
		  }
               }
	       if (ntriple == 1)
		 {
		   mmfftype = 43;  // N-c%N
		   goto L_10;
		 }
               if (ndouble >= 1)
               {
                   mmxtype = 9;
                   mmfftype = 8;  // looking for N-x=y or N-x=N+
                   for (j=0; j < jji; j++)
                   {
		     if (atom.atomnum[atom.iat[i][j]] == 6)  // adjacent carbon
                       {
                           jatm = atom.iat[i][j];
                           jjk = atom.nconnect[jatm];
                           if (jjk == 3) // sp2 carbon
                           {
                               if (atom.iat[jatm][0] == i)
                               {
                                   if (atom.bo[jatm][1] == 2)
                                   {
                                      katm = atom.iat[jatm][1];
                                      latm = atom.iat[jatm][2];
                                   }else
                                   {
                                      latm = atom.iat[jatm][1];
                                      katm = atom.iat[jatm][2];
                                   }
                               } else if (atom.iat[jatm][1] == i)
                               {
                                   if (atom.bo[jatm][0] == 2)
                                   {
                                      katm = atom.iat[jatm][0];
                                      latm = atom.iat[jatm][2];
                                   }else
                                   {
                                      latm = atom.iat[jatm][0];
                                      katm = atom.iat[jatm][2];
                                   }
                               }else if (atom.iat[jatm][2] == i)
                               {
                                   if (atom.bo[jatm][0] == 2)
                                   {
                                      katm = atom.iat[jatm][0];
                                      latm = atom.iat[jatm][1];
                                   }else
                                   {
                                      latm = atom.iat[jatm][0];
                                      katm = atom.iat[jatm][1];
                                   }
                               }
                               if (atom.atomnum[katm] == 7 && atom.atomnum[latm] == 7) // N-c(N)=N
                               {
                                   jbo = atom.tbo[katm];
                                   kbo = atom.tbo[latm];
                                   adjn = FALSE;
                                   for(k=0; k < atom.nconnect[latm]; k++)
                                   {
                                          if (atom.bo[latm][k] == 2)
                                             adjn = TRUE;
                                   }
                                   if ( jbo == 4 && kbo == 3 && adjn == FALSE)
                                   {
                                       mmfftype =56;
                                       goto L_10;
                                   }
                               }
                               if (atom.atomnum[katm] == 7) // N-C=N+
                               {
                                  jj_bo = atom.tbo[katm];
                                  noxide = FALSE;
                                  for (l=0; l < atom.nconnect[katm]; l++)
                                  {
                                      if (atom.atomnum[atom.iat[katm][l]] == 8)
                                      {
                                          latm = atom.iat[katm][l];
                                          jjk = atom.nconnect[latm];
                                          if (jjk == 1)
                                            noxide = TRUE;
                                      }
                                  }
                                  if (jj_bo == 4  && noxide == FALSE)   // N-C=N+ 
                                  {
				    if ((atom.flags[katm] & aromatic_mask) && is_ring61(katm) ) // pyridinium type
				      {
					mmfftype = 40;
					goto L_10;
				      } else
				      {
					mmfftype = 55;
					goto L_10;
				      }
                                  }
                               }
                           }
                       }
                       if (atom.iat[i][j] != 0 && atom.atomnum[atom.iat[i][j]] == 7) // n-n=n
                       {
                           jatm = atom.iat[i][j];
                           for (k=0; k < atom.nconnect[jatm]; k++)
                           {
                               if (atom.iat[jatm][k] != 0 && atom.iat[jatm][k] != i)
                               {
                                   if (atom.bo[jatm][k] == 2 && atom.atomnum[atom.iat[jatm][k]] == 7)
                                   {
                                       jbo = atom.tbo[jatm];
                                       katm = atom.iat[jatm][k];
                                       kbo = atom.tbo[katm];
                                       if (jbo == 3 && kbo == 3)
                                       {
                                          mmfftype = 10;
                                          goto L_10;
                                       }
                                   }
                               }
                           }
                       }
                       
                   }
		   // sulfonamide
                   for(j=0; j < atom.nconnect[i]; j++)
                   {
                       jatm = atom.iat[i][j];
                       if (atom.atomnum[jatm] == 16)   // sulfonamide - need to check this further
                       {
                           jjk = 0;
                           for (k=0; k < atom.nconnect[jatm]; k++)
                           {
                               if (atom.atomnum[atom.iat[jatm][k]] == 8 && atom.bo[jatm][k] == 2)
                                  jjk++;
                           }
                           if (jjk >= 2)
                           {
                              mmxtype = 9;
                              mmfftype = 43;
                              goto L_10;
                           }
                         }
                         if (atom.atomnum[jatm] == 5) // boron
                         {
                           mmxtype = 9;
                           mmfftype = 9;
                           goto L_10;
                         }
                   }
		   // amide and thioamide
                   for(j=0; j < atom.nconnect[i]; j++)
                   {
                     if (atom.iat[i][j] != 0)
                     {
                       jatm = atom.iat[i][j];
                         for (k=0; k < atom.nconnect[jatm]; k++)
                         {
                            if (atom.iat[jatm][k] != 0 && atom.bo[jatm][k] == 2)
                            {
                               if (atom.atomnum[jatm] == 6 && atom.atomnum[atom.iat[jatm][k]] == 8) // amide
                               {
                                   mmxtype = 9;
                                   mmfftype = 10;
                                   goto L_10;
                               }
                               if (atom.atomnum[jatm] == 6 && atom.atomnum[atom.iat[jatm][k]] == 16) // thioamide
                               {
                                   mmxtype = 9;
                                   mm3type = 9;
                                   mmfftype = 10;
                                   goto L_10;
                               }
                            }
                         }
                     }
                   }
                   // eneamine
                       for (j=0; j < jji; j++)
                       {
                           if (atom.iat[i][j] != 0 && atom.bo[i][j] != 9)
                           {
                               if (atom.atomnum[atom.iat[i][j]] == 6)
                               {
                                   jatm = atom.iat[i][j];
                                   for (k=0; k < atom.nconnect[jatm]; k++)
                                   {
                                       if (atom.iat[jatm][k] != 0 && atom.bo[jatm][k] == 2)
                                       {
                                          mmxtype = 9;
                                          mmfftype = 40;
                                          goto L_10;
                                       }
                                   }
                               }
                           }
                       }
		       // n-p=x
                   for(j=0; j < jji; j++)
                   {
                     if (atom.iat[i][j] != 0)
                     {
                       jatm = atom.iat[i][j];
                         for (k=0; k < atom.nconnect[jatm]; k++)
                         {
                            if (atom.iat[jatm][k] != 0 && atom.bo[jatm][k] == 2)
                            {                         
                               if (atom.atomnum[jatm] == 15) // n-p=x
                               {
                                   jjk = atom.nconnect[jatm];
                                   if (jjk > 3)
                                   {
                                       mmxtype = 8;
                                       mmfftype = 8;
                                       goto L_10;
                                   }
                               }
                            }
                          }
                       }
                   }
               }
// end of ndouble test
	       // n-s and n-b
                
		   // sulfonamide
               for (j=0; j < jji; j++)
               {
                   if (atom.iat[i][j] != 0 && atom.bo[i][j] != 9)
                   {             
                       if (atom.atomnum[atom.iat[i][j]] == 6) // carbon
                       {
                           for (k=0; k < atom.nconnect[atom.iat[i][j]]; k++)
                           {
                               if (atom.iat[atom.iat[i][j]][k] != 0 && atom.iat[atom.iat[i][j]][k] != i)
                               {
                                   if (atom.bo[atom.iat[i][j]][k] == 3 && atom.atomnum[atom.iat[atom.iat[i][j]][k]] == 7)
                                   {
                                       mmfftype = 43;
                                       goto L_10;
                                   }
                               }
                           }
                       } else if (atom.atomnum[atom.iat[i][j]] == 7) // nitrogen - n-n bond  - n-n-c=o
                       {
                           jatm = atom.iat[i][j];
                           for (k=0; k < atom.nconnect[jatm]; k++)
                           {
                               if (atom.iat[jatm][k] != 0 && atom.iat[jatm][k] != i)
                               {
                                   katm = atom.iat[jatm][k];
                                   if (atom.atomnum[katm] == 6)
                                   {
                                       for (l=0; l < atom.nconnect[katm]; l++)
                                       {
                                           if (atom.iat[katm][l] != 0 && atom.iat[katm][l] != jatm)
                                           {
                                               if (atom.atomnum[atom.iat[katm][l]] == 8 && atom.bo[katm][l] == 2)
                                               {
                                                   mmxtype = 9;
                                                   goto L_10;
                                               }
                                           }
                                       }
                                   }
                               }
                           }
                       }
                   }
               }
               for (j=0; j < jji; j++)  // guanadinium
               {
                   if (atom.iat[i][j] != 0 && atom.atomnum[atom.iat[i][j]] == 6)
                   {
                       jatm = atom.iat[i][j];
                       jjk = atom.nconnect[jatm];
                       if (jjk == 3 && (atom.atomnum[atom.iat[jatm][0]] == 7 &&
                             atom.atomnum[atom.iat[jatm][1]] == 7 && atom.atomnum[atom.iat[jatm][2]] == 7) )
                       {
                           mmfftype = 56;
                           goto L_10;
                       }
                   }
               }
               goto L_10;           
          } else if (atom.atomnum[i] == 8) // oxygen
          {
/* ========================  Oxygen  ==================================== */
              if (atom.mmx_type[i] == 53)  // TS Oxygen
              {
                  mmxtype = 53;
                  mm3type = 0;
                  mmfftype = 0;
                  goto L_10;
              }
	      jji = atom.nconnect[i];
              jjbo = atom.tbo[i];
              if (jjbo == 3)  // o+
              {
                  mmxtype = 46;
                  mmfftype = 49;
                  mm3type = 0;
                  for (j=0; j < jji; j++)
                  {
                      if (atom.iat[i][j] != 0 && atom.bo[i][j] == 2)
                      {
                          mmfftype = 51;
                          goto L_10;
                      }
                  }
                  goto L_10;
              }
              mmxtype = atom.mmx_type[i];
              if (mmxtype == 66)
              {
                mm3type = 47;
                mmfftype = 32;
		gafftype = 14;
              }
              if (mmxtype == 6 || mmxtype == 7 || mmxtype == 66)
              {
                  for (j=0; j < jji; j++)
                  {
                      if (atom.bo[i][j] == 2)  // x=0
                      {
                          mmxtype = 7;
                          mm3type = 7;
                          mmfftype = 7;
			  gafftype = 14;
                          jatm = atom.iat[i][j];
                          if (atom.atomnum[jatm] == 15)   // p=o
                          {
                              for (k=0; k < atom.nconnect[jatm]; k++)
                              {
                                  jjk = 0;
                                  if (atom.iat[jatm][k] != 0 && atom.iat[jatm][k] != i && atom.atomnum[atom.iat[jatm][k]] == 8)
                                  {
                                      katm = atom.iat[jatm][k];
				      jjk = atom.nconnect[katm];
                                      if (jjk == 1 && atom.mmx_type[katm] == 42)
                                      {
                                         mmxtype = 66;
                                         mm3type = 47;
                                      }
                                  }
                              }
                              mmfftype = 32;
                              goto L_10;
                          }
                          if (atom.atomnum[jatm] == 16)   // s=o
                          {
                              for (k=0; k < atom.nconnect[jatm]; k++)
                              {
                                  jjk = 0;
                                  jj_bo = 0;
                                  if (atom.iat[jatm][k] != 0 && atom.iat[jatm][k] != i && atom.atomnum[atom.iat[jatm][k]] == 8)
                                  {
                                      katm = atom.iat[jatm][k];
				      jjk = atom.nconnect[katm];
				      jj_bo = atom.tbo[katm];
                                      if (jjk == 1 && jj_bo == 1 && atom.mmx_type[katm] == 42)
                                      {
                                         mmxtype = 66;
                                         mm3type = 47;
                                         mmfftype = 32;
                                         goto L_10;
                                      }
                                  }
                              }
                          }
                          if (atom.atomnum[jatm] == 6)    // c=o
                          {
                              mmfftype = 7;
                              for (k=0; k < atom.nconnect[jatm]; k++)
                              {
                                  if (atom.iat[jatm][k] != 0 && atom.iat[jatm][k] != i && atom.bo[jatm][k] != 9)
                                  {
                                      if (atom.atomnum[atom.iat[jatm][k]] == 8)
                                      {
                                          katm = atom.iat[jatm][k];
                                          jjk = atom.nconnect[katm];
                                          if (jjk == 1)   // carboxylate
                                          {
                                              mmxtype = 66;
                                              mm3type = 47;
                                              mmfftype = 32;
					      gafftype = 14;
                                              goto L_10;
                                          }
                                          if (jjk == 2)  // ester and acids
                                          {
                                              if ( atom.atomnum[atom.iat[katm][0]] != 1 && atom.atomnum[atom.iat[katm][1]] != 1) // ester
                                              {
                                                  mm3type = 78;
                                                  goto L_10;
                                              } else if ( atom.atomnum[atom.iat[katm][0]] == 1 || atom.atomnum[atom.iat[katm][1]] == 1) // acid
                                              {
                                                  mm3type = 77;
                                                  goto L_10;
                                              }
                                          }
                                      }
                                  }
                              }
                              for (k=0; k < atom.nconnect[jatm]; k++)
                              {
                                  if (atom.iat[jatm][k] != 0 && atom.iat[jatm][k] != i && atom.bo[jatm][k] != 9)
                                  {
                                      if (atom.atomnum[atom.iat[jatm][k]] == 7)  // amide
                                      {
                                          mm3type = 7;  // should be 79 but database only has 58-79
                                          goto L_10;
                                      }
                                  }
                              }
                              for (k=0; k < atom.nconnect[jatm]; k++)
                              {
                                  if (atom.iat[jatm][k] != 0 && atom.iat[jatm][k] != i && atom.bo[jatm][k] != 9)
                                  {
                                      if (atom.atomnum[atom.iat[jatm][k]] == 6) // vinyl ketone
                                      {
                                          katm = atom.iat[jatm][k];
                                          for (l=0; l < atom.nconnect[katm]; l++)
                                          {
                                              if (atom.iat[katm][l] != 0 && atom.bo[katm][l] == 2)
                                              {
                                                  mm3type = 81;
                                                  goto L_10;
                                              }
                                          }
                                      }
                                  }
                              }
                          }           
		                     
                          for (k=0; k < atom.nconnect[jatm]; k++)
                          {
                              if (atom.atomnum[atom.iat[jatm][k]] == 8 && atom.iat[jatm][k] != i && atom.bo[jatm][k] == 1)
                              {
                                  latm = atom.iat[jatm][k];
				  jjk = atom.nconnect[latm];
                                 if (jjk == 1)
                                  {
                                      mmxtype = 66;
                                      mm3type =  47;
                                      mmfftype = 32;
				      gafftype = 14;
                                      goto L_10;
                                  }
                              } else if (atom.atomnum[atom.iat[jatm][k]] == 8 && atom.iat[jatm][k] != i && atom.bo[jatm][k] == 2)
                              {
                                  mmfftype = 32;
                                  goto L_10;
                              } else if (atom.atomnum[atom.iat[jatm][k]] == 16 && atom.iat[jatm][k] != i && atom.bo[jatm][k] == 2)
                              {
                                  jjk = 0;
                                  latm = atom.iat[jatm][k];
				  jjk = atom.nconnect[latm];
                                  if (jjk == 1)
                                  {
                                      mmxtype = 66;
                                      mm3type = 47;
                                  }
                                  mmfftype = 32;
                                  goto L_10;
                              } else if (atom.atomnum[atom.iat[jatm][k]] == 7 && atom.iat[jatm][k] != i && atom.bo[jatm][k] == 2)
                              {
                                  mmfftype = 32;
                                  goto L_10;
                              }
                          }    
                          goto L_10;
		      } else if (atom.bo[i][j] == 3)
                      {
                          mmxtype = 46;
                          mmfftype = 53;
                          goto L_10;
                      }
                  } // got here with only single bonds - don't reset 66
                  jji = atom.nconnect[i];
                  nh = 0;
		  gafftype = 14;
                  for (j=0; j < jji; j++)
                  {
                      if (atom.atomnum[atom.iat[i][j]] == 1)
                          nh++;
                  }
                  if (jji == 1) // only one bond
                  {
                      if (atom.atomnum[atom.iat[i][0]] == 6)
                      {
                          jatm = atom.iat[i][0];
                          for (j=0; j < atom.nconnect[jatm]; j++)
                          {
                              if (atom.iat[jatm][j] != 0 && atom.bo[jatm][j] == 2)
                              {
                                  if (atom.atomnum[atom.iat[jatm][j]] == 6) // c=c
                                  {
                                     mmfftype = 35;
                                     goto L_10;
                                  }
                              }
                          }
                      }
                      if (atom.atomnum[atom.iat[i][0]] == 7)  //  n-o
                      {
                          mmxtype = 66;
                          mm3type = 69;
                          mmfftype = 32;
                          goto L_10;
                      }
                      if (atom.atomnum[atom.iat[i][0]] == 15 || atom.atomnum[atom.iat[i][0]] == 16)  //  p-o
                      {
                          mmxtype = 66;
                          mm3type = 7;
                          mmfftype = 32;
                          goto L_10;
                      }
                      
                  }
                  // test for epoxides
                  if (is_ring31(i))
                  {
		      gafftype = 16;
                      mm3type = 49;
                      mmfftype = 6;
                      goto L_10;
                  } 
                  if (is_ring51(i))
                  {
		    gafftype = 16;
                     if (atom.flags[i] & aromatic_mask)
                     {
                         mmfftype = 59;
                         mm3type = 41;
                         goto L_10;
                     }
                  }
                  if (jji == 2)
                  {
                      if (nh == 2) // mmff water
                      {
                          mmfftype = 70;
                          goto L_10;
                      } else if (nh == 1) // ROH and RCOOH
                      {
                          mm3type = 6;
                          mmfftype = 6;
		          gafftype = 15;
                          jatm = atom.iat[i][0];
                          for (j=0; j < atom.nconnect[jatm]; j++)
                          {
                              if (atom.iat[jatm][j] != 0 && atom.bo[jatm][j] == 2)
                              {
                                  if (atom.atomnum[atom.iat[jatm][j]] == 8)
                                  {
                                      mm3type = 75;
                                      break;
                                  }
                              }
                          }
                          goto L_10;
                      } else
                      {
                          mmfftype = 6;
//                          goto L_10;                             
                      }
                      jatm = atom.iat[i][0];
                      jtype = 0;
                      katm = atom.iat[i][1];
                      ktype = 0;
                      for (j=0; j < atom.nconnect[jatm]; j++)
                      {
                          if (atom.iat[jatm][j] != 0 && atom.bo[jatm][j] == 2)
                          {
                              if (atom.atomnum[atom.iat[jatm][j]] == 8)
                              {
                                 jtype = 3;
                                 break;
                              } else if (atom.atomnum[atom.iat[jatm][j]] == 6)
                              {
                                  jtype = 2;
                                  break;
                              }
                          }
                      }
                      for (j=0; j < atom.nconnect[katm]; j++)
                      {
                          if (atom.iat[katm][j] != 0 && atom.bo[katm][j] == 2)
                          {
                              if (atom.atomnum[atom.iat[katm][j]] == 8)
                              {
                                 ktype = 3;
                                 break;
                              }else if (atom.atomnum[atom.iat[katm][j]] == 6)
                              {
                                  ktype = 2;
                                  break;
                              }
                          }
                      }
                      if (jtype == 3 && ktype == 3) // anhydrides
                         mm3type = 148;
                      else if (jtype == 3 || ktype == 3)  // carboxyl & ester
                         mm3type = 75;
                      else if (jtype == 2 && ktype == 2) // furan type ??
                         mm3type = 41;
                      else if (jtype == 2 || ktype == 2) // vinyl type ??
                         mm3type = 41;
                      else
                         mm3type = 6;
                      goto L_10;
                  }
                  if (mmxtype != 66)
                  {
                     mmxtype = 6;
                     mm3type = 6;
                     mmfftype = 6;
		     gafftype = 16;
                  }
                  goto L_10;
              } else if (mmxtype == 42)
              {
                  mmfftype = 32;
		  gafftype = 14;
                  for (j=0; j < atom.nconnect[i]; j++)
                  {
                      if (atom.atomnum[atom.iat[i][j]] == 16 || atom.atomnum[atom.iat[i][j]] == 15)
                      {
                          mmxtype = 66;
                          mm3type = 47;
                          mmfftype = 32;
                          goto L_10;
                      }
                      if (atom.atomnum[atom.iat[i][j]] == 7)
                      {
                          mmxtype = 42;
                          mm3type = 69;
			  jatm = atom.iat[i][j];
                          jj_bo = atom.tbo[jatm];
                          if (jj_bo >= 4)
                          {
                              mmxtype = 66;
                              mm3type = 69;
                              mmfftype = 32;
                          } else
                            mmfftype = 35;
                         goto L_10; 
                      }
                      if (atom.atomnum[atom.iat[i][j]] == 6 )
                      {
                          jatm = atom.iat[i][j];
                          mmfftype = 35;
                          for (k=0; k < atom.nconnect[jatm]; k++)
                          {
                              if (atom.iat[jatm][k] != i && atom.atomnum[atom.iat[jatm][k]] == 8 && atom.bo[jatm][k] == 2)
                              {
                                  mmxtype = 66;
                                  mm3type = 47;
                                  mmfftype = 32;
                                  goto L_10;
                              }
                          }
                      }
                  }
                  goto L_10;
              } else
                 goto L_10;
          } else if (atom.atomnum[i] == 13) // aluminum
          {
              mmxtype = 44;
              mm3type = 0;
              jji = atom.nconnect[i];
              if (jji == 4)
                 mmxtype = 58;
              goto L_10;
          } else if (atom.atomnum[i] == 15) // phosphorus
          {
              mmxtype = 25;
              mm3type = 25;
              mmfftype = 25;
              jji = atom.nconnect[i];
              for (j=0; j < atom.nconnect[i]; j++)
              {
                   if (atom.atomnum[atom.iat[i][j]] == 6 && atom.bo[i][j] == 2)
                   {
                      mmxtype = 67;
                      mmfftype = 75;
                      goto L_10;
                   }
                   if (atom.atomnum[atom.iat[i][j]] == 8 && atom.bo[i][j] == 2)
                      mm3type = 153;
              }
              if (jji >= 5)
              {
                 mmxtype = 47;
                 mm3type = 60;
              }else if (jji == 3)
              {
                  mmxtype = 25;
                  mmfftype = 26;
              }
              goto L_10;
          } else if (atom.atomnum[i] == 16) // sulfur
          {
              mmxtype = 15;
              mm3type = 15;
              mmfftype = 15;
              jji = atom.nconnect[i];
              ndouble = 0;
              for (j=0; j < jji; j++)
              {
                  if (atom.iat[i][j] != 0 && atom.bo[i][j] != 9 && atom.mmx_type[atom.iat[i][j]] != 20 &&
                  (atom.mmx_type[atom.iat[i][j]] < 300))
                  {
                      if (atom.bo[i][j] == 2)
                         ndouble++;
                  }
              }
              if (jji == 1 && ndouble == 1)
              {
                  if (atom.atomnum[atom.iat[i][0]] == 16)  // s=s-
                  {
                      mmfftype = 72;
                      goto L_10;
                  }
                  if (atom.atomnum[atom.iat[i][0]] == 6)  // s=c
                  {
                      mmxtype = 38;
                      mm3type = 74;
                      if (atom.bo[i][0] == 2)
                      {
                         mmfftype = 16;
                         jatm = atom.iat[i][0];
                          for (k=0; k < atom.nconnect[jatm]; k++)
                          {
                              if (atom.iat[jatm][k] != 0 && atom.bo[jatm][k] != 9)
                              {
                                  if (atom.atomnum[atom.iat[jatm][k]] == 16 && atom.iat[jatm][k] != i)
                                  {
                                      katm = atom.iat[jatm][k];
                                      jj_bo = atom.tbo[katm];
                                      if (jj_bo == 1)  //  s=c-s    i-jatm-katm
                                      {
                                          mmfftype = 72;
                                          goto L_10;
                                      }
                                  }
                              }
                          }
                          goto L_10;
                      }else if (atom.bo[i][0] == 1)
                      {
                          mmfftype = 72;
                          goto L_10;
                      }
                      goto L_10;
                  }
                  if (atom.atomnum[atom.iat[i][0]] == 15 ) // s=p
                  {
                      mmxtype = 38;
                      mm3type = 74;
                      mmfftype = 72;
                      goto L_10;
                  }                  
              } else if (jji == 1 && ndouble == 0)  // s-c=s
	      {
		if (atom.input_charge[i] == -1) // negative s-c
		  {
		    mmfftype = 72;
		    goto L_10;
		  }
		  jatm = atom.iat[i][0];
		  for (k=0; k < atom.nconnect[jatm]; k++)
		    {
		      if (atom.bo[jatm][k] == 2)
			{
			  if (atom.atomnum[atom.iat[jatm][k]] == 16)
			    {
			      mmfftype = 72;
			      goto L_10;
			    }
			}
		    }
		  goto L_10;
	      }
              if (jji == 2)
              {
                  mmxtype = 15;
                  mm3type = 15;
                  mmfftype = 15;
                  if (is_ring51(i))  // thiophene
                  {
                      if (atom.flags[i] & aromatic_mask)
                      {
                          mmfftype = 44;
                          mm3type = 42;
                          goto L_10;
                      }
                  }
                  if (ndouble == 2)
                      mmfftype = 74;
                  goto L_10;
              }
              if (jji == 3)
              {
                  mmxtype = 17;
                  mm3type = 17;
                  mmfftype = 17;
                  if (ndouble == 2)
                    mmfftype = 73;
                  if (ndouble == 3)
                    mmfftype = 18;
                  goto L_10;
              }
              if (jji == 4)
              {
                  mmxtype = 18;
                  mm3type = 18;
                  for (j=0; j < jji; j++)
                  {
                      if (atom.atomnum[atom.iat[i][j]] == 7 ) // sulfamide
                         mm3type = 154;
                  }
                  mmfftype = 18;
                  goto L_10;
              }
              goto L_10;
          } else if (atom.atomnum[i] == 17) // chlorine
          {
              mm3type = 12;
              mmxtype = 12;
              mmfftype = 12;
	      if (atom.input_charge[i] == -1)
		mmfftype = 90;
	      gafftype = 33;
              jji = atom.nconnect[i];
              if (jji == 2)
                 mmxtype = 74;  // bridging chlorine
              else if (jji == 4) // perchlorate
                 mmfftype = 77;
              goto L_10;         
          } else if (atom.atomnum[i] == 34) // selenium
          {
              mm3type = 34;
              mmxtype = 34;
              for (j=0; j < jji; j++)
              {
                  if (atom.bo[i][j] == 2)
                     mmxtype = 39;
              }
              goto L_10;
          } else if (atom.atomnum[i] == 9)  // Florine
          {
              mmxtype = 11;
              mm3type = 11;
              mmfftype = 11;
	      gafftype = 32;
	      if (atom.input_charge[i] == -1)
		mmfftype = 89;
              goto L_10;
          } else if (atom.atomnum[i] == 35)   // Bromine
          {
              mmxtype = 13;
              mm3type = 13;
              mmfftype = 13;
	      gafftype = 34;
	      if (atom.input_charge[i] == -1)
		mmfftype = 91;
              goto L_10;
          } else if (atom.atomnum[i] == 53)    // Iodine
          {
              if (atom.mmx_type[i] == 54)  // Sn2 I 
              {
                 mmxtype = 54;
                 mm3type = 0;
                 mmfftype = 0;
                 goto L_10;
              }
              mmxtype = 14;
              mm3type = 14;
              mmfftype = 14;
	      gafftype = 35;
              goto L_10;
          } else if (atom.atomnum[i] == 14 )  // Silicon
          {
              mmxtype = 19;
              mm3type = 19;
              mmfftype = 19;
              goto L_10;
          } else if (atom.atomnum[i] == 50 ) // Tin
          {
              mmxtype = 32;
              mm3type = 32;
              mmfftype = 0;
              goto L_10;
          } else if (mmxtype >= 300) // metal atom - try to assign MMFF type
          {
              mm3type = mmxtype;
              if (mmfftype == 0)
                 mmfftype = mmxtype;
          }
L_10:
          set_atomtype(i,mmxtype,mm3type,mmfftype);
L_20:
          continue;  // do nothing
      }
   set_atomtypes(get_field());
   }
}
/* --------------------------------------------------------- */
/* ==================================================  */
// look for ionic types input with bonds and adjust to make ionic
//
void adjust_mmfftypes()
{
  int i,j, iatt,jji;

    for (i=1; i <= natom; i++)
    {
        jji = atom.nconnect[i];
        if (atom.mmff_type[i] == 89 ) // F-
        {
            for (j=0; j < jji; j++)
            {
                if (atom.iat[i][j] != 0)
                    deletebond(i,atom.iat[i][j]);
            }
        } else if (atom.mmff_type[i] == 90 ) // CL-
        {
            for (j=0; j < jji; j++)
            {
                if (atom.iat[i][j] != 0)
                    deletebond(i,atom.iat[i][j]);
            }
        } else if (atom.mmff_type[i] == 91 ) // Br-
        {
            for (j=0; j < jji; j++)
            {
                if (atom.iat[i][j] != 0)
                    deletebond(i,atom.iat[i][j]);
            }
        } else if (atom.mmff_type[i] == 92 ) // Li+
        {
            for (j=0; j < jji; j++)
            {
                if (atom.iat[i][j] != 0)
                {
                    iatt = atom.iat[i][j];
                    deletebond(i,iatt);
                    if (atom.atomnum[iatt] == 8)
                       set_atomtype(iatt,42,0,92);
                }
            }
        } else if (atom.mmff_type[i] == 93 ) // Na+
        {
            for (j=0; j < jji; j++)
            {
                if (atom.iat[i][j] != 0)
                {
                    iatt = atom.iat[i][j];
                    deletebond(i,iatt);
                    if (atom.atomnum[iatt] == 8)
                       set_atomtype(iatt,42,0,93);
                }
            }
        } else if (atom.mmff_type[i] == 94 ) // K+
        {
            for (j=0; j < jji; j++)
            {
                if (atom.iat[i][j] != 0)
                {
                    iatt = atom.iat[i][j];
                    deletebond(i,iatt);
                    if (atom.atomnum[iatt] == 8)
                       set_atomtype(iatt,42,0,94);
                }
            }
        } else if (atom.mmff_type[i] == 95 ) // Zn+2
        {
            for (j=0; j < jji; j++)
            {
                if (atom.iat[i][j] != 0)
                {
                    iatt = atom.iat[i][j];
                    deletebond(i,iatt);
                    if (atom.atomnum[iatt] == 8)
                       set_atomtype(iatt,42,0,95);
                }
            }
        } else if (atom.mmff_type[i] == 96 ) // Ca+2
        {
            for (j=0; j < jji; j++)
            {
                if (atom.iat[i][j] != 0)
                {
                    iatt = atom.iat[i][j];
                    deletebond(i,iatt);
                    if (atom.atomnum[iatt] == 8)
                       set_atomtype(iatt,42,0,96);
                }
            }
        }
    }
}
// ==============================
// coung number of electron withdrawing groups attached to carbon
int count_ewg(int jatm) 
{
  int i, newg;

  newg = 0;
  for (i=0; i < atom.nconnect[jatm]; i++)
    {
      if (atom.atomnum[jatm] == 7 || atom.atomnum[jatm] == 8 || atom.atomnum[jatm] == 15 || atom.atomnum[jatm] == 16)
	newg++;
      if (atom.atomnum[jatm] == 9 || atom.atomnum[jatm] == 17 || atom.atomnum[jatm] == 35 ||  atom.atomnum[jatm] == 53)
	newg++;
    }
  return (newg);
}
    
