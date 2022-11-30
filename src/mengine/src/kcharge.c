#define EXTERN extern

#include "pcwin.h"
#include "bonds_ff.h"
#include "atom_k.h"
#include "utility.h"
#include "job_control.h"

#define MAXOPLS  600

void numeral(int, char *, int);
int find_rsize(int,int);
void search_rings(int);
void message_alert(char *,char *);
void get_rsize(int,int,int, int *);
int get_field(void);
void compute_fcharge(int natom,int *type,int *atomnum,int **iat,int **bo,double *charge,double *formal_charge);
void kcharge(int natom,int *type,int *atomnum,long int *flags,int **iat,int **bo,double *charge,double *sigma_charge,double *formal_charge);

EXTERN struct ElementType { 
                        char symbol[3];
                        int   atomnum;
                        float weight, covradius, vdwradius;
                        int s,p,d,f, type;
                   } Elements[] ;

        
EXTERN struct t_charge_k {
         int ncharge, nbndchrg, nbndchrgdel;
         int type[MAXATOMTYPE], btype[MAXBONDCONST], btypedel[MAXBONDCONST];
         float charge[MAXATOMTYPE], bcharge[MAXBONDCONST], formchrg[MAXATOMTYPE], bchargedel[MAXBONDCONST];
         float typechrg[MAXATOMTYPE];
         } charge_k;
         
EXTERN struct  t_dipole_k {
        int ndipole,ndipole3,ndipole4,ndipole5;
        char   kb[MAXBONDCONST][7], kb3[MAXBOND3CONST][7], kb4[MAXBOND4CONST][7], kb5[MAXBOND5CONST][7];
        float bmom[MAXBONDCONST],bmom3[MAXBOND3CONST],bmom4[MAXBOND4CONST],bmom5[MAXBOND5CONST];
         } dipole_k;

void kcharge(int natom,int *type,int *atomnum,long int *flags,int **iat,int **bo,double *charge,double *sigma_charge,double *formal_charge)
{
  int i, j, ia, ib, it, iit, kit, jji,field;
   long int mask;
   char pa[4],pb[4],pt[7];
   float bmomj;
   
   mask = 1L << 0;
   field = get_field();
   for (i=1; i <= natom; i++)
   {
      sigma_charge[i] = 0.0;
      if (field == MMX && (type[i] < 300))
        charge[i] = 0.0;
      else if (field != MMX && type[i] < 300)
        charge[i] = 0.0;
   }

   if (field == MMFF94)
   {
     compute_fcharge(natom,type,atomnum,iat,bo,charge,formal_charge);
       for(i=0; i < bonds_ff.nbnd; i++)
       {
           ia = bonds_ff.i12[i][0];
           ib = bonds_ff.i12[i][1];
           iit = type[ia];
           kit = type[ib];
           if (iit < kit)
             it = iit*100+kit;
           else
             it = kit*100+iit;
           if (bonds_ff.index[i] == 1)
           {
              for(j=0; j < charge_k.nbndchrgdel; j++)
              {
                 if (it == charge_k.btypedel[j])
                 {
                   if (iit < kit)
                   {
                     charge[ia] -= charge_k.bchargedel[j];
                     charge[ib] += charge_k.bchargedel[j];
                   }else
                   {
                     charge[ib] -= charge_k.bchargedel[j];
                     charge[ia] += charge_k.bchargedel[j];
                   }
                   break;
                 }
              }
           } else
           {
              for(j=0; j < charge_k.nbndchrg; j++)
              {
                 if (it == charge_k.btype[j])
                 {
                   if (iit < kit)
                   {
                     charge[ia] -= charge_k.bcharge[j];
                     charge[ib] += charge_k.bcharge[j];
                   }else
                   {
                     charge[ib] -= charge_k.bcharge[j];
                     charge[ia] += charge_k.bcharge[j];
                   }
                   break;
                 }
             }
           }
       }
       // add in formal charges
       for (i=1; i <= natom; i++)
       {
           if (atomnum[i] == 15)
               formal_charge[i] = 0.0;
           if (atomnum[i] == 16 && type[i] != 72)
               formal_charge[i] = 0.0;             
       }              
       for (i=1; i <= natom; i++)
       {
           charge[i] += (1.0 - atom_k.ligands[type[i]]*charge_k.formchrg[type[i]])*formal_charge[i] ;
           for(j=0; j < MAXIAT; j++)
           {
               if (iat[i][j] != 0)
               {
                  if (formal_charge[iat[i][j]] < 0.0)
		    {
                           charge[i] += charge_k.formchrg[type[iat[i][j]]]*formal_charge[iat[i][j]];
		    }
               }
           }
       }
   } else if (field == MMX)
   {
      // compute charges from bond moments
      for (i=0; i < bonds_ff.nbnd; i++)
      {
          ia = bonds_ff.i12[i][0];
          ib = bonds_ff.i12[i][1];
          iit = type[ia];
          kit = type[ib];
         if( field == MMX)
         {
             if (iit >= 300)
                iit = 300;
             if( kit >= 300)
                kit = 300;
             if (iit == 40 && ( kit != 2 && kit != 3 && kit != 4 && kit != 40) )
                iit = 2;
             if (kit == 40 && ( iit != 2 && iit != 3 && iit != 4 && iit != 40) )
                kit = 2;
         }
          numeral(iit,pa,3);
          numeral(kit,pb,3);
          
          if (iit < kit)
          {
              strcpy(pt,pa);
              strcat(pt,pb);
          } else
          {
              strcpy(pt,pb);
              strcat(pt,pa);
          }

              bmomj = 0.0;
              for (j=0; j < dipole_k.ndipole; j++)
              {
                  if (strcmp(pt,dipole_k.kb[j]) == 0)
                  {
                      bmomj = dipole_k.bmom[j];
                      break;
                  }
              }
          
          if (bmomj != 0.0)
          {
              bmomj = 0.42*bmomj/bonds_ff.bl[i];
              if (iit <= kit)
              {
                  charge[ia] += 0.5*bmomj;
                  charge[ib] -= 0.5*bmomj;
              } else
              {
                  charge[ia] -= 0.5*bmomj;
                  charge[ib] += 0.5*bmomj;
              }
          }
          if (kit == 66)
          {
              if (iit == 3)
                charge[ib] -= 0.5;
              if (iit == 18)
                charge[ib] -= 0.33;
          } else if (iit == 66)
          {
              if (kit == 3)
                charge[ia] -= 0.5;
              if (kit == 18)
                charge[ia] -= 0.33;
          }
      }
      // check for isolated charged atoms
      for (i=1; i <= natom; i++)
      {
          if ( !(flags[i] & mask))
          {
              if (type[i] == 16 || type[i] == 30)
                  charge[i] += 1.0;
              if (type[i] == 41)  // N+
              {
                  for(j=0; j < MAXIAT; j++)
                  {
                      if (iat[i][j] != 0)
                      {
                          if (type[iat[i][j]] == 42 || type[iat[i][j]] == 4
                          || type[iat[i][j]] == 66)
                             goto L_10;
                      }
                  }
                  charge[i] += 1.0;
              }
              if (type[i] == 46) // O+
              {
                  for (j=0; j < MAXIAT; j++)
                  {
                      if (bo[i][j] == 3) // co in metal carbonyl
                        goto L_10;
                  }                 
                  charge[i] += 1.0;
              }
              if (type[i] == 27 || type[i] == 42 || type[i] == 48)
                 charge[i] -= 1.0;
              if (type[i] >= 11 && type[i] <= 14)
              {
                  jji = 0;
                  for (j=0; j < MAXIAT; j++)
                  {
                      if (iat[i][j] != 0)
                        jji++;
                  }
                  if (jji == 0)
                    charge[i] -= 1.0;
              }
              if (type[i] == 25 || type[i] == 47)
              {
                  jji = 0;
                  for (j=0; j < MAXIAT; j++)
                  {
                      if (iat[i][j] != 0)
                        jji++;
                  }
                  if (jji == 6)
                    charge[i] -= 1.0;
              }
              if (type[i] == 58) // Al
                 charge[i] += 1.0;
                     
          }
L_10:
          continue;
      }
      // search list of charged atoms read in
      for (i=1; i <= natom; i++)
      {
        for(j=0; j < charge_k.ncharge; j++)
        {
          if (type[i] == charge_k.type[j])
          {
              charge[i] = charge_k.charge[j];
              break;
          }                   
        }
      }     
   } else
   {
      // search list of charges
       for (i=1; i <= natom; i++)
       {
           for(j=0; j < charge_k.ncharge; j++)
           {
               if (type[i] == charge_k.type[j])
               {
                   charge[i] = charge_k.charge[j];
                   break;
               }                   
           }
       }
   }
   for (i=1;i <= natom; i++)
      sigma_charge[i] = charge[i];

   // scale or turn off charges based on control from PYMOL
   if (job_control.use_charge)
     {
       for (i=1; i <= natom; i++)
	 charge[i] *= job_control.scale;
     }
}
/* -------------------------------------------------------- */
void compute_fcharge(int natom,int *type,int *atomnum,int **iat,int **bo,double *charge,double *formal_charge)
{
    int i, j, jjbo;
    int jatm, k, jji, attached[4];
    int *used;
    float temp_charge;
    
    used = ivector(1,natom+1);
    for (i=1; i <= natom; i++)
       used[i] = FALSE;
       
    for (i=1; i <= natom; i++)
    {
        formal_charge[i] = charge_k.typechrg[type[i]];
    }
//  find variable charges: types 32, 72, 76, 81
    for (i=1; i <= natom; i++)
    {
        if (used[i] == FALSE)
        {
            if (formal_charge[i] != 0.0)
            {
               if (type[i] == 32)
               {
                   jatm = iat[i][0];
                   jji = 0;
                   jjbo = 0;
                   for (k=0; k < MAXIAT; k++)
                   {
                       if (iat[jatm][k] != 0 && bo[jatm][k] != 9)
                       {
                           if (type[iat[jatm][k]] == 32)
                           {
                               if (bo[jatm][k] == 1)
                               {
                                   jjbo++;
                                   attached[jji] = iat[jatm][k];
                                   jji++;
                               } else
                               {
                                   attached[jji] = iat[jatm][k];
                                   jji++;
                               }
                           }
                       }
                   }
                   temp_charge = -(double)jjbo/jji;
                   if (jji == 1)  // only 1 type 32
                     formal_charge[i] = 0.0;
                   else if (jji == 2)
                   {
                       if (atomnum[jatm] == 15)  // p-o  treat as single negative charge;
                          temp_charge = -1.00/jji;
                       if (type[jatm] == 73)
                       {
                           temp_charge = -0.500;
                           formal_charge[attached[0]] = temp_charge;
                           formal_charge[attached[1]] = temp_charge;
                           used[attached[0]] = TRUE;                           
                           used[attached[1]] = TRUE;
                       } else if (atomnum[jatm] != 7)
                       {
                           formal_charge[attached[0]] = temp_charge;
                           formal_charge[attached[1]] = temp_charge;
                           used[attached[0]] = TRUE;                           
                           used[attached[1]] = TRUE;
                       }else
                       {
                           formal_charge[attached[0]] = 0.0;
                           formal_charge[attached[1]] = 0.0;
                           used[attached[0]] = TRUE;                           
                           used[attached[1]] = TRUE;
                       }                         
                   } else if (jji == 3)
                   {
                        formal_charge[attached[0]] = temp_charge;
                        formal_charge[attached[1]] = temp_charge;
                        formal_charge[attached[2]] = temp_charge;
                        used[attached[0]] = TRUE;                           
                        used[attached[1]] = TRUE;                           
                        used[attached[2]] = TRUE;                           
                   } else if (jji == 4)
                   {
                        formal_charge[attached[0]] = temp_charge;
                        formal_charge[attached[1]] = temp_charge;
                        formal_charge[attached[2]] = temp_charge;
                        formal_charge[attached[3]] = temp_charge;
                        used[attached[0]] = TRUE;                           
                        used[attached[1]] = TRUE;                           
                        used[attached[2]] = TRUE;                           
                        used[attached[3]] = TRUE;                           
                   }
               }else if (type[i] == 72)
               {
                   jatm = iat[i][0];
                   jji = 0;
                   jjbo = 0;
                   for (k=0; k < MAXIAT; k++)
                   {
                       if (iat[jatm][k] != 0 && bo[jatm][k] != 9)
                       {
                           if (type[iat[jatm][k]] == 72)
                           {
                               if (bo[jatm][k] == 1)
                               {
                                   jjbo++;
                                   attached[jji] = iat[jatm][k];
                                   jji++;
                               } else
                               {
                                   attached[jji] = iat[jatm][k];
                                   jji++;
                               }
                           }
                       }
                   }
                   if (jji == 1)
                   {
                       if (atomnum[jatm] == 15)
                          formal_charge[i] = 0.0;
                       else
                          formal_charge[i] = -1.00;
                   }else
                   {
                       formal_charge[attached[0]] = -0.5;
                       formal_charge[attached[1]] = -0.5;
                       used[attached[0]] = TRUE;
                       used[attached[1]] = TRUE;
                   }
	       } else if (type[i] == 76)
	       {
		 jji = 1;
		 attached[0] = i;
		 for (k=i+1; k < natom; k++)
		   {
		     if (type[k] == 76)
		       {
			 attached[jji] = k;
			 jji++;
		       }
		   }
		 if (jji == 1)
		   {
		     temp_charge = -1.00;
		     formal_charge[i] = temp_charge;
		     used[i] = TRUE;
		   } else
		   {
		     temp_charge = -1.00/jji;
		     for (j=0; j < jji; j++)
		       {
			 formal_charge[attached[j]] = temp_charge;
			 used[attached[j]] = TRUE;
		       }
		   }
               } else if (type[i] == 81)
               {
                   jji = 0;
                   for (k=0; k < MAXIAT ; k++)
                   {
                       if (iat[i][k] != 0 && type[iat[i][k]] == 80)
                       {
                           jatm = iat[i][k];
                           jji = 0;
                           for (j=0; j < MAXIAT; j++)
                           {
                               if (iat[jatm][j] != 0 && bo[jatm][j] != 9)
                               {
                                   if (atomnum[iat[jatm][j]] == 7)
                                   {
                                       attached[jji] = iat[jatm][j];
                                       jji++;
                                   }
                               }
                           }
                       }
                   }
                   if (jji == 0)
                   {
                     temp_charge = 1.00;
                     formal_charge[i] = temp_charge;
                     used[i] = TRUE;
                   } else
                   {
                      temp_charge = 1.00/jji;
                      for (j=0; j < jji; j++)
                      {
                          formal_charge[attached[j]] = temp_charge;
                          used[attached[j]] = TRUE;
                      }
                   }
               }
            }
        }
    }
    free_ivector(used,1, natom+1);    
}
