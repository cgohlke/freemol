#define EXTERN extern

#include "pcwin.h"
#include "torsions.h"
#include "atom_k.h"
#include "bonds_ff.h"

int isbond(int,int);
int is_delocalbond(int,int);
void numeral(int,char *,int);
void four(char *,char *,char *,char *, char *);
int find_bond(int,int);
int is_ring42(int,int);
int is_ring54(int, int, int, int);
void message_alert(char *, char *);
int have_ring3(void);
int have_ring4(void);
int have_ring5(void);
int have_ring6(void);
int get_field(void);
void ktorsion(int *type,int *tclass,int **iat,int **bo);
char *get_structure_title(void);
        
EXTERN struct t_allene {
    int nallene, ntor[10];
     } allene;
        
        
EXTERN struct t_torkn1 {
        int  use_tor4, use_tor5;
        int  ntor, ntor4, ntor5, ntordel, torindex[MAXTORDEL];
        char  kv[MAXTORCONST][13], kv4[MAXTOR4CONST][13], kv5[MAXTOR5CONST][13];
        char  kvdel[MAXTORDEL][13];
        float tv1[MAXTORCONST], tv2[MAXTORCONST], tv3[MAXTORCONST];
        float tv4[MAXTORCONST], tv5[MAXTORCONST], tv6[MAXTORCONST];
        int   phase1[MAXTORCONST], phase2[MAXTORCONST], phase3[MAXTORCONST];
        int   phase4[MAXTORCONST], phase5[MAXTORCONST], phase6[MAXTORCONST];
        float tv41[MAXTOR4CONST], tv42[MAXTOR4CONST], tv43[MAXTOR4CONST];
        int   phase41[MAXTOR4CONST], phase42[MAXTOR4CONST], phase43[MAXTOR4CONST];
        float tv51[MAXTOR5CONST], tv52[MAXTOR5CONST], tv53[MAXTOR5CONST];
        int   phase51[MAXTOR5CONST], phase52[MAXTOR5CONST], phase53[MAXTOR5CONST];
        float tvdel1[MAXTORDEL], tvdel2[MAXTORDEL], tvdel3[MAXTORDEL];
        int   phasedel1[MAXTORDEL], phasedel2[MAXTORDEL], phasedel3[MAXTORDEL];
        } torkn1;
        
EXTERN int Missing_constants;

void ktorsion(int *type,int *tclass,int **iat,int **bo)
{
  int i, j, ierr, itor,field;
    int ia, ib, ic, id;
    int ita, itb, itc, itd;
    int cl_a, cl_b, cl_c, cl_d;
    int use_ring4, use_ring5;
    int nb1,nb2,nb3, nbo1,nbo2,nbo3;
    char izero[4];
    char pa[4],pb[4],pc[4],pd[4],pt[13], kv1[13];
    char k1[13], k2[13], k3[13], k4[13];

    itor = 0;
    nb1 = nb2 = nb3 = 0;
    nbo1 = nbo2 = nbo3 = 0;
    use_ring4 = FALSE;
    if ( (have_ring4()) && torkn1.ntor4 > 0)
       use_ring4 = TRUE;

    use_ring5 = FALSE;
    if ( (have_ring5()) && torkn1.ntor5 > 0)
       use_ring5 = TRUE;
    field = get_field();
    for (i=0; i < torsions.ntor; i++)
    {
        ia = torsions.i14[i][0];
        ib = torsions.i14[i][1];
        ic = torsions.i14[i][2];
        id = torsions.i14[i][3];

        ita = type[ia];
        itb = type[ib];
        itc = type[ic];
        itd = type[id];

        cl_a = tclass[ia];
        cl_b = tclass[ib];
        cl_c = tclass[ic];
        cl_d = tclass[id];

        numeral(ita,pa,3);
        numeral(itb,pb,3);
        numeral(itc,pc,3);
        numeral(itd,pd,3);
        
        if (itb < itc )
           four(pa,pb, pc, pd, pt);
        else if (itb > itc)
           four(pd,pc, pb, pa, pt);
        else if (ita < itd)
           four(pa,pb, pc, pd, pt);
        else
           four(pd,pc, pb, pa, pt);
        strcpy(izero,"  0");

        if (field == MMFF94)
        {

/* get bond index each bond for MMFF94  */
           nb1 = find_bond(ia,ib);
           nb2 = find_bond(ib,ic);
           nb3 = find_bond(ic,id);
           nbo1 = bonds_ff.index[nb1];
           nbo2 = bonds_ff.index[nb2];
           nbo3 = bonds_ff.index[nb3];

           // pt is specific angle
           // stage 2 in step down
             numeral(atom_k.tclass[ita],pa,3);
             numeral(atom_k.tclass[itb],pb,3);
             numeral(atom_k.tclass[itc],pc,3);
             numeral(atom_k.tclass[itd],pd,3);
           if (atom_k.tclass[itb] < atom_k.tclass[itc])
             {
               four(pa,pb,pc,pd,k1);
             } else if (atom_k.tclass[itc] < atom_k.tclass[itb])
             {
               four(pd,pc,pb,pa,k1);
             }else
             {
                 if (atom_k.tclass[ita] < atom_k.tclass[itd])
                    four(pa,pb,pc,pd,k1);
                 else
                   four(pd,pc,pb,pa,k1);
             }
           // stage 3 in step down
             numeral(atom_k.tclass1[ita],pa,3);
             numeral(atom_k.tclass[itb],pb,3);
             numeral(atom_k.tclass[itc],pc,3);
             strcpy(pd,izero);
           if (atom_k.tclass[itb] < atom_k.tclass[itc])
             {
               four(pa,pb,pc,pd,k2);
             } else if (atom_k.tclass[itc] < atom_k.tclass[itb])
             {
               four(pd,pc,pb,pa,k2);
             }else
             {
                 if (atom_k.tclass1[ita] < 0)
                    four(pa,pb,pc,pd,k2);
                 else
                   four(pd,pc,pb,pa,k2);
             }
           // stage 4 in step down
             numeral(atom_k.tclass1[itd],pd,3);
             numeral(atom_k.tclass[itb],pb,3);
             numeral(atom_k.tclass[itc],pc,3);
             strcpy(pa,izero);
           if (atom_k.tclass[itb] < atom_k.tclass[itc])
             {
               four(pa,pb,pc,pd,k3);
             } else if (atom_k.tclass[itc] < atom_k.tclass[itb])
             {
               four(pd,pc,pb,pa,k3);
             }else
             {
                 if (atom_k.tclass1[itd] < 0)
                    four(pd,pc,pb,pa,k3);
                 else
                   four(pa,pb,pc,pd,k3);
             }
           // stage 5
             strcpy(pa,izero);
             strcpy(pd,izero);
           if (atom_k.tclass[itb] < atom_k.tclass[itc])
             {
               four(pa,pb,pc,pd,k4);
             } else if (atom_k.tclass[itc] < atom_k.tclass[itb])
             {
               four(pd,pc,pb,pa,k4);
             } else
                four(pa,pb,pc,pd,k4);

        }
        
        ierr = FALSE;
//  check for four membered rings
        if ( isbond(ia,id) )
        {
          for(j=0; j < torkn1.ntor4; j++)
          {
             strcpy(kv1,torkn1.kv4[j]);
             if (strcmp(pt,kv1) == 0 )
             {
               torsions.v1[i] =  torkn1.tv41[j];
               torsions.ph1[i] = torkn1.phase41[j];     
               torsions.v2[i] =  torkn1.tv42[j];
               torsions.ph2[i] = torkn1.phase42[j];     
               torsions.v3[i] =  torkn1.tv43[j];
               torsions.ph3[i] = torkn1.phase43[j];
               ierr = TRUE;
               goto L_10;
             }  
          }
          for(j=0; j < torkn1.ntor4; j++)
          {
             strcpy(kv1,torkn1.kv4[j]);
             if (strcmp(k1,kv1) == 0 )
             {
               torsions.v1[i] =  torkn1.tv41[j];
               torsions.ph1[i] = torkn1.phase41[j];     
               torsions.v2[i] =  torkn1.tv42[j];
               torsions.ph2[i] = torkn1.phase42[j];     
               torsions.v3[i] =  torkn1.tv43[j];
               torsions.ph3[i] = torkn1.phase43[j];
               ierr = TRUE;
               goto L_10;
             }  
          }
          for(j=0; j < torkn1.ntor4; j++)
          {
             strcpy(kv1,torkn1.kv4[j]);
             if (strcmp(k2,kv1) == 0 )
             {
               torsions.v1[i] =  torkn1.tv41[j];
               torsions.ph1[i] = torkn1.phase41[j];     
               torsions.v2[i] =  torkn1.tv42[j];
               torsions.ph2[i] = torkn1.phase42[j];     
               torsions.v3[i] =  torkn1.tv43[j];
               torsions.ph3[i] = torkn1.phase43[j];
               ierr = TRUE;
               goto L_10;
             }  
          }
          for(j=0; j < torkn1.ntor4; j++)
          {
             strcpy(kv1,torkn1.kv4[j]);
             if (strcmp(k3,kv1) == 0 )
             {
               torsions.v1[i] =  torkn1.tv41[j];
               torsions.ph1[i] = torkn1.phase41[j];     
               torsions.v2[i] =  torkn1.tv42[j];
               torsions.ph2[i] = torkn1.phase42[j];     
               torsions.v3[i] =  torkn1.tv43[j];
               torsions.ph3[i] = torkn1.phase43[j];
               ierr = TRUE;
               goto L_10;
             }  
          }
          for(j=0; j < torkn1.ntor4; j++)
          {
             strcpy(kv1,torkn1.kv4[j]);
             if (strcmp(k4,kv1) == 0 )
             {
               torsions.v1[i] =  torkn1.tv41[j];
               torsions.ph1[i] = torkn1.phase41[j];     
               torsions.v2[i] =  torkn1.tv42[j];
               torsions.ph2[i] = torkn1.phase42[j];     
               torsions.v3[i] =  torkn1.tv43[j];
               torsions.ph3[i] = torkn1.phase43[j];
               ierr = TRUE;
               goto L_10;
             }  
          }
        }
//   delocalized torsions 
        if ( nbo1 ==1 || nbo2 == 1 || nbo3 == 1 )
        {
          for(j=0; j < torkn1.ntordel; j++)
          {
             strcpy(kv1,torkn1.kvdel[j]);
             if (strcmp(pt,kv1) == 0  )
             {
               torsions.v1[i] =  torkn1.tvdel1[j];
               torsions.ph1[i] = torkn1.phasedel1[j];     
               torsions.v2[i] =  torkn1.tvdel2[j];
               torsions.ph2[i] = torkn1.phasedel2[j];     
               torsions.v3[i] =  torkn1.tvdel3[j];
               torsions.ph3[i] = torkn1.phasedel3[j];
               ierr = TRUE;
               nbo1 = bonds_ff.index[nb1];
               nbo2 = bonds_ff.index[nb2];
               nbo3 = bonds_ff.index[nb3];

               //              printf("delocal tor: %d %d %d %d = %d %d %d ; %d\n",ia,ib,ic,id,nbo1,nbo2,nbo3,torkn1.torindex[j]);
               if (nbo2 == 1 && torkn1.torindex[j] == 1)
                 goto L_10;
               if (nbo2 == 0 && (nbo1 == 1 || nbo3 == 1) && torkn1.torindex[j] == 2)
                 goto L_10;
             }  
          }
          for(j=0; j < torkn1.ntordel; j++)
          {
             strcpy(kv1,torkn1.kvdel[j]);
             if (strcmp(k1,kv1) == 0  )
             {
               torsions.v1[i] =  torkn1.tvdel1[j];
               torsions.ph1[i] = torkn1.phasedel1[j];     
               torsions.v2[i] =  torkn1.tvdel2[j];
               torsions.ph2[i] = torkn1.phasedel2[j];     
               torsions.v3[i] =  torkn1.tvdel3[j];
               torsions.ph3[i] = torkn1.phasedel3[j];
               ierr = TRUE;
               nbo1 = bonds_ff.index[nb1];
               nbo2 = bonds_ff.index[nb2];
               nbo3 = bonds_ff.index[nb3];

               //              printf("delocal tor: %d %d %d %d = %d %d %d ; %d\n",ia,ib,ic,id,nbo1,nbo2,nbo3,torkn1.torindex[j]);
               if (nbo2 == 1 && torkn1.torindex[j] == 1)
                 goto L_10;
               if (nbo2 == 0 && (nbo1 == 1 || nbo3 == 1) && torkn1.torindex[j] == 2)
                 goto L_10;
             }  
          }
          for(j=0; j < torkn1.ntordel; j++)
          {
             strcpy(kv1,torkn1.kvdel[j]);
             if (strcmp(k2,kv1) == 0  )
             {
               torsions.v1[i] =  torkn1.tvdel1[j];
               torsions.ph1[i] = torkn1.phasedel1[j];     
               torsions.v2[i] =  torkn1.tvdel2[j];
               torsions.ph2[i] = torkn1.phasedel2[j];     
               torsions.v3[i] =  torkn1.tvdel3[j];
               torsions.ph3[i] = torkn1.phasedel3[j];
               ierr = TRUE;
               nbo1 = bonds_ff.index[nb1];
               nbo2 = bonds_ff.index[nb2];
               nbo3 = bonds_ff.index[nb3];

               //              printf("delocal tor: %d %d %d %d = %d %d %d ; %d\n",ia,ib,ic,id,nbo1,nbo2,nbo3,torkn1.torindex[j]);
               if (nbo2 == 1 && torkn1.torindex[j] == 1)
                 goto L_10;
               if (nbo2 == 0 && (nbo1 == 1 || nbo3 == 1) && torkn1.torindex[j] == 2)
                 goto L_10;
             }  
          }
          for(j=0; j < torkn1.ntordel; j++)
          {
             strcpy(kv1,torkn1.kvdel[j]);
             if (strcmp(k3,kv1) == 0  )
             {
               torsions.v1[i] =  torkn1.tvdel1[j];
               torsions.ph1[i] = torkn1.phasedel1[j];     
               torsions.v2[i] =  torkn1.tvdel2[j];
               torsions.ph2[i] = torkn1.phasedel2[j];     
               torsions.v3[i] =  torkn1.tvdel3[j];
               torsions.ph3[i] = torkn1.phasedel3[j];
               ierr = TRUE;
               nbo1 = bonds_ff.index[nb1];
               nbo2 = bonds_ff.index[nb2];
               nbo3 = bonds_ff.index[nb3];

               //              printf("delocal tor: %d %d %d %d = %d %d %d ; %d\n",ia,ib,ic,id,nbo1,nbo2,nbo3,torkn1.torindex[j]);
               if (nbo2 == 1 && torkn1.torindex[j] == 1)
                 goto L_10;
               if (nbo2 == 0 && (nbo1 == 1 || nbo3 == 1) && torkn1.torindex[j] == 2)
                 goto L_10;
             }  
          }
          for(j=0; j < torkn1.ntordel; j++)
          {
             strcpy(kv1,torkn1.kvdel[j]);
             if (strcmp(k4,kv1) == 0  )
             {
               torsions.v1[i] =  torkn1.tvdel1[j];
               torsions.ph1[i] = torkn1.phasedel1[j];     
               torsions.v2[i] =  torkn1.tvdel2[j];
               torsions.ph2[i] = torkn1.phasedel2[j];     
               torsions.v3[i] =  torkn1.tvdel3[j];
               torsions.ph3[i] = torkn1.phasedel3[j];
               ierr = TRUE;
               nbo1 = bonds_ff.index[nb1];
               nbo2 = bonds_ff.index[nb2];
               nbo3 = bonds_ff.index[nb3];

               //              printf("delocal tor: %d %d %d %d = %d %d %d ; %d\n",ia,ib,ic,id,nbo1,nbo2,nbo3,torkn1.torindex[j]);
               if (nbo2 == 1 && torkn1.torindex[j] == 1)
                 goto L_10;
               if (nbo2 == 0 && (nbo1 == 1 || nbo3 == 1) && torkn1.torindex[j] == 2)
                 goto L_10;
             }  
          }
        }
//   check for five membered rings
         if (is_ring54(ia,ib,ic,id) )
         {
          for(j=0; j < torkn1.ntor5; j++)
          {
             strcpy(kv1,torkn1.kv5[j]);
             if (strcmp(pt,kv1) == 0 )
             {
               torsions.v1[i] =  torkn1.tv51[j];
               torsions.ph1[i] = torkn1.phase51[j];     
               torsions.v2[i] =  torkn1.tv52[j];
               torsions.ph2[i] = torkn1.phase52[j];     
               torsions.v3[i] =  torkn1.tv53[j];
               torsions.ph3[i] = torkn1.phase53[j];
               ierr = TRUE;
               goto L_10;
             }  
          }
          for(j=0; j < torkn1.ntor5; j++)
          {
             strcpy(kv1,torkn1.kv5[j]);
             if (strcmp(k1,kv1) == 0 )
             {
               torsions.v1[i] =  torkn1.tv51[j];
               torsions.ph1[i] = torkn1.phase51[j];     
               torsions.v2[i] =  torkn1.tv52[j];
               torsions.ph2[i] = torkn1.phase52[j];     
               torsions.v3[i] =  torkn1.tv53[j];
               torsions.ph3[i] = torkn1.phase53[j];
               ierr = TRUE;
               goto L_10;
             }  
          }
          for(j=0; j < torkn1.ntor5; j++)
          {
             strcpy(kv1,torkn1.kv5[j]);
             if (strcmp(k2,kv1) == 0 )
             {
               torsions.v1[i] =  torkn1.tv51[j];
               torsions.ph1[i] = torkn1.phase51[j];     
               torsions.v2[i] =  torkn1.tv52[j];
               torsions.ph2[i] = torkn1.phase52[j];     
               torsions.v3[i] =  torkn1.tv53[j];
               torsions.ph3[i] = torkn1.phase53[j];
               ierr = TRUE;
               goto L_10;
             }  
          }
          for(j=0; j < torkn1.ntor5; j++)
          {
             strcpy(kv1,torkn1.kv5[j]);
             if (strcmp(k3,kv1) == 0 )
             {
               torsions.v1[i] =  torkn1.tv51[j];
               torsions.ph1[i] = torkn1.phase51[j];     
               torsions.v2[i] =  torkn1.tv52[j];
               torsions.ph2[i] = torkn1.phase52[j];     
               torsions.v3[i] =  torkn1.tv53[j];
               torsions.ph3[i] = torkn1.phase53[j];
               ierr = TRUE;
               goto L_10;
             }  
          }
          for(j=0; j < torkn1.ntor5; j++)
          {
             strcpy(kv1,torkn1.kv5[j]);
             if (strcmp(k4,kv1) == 0 )
             {
               torsions.v1[i] =  torkn1.tv51[j];
               torsions.ph1[i] = torkn1.phase51[j];     
               torsions.v2[i] =  torkn1.tv52[j];
               torsions.ph2[i] = torkn1.phase52[j];     
               torsions.v3[i] =  torkn1.tv53[j];
               torsions.ph3[i] = torkn1.phase53[j];
               ierr = TRUE;
               goto L_10;
             }  
          }
         }
//   check regular parameters
         if (field == MMFF94)
         {
             for(j=0; j < torkn1.ntor; j++)   // specific constants
             {
                 strcpy(kv1,torkn1.kv[j]);
                 if (strcmp(pt,kv1) == 0)
                 {
                     torsions.v1[i] = torkn1.tv1[j];
                     torsions.ph1[i] = torkn1.phase1[j];     
                     torsions.v2[i] = torkn1.tv2[j];
                     torsions.ph2[i] = torkn1.phase2[j];     
                     torsions.v3[i] = torkn1.tv3[j];
                     torsions.ph3[i] = torkn1.phase3[j];     
                     torsions.v4[i] = torkn1.tv4[j];
                     torsions.ph4[i] = torkn1.phase4[j];     
                     torsions.v5[i] = torkn1.tv5[j];
                     torsions.ph5[i] = torkn1.phase5[j];     
                     torsions.v6[i] = torkn1.tv6[j];
                     torsions.ph6[i] = torkn1.phase6[j];
                     ierr = TRUE;
                     goto L_10;
                 }
             }
             // no specific step down to class types
 
             for(j=0; j < torkn1.ntor; j++)  //  class 2 in MMFF
             {
                 strcpy(kv1,torkn1.kv[j]);
                 if (strcmp(k1,kv1) == 0 )
                 {
                     torsions.v1[i] = torkn1.tv1[j];
                     torsions.ph1[i] = torkn1.phase1[j];     
                     torsions.v2[i] = torkn1.tv2[j];
                     torsions.ph2[i] = torkn1.phase2[j];     
                     torsions.v3[i] = torkn1.tv3[j];
                     torsions.ph3[i] = torkn1.phase3[j];     
                     torsions.v4[i] = torkn1.tv4[j];
                     torsions.ph4[i] = torkn1.phase4[j];     
                     torsions.v5[i] = torkn1.tv5[j];
                     torsions.ph5[i] = torkn1.phase5[j];     
                     torsions.v6[i] = torkn1.tv6[j];
                     torsions.ph6[i] = torkn1.phase6[j];
                     ierr = TRUE;
                     goto L_10;
                 }
             }
             //

             for(j=0; j < torkn1.ntor; j++)  // stage 3
             {
                 strcpy(kv1,torkn1.kv[j]);
                 if (strcmp(k2,kv1) == 0  )
                 {
                     torsions.v1[i] = torkn1.tv1[j];
                     torsions.ph1[i] = torkn1.phase1[j];     
                     torsions.v2[i] = torkn1.tv2[j];
                     torsions.ph2[i] = torkn1.phase2[j];     
                     torsions.v3[i] = torkn1.tv3[j];
                     torsions.ph3[i] = torkn1.phase3[j];     
                     torsions.v4[i] = torkn1.tv4[j];
                     torsions.ph4[i] = torkn1.phase4[j];     
                     torsions.v5[i] = torkn1.tv5[j];
                     torsions.ph5[i] = torkn1.phase5[j];     
                     torsions.v6[i] = torkn1.tv6[j];
                     torsions.ph6[i] = torkn1.phase6[j];
                     ierr = TRUE;
                     goto L_10;
                 }
             }


             for(j=0; j < torkn1.ntor; j++)  // stage 4
             {
                 strcpy(kv1,torkn1.kv[j]);
                 if (strcmp(k3,kv1) == 0  )
                 {
                     torsions.v1[i] = torkn1.tv1[j];
                     torsions.ph1[i] = torkn1.phase1[j];     
                     torsions.v2[i] = torkn1.tv2[j];
                     torsions.ph2[i] = torkn1.phase2[j];     
                     torsions.v3[i] = torkn1.tv3[j];
                     torsions.ph3[i] = torkn1.phase3[j];     
                     torsions.v4[i] = torkn1.tv4[j];
                     torsions.ph4[i] = torkn1.phase4[j];     
                     torsions.v5[i] = torkn1.tv5[j];
                     torsions.ph5[i] = torkn1.phase5[j];     
                     torsions.v6[i] = torkn1.tv6[j];
                     torsions.ph6[i] = torkn1.phase6[j];
                     ierr = TRUE;
                     goto L_10;
                 }
             }

             for(j=0; j < torkn1.ntor; j++)  // stage 5
             {
                 strcpy(kv1,torkn1.kv[j]);
                 if (strcmp(k4,kv1) == 0  )
                 {
                     torsions.v1[i] = torkn1.tv1[j];
                     torsions.ph1[i] = torkn1.phase1[j];     
                     torsions.v2[i] = torkn1.tv2[j];
                     torsions.ph2[i] = torkn1.phase2[j];     
                     torsions.v3[i] = torkn1.tv3[j];
                     torsions.ph3[i] = torkn1.phase3[j];     
                     torsions.v4[i] = torkn1.tv4[j];
                     torsions.ph4[i] = torkn1.phase4[j];     
                     torsions.v5[i] = torkn1.tv5[j];
                     torsions.ph5[i] = torkn1.phase5[j];     
                     torsions.v6[i] = torkn1.tv6[j];
                     torsions.ph6[i] = torkn1.phase6[j];
                     ierr = TRUE;
                     goto L_10;
                 }
             }
         }              
//  missing parameters
          if (ierr == FALSE)
          {
//              Missing_constants = TRUE;
            fprintf(pcmlogfile,"Torsion constants missing: angle: %d %d %d %d  types: %d %d %d %d  %s\n",ia,ib,ic,id,
		    ita, itb, itc, itd,get_structure_title());
          }
L_10:
    continue;
    }

    if (allene.nallene > 0)
    {
        for (i = 0; i < allene.nallene; i++)
        {
            torsions.v1[allene.ntor[i]] = 0;
            torsions.v2[allene.ntor[i]] = -11.5;
            torsions.v3[allene.ntor[i]] = 0;
        }
    }
}
/*  --------------------------------------------   */
void four(char *pa, char *pb, char *pc, char *pd, char *pt)
{
        strcpy(pt,pa);
        strcat(pt,pb);
        strcat(pt,pc);
        strcat(pt,pd);
} 
