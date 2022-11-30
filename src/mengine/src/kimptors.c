#define EXTERN extern
#include "pcwin.h"
#include "pcmod.h"
#include "utility.h"
#include "angles.h"

void four(char *,char *,char *,char *,char *);
void numeral(int,char *,int);

     
EXTERN struct t_improptor_k {
        int nimptor, nimprop;
        char kv[MAXIMP][13];
        float cimptor[MAXIMP], tdi[MAXIMP];
        float v1[MAXIMP], v2[MAXIMP], v3[MAXIMP];
        int   ph1[MAXIMP], ph2[MAXIMP], ph3[MAXIMP];
        } improptor_k;

struct t_improp {
        int nimptors, **iiprop;
        float *v1, *v2, *v3;
        int   *ph1, *ph2, *ph3;        
        } improp; 
                
void kimptors(void);
        

void kimptors()
{
   int i, j, ia, ib, ic, id, k;
   int ita, itb, itc, itd;
   char pa[4],pb[4],pc[4],pd[4],it[6][13];
 
   improp.nimptors = 0;  
   improp.iiprop = imatrix(0,MAXOOP,0,4);
   improp.v1 = vector(0,MAXOOP);
   improp.v2 = vector(0,MAXOOP);
   improp.v3 = vector(0,MAXOOP);
   improp.ph1 = ivector(0,MAXOOP);
   improp.ph2 = ivector(0,MAXOOP);
   improp.ph3 = ivector(0,MAXOOP);

   for (i=0; i < angles.nang; i++)
   {
      ia = angles.i13[i][0];
      ib = angles.i13[i][1];
      ic = angles.i13[i][2];
 
      if (angles.i13[i][3] != 0)
      {
         id = angles.i13[i][3];
         ita = atom.tclass[ia];
         itb = atom.tclass[ib];
         itc = atom.tclass[ic];
         itd = atom.tclass[id];
         numeral(ita,pa,3);
         numeral(itb,pb,3);
         numeral(itc,pc,3);
         numeral(itd,pd,3);
         
         four( pa, pc, pb, pd, it[0] );
         four( pa, pd, pb, pc, it[1] );
         four( pc, pa, pb, pd, it[2] );
         four( pc, pd, pb, pa, it[3] );
         four( pd, pa, pb, pc, it[4] );
         four( pd, pc, pb, pa, it[5] );
         
         improp.iiprop[improp.nimptors][0] = ib;
         for (j=0; j < improptor_k.nimptor; j++)
         {
            for(k=0; k < 6; k++)
            {
               if (strcmp(it[k],improptor_k.kv[j]) == 0)
               {
                 if (k == 0 || k == 1) improp.iiprop[improp.nimptors][1] = ia;
                 if (k == 2 || k == 3) improp.iiprop[improp.nimptors][1] = ic;
                 if (k == 4 || k == 5) improp.iiprop[improp.nimptors][1] = id;
                 
                 if (k == 2 || k == 4) improp.iiprop[improp.nimptors][2] = ia;
                 if (k == 0 || k == 5) improp.iiprop[improp.nimptors][2] = ic;
                 if (k == 1 || k == 3) improp.iiprop[improp.nimptors][2] = id;
                                 
                 if (k == 3 || k == 5) improp.iiprop[improp.nimptors][3] = ia;
                 if (k == 1 || k == 4) improp.iiprop[improp.nimptors][3] = ic;
                 if (k == 0 || k == 2) improp.iiprop[improp.nimptors][3] = id;
                 improp.v1[improp.nimptors] = improptor_k.v1[j];
                 improp.v2[improp.nimptors] = improptor_k.v2[j];
                 improp.v3[improp.nimptors] = improptor_k.v3[j];
                 improp.ph1[improp.nimptors]  = improptor_k.ph1[j];
                 improp.ph2[improp.nimptors]  = improptor_k.ph2[j];
                 improp.ph3[improp.nimptors]  = improptor_k.ph3[j];
                 improp.nimptors++;
                 break;
               }
            }
         }
      }
  }
}
