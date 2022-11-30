#define EXTERN extern
#include "pcwin.h"
#include "angles.h"

EXTERN struct t_ooplane_k {
         int nopbend;
         char iopb[MAXOOP][13];
         float copb[MAXOOP];
        } ooplane_k;

void numeral(int,char *,int);
int get_field(void);
void kopbend(int *type,int *tclass);
EXTERN int Missing_constants;

void kopbend(int *type,int *tclass)
{
  int i, j, ia,ib,ic,id, ierr,field;
   int ita,itb, itc,itd, itd_class;
   int it, it1;
   char pa[4],pb[4],pc[4],pd[4],pdc[4],pt[13],pt1[13],pt2[13];
 
   angles.nopb = 0;  
   field = get_field();
   for (i=0; i < angles.nang; i++)
   {
      ia = angles.i13[i][0];
      ib = angles.i13[i][1];
      ic = angles.i13[i][2];
      ita = type[ia];
      itb = type[ib];
      itc = type[ic];
      ierr = FALSE;
      if (angles.i13[i][3] != 0 )
      {
             // angle ia-ib-id-ic
             id = angles.i13[i][3];
             itd = type[id];
             itd_class = tclass[id];
             numeral(ita,pa,3);
             numeral(itb,pb,3);
             numeral(itc,pc,3);
             numeral(itd,pd,3);
             numeral(itd_class,pdc,3);
             if (field == MMX)
             {
                 strcpy(pt,"  0"); strcat(pt,pb); strcat(pt,pd); strcat(pt,"  0");
             }else if (field == MMFF94)
             {
               if (ita <= itc && ita <= itd)
               {
                 if (itc <= itd)
                 {
                   strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                 } else
                 {
                   strcpy(pt,pa); strcat(pt,pb); strcat(pt,pd); strcat(pt,pc);
                 }
               } else if ( itc <= ita && itc <= itd)
               {
                 if (ita <= itd)
                 {
                    strcpy(pt,pc); strcat(pt,pb); strcat(pt,pa); strcat(pt,pd);
                 } else
                 {
                    strcpy(pt,pc); strcat(pt,pb); strcat(pt,pd); strcat(pt,pa);
                 }
               } else if (itd <= ita && itd <= itc)
               {
                 if ( ita <= itc)
                 {
                   strcpy(pt,pd); strcat(pt,pb); strcat(pt,pa); strcat(pt,pc);
                 } else
                 {
                   strcpy(pt,pd); strcat(pt,pb); strcat(pt,pc); strcat(pt,pa);
                 }
               }
             }
             strcpy(pt1,"  0"); strcat(pt1,pb); strcat(pt1,pdc); strcat(pt1,"  0");
             strcpy(pt2,"  0"); strcat(pt2,pb); strcat(pt2,"  0"); strcat(pt2,"  0");
             it = itb*100 + itd;
             it1 = itb*100;
             ierr = FALSE;
             for (j=0; j < ooplane_k.nopbend; j++)
             {
                if (strcmp(pt,ooplane_k.iopb[j]) == 0 || strcmp(pt1,ooplane_k.iopb[j]) == 0 || strcmp(pt2,ooplane_k.iopb[j]) == 0)
                {
                     angles.copb[i] = ooplane_k.copb[j];
                     angles.angin[i] = TRUE;
                     angles.nopb++;
                     ierr = TRUE;
                     break;
                }
             }
      }
   }
}
