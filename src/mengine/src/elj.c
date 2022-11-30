#define EXTERN extern

#include "pcwin.h"

void elj(int natom,int *type,int *use,double *x,double *y,double *z,double vdwcut,int **skip,double **vrad,double **veps,double *evdw,double *e14);
void elj1(int natom,int *type,int *use,double *x,double *y,double *z,double vdwcut,int **skip,double **vrad,double **veps,double *evdw,double *e14,
           double **devdw,double **de14);
void elj2(int iatom,int natom,int *type,int *use,double *x,double *y,double *z,double vdwcut,int **skip,double **vrad,double **veps,
           float **hessx,float **hessy,float **hessz);
       
EXTERN struct t_minim_values {
        int iprint, ndc, nconst;
        float dielc;
        } minim_values;
            
void elj(int natom,int *type,int *use,double *x,double *y,double *z,double vdwcut,int **skip,double **vrad,double **veps,double *evdw,double *e14)
{
      int i,ia, ib,j, ita, itb;
      double xi,yi,zi,xr,yr,zr;
      double xk,yk,zk;
      double e,p2,p6,p12,rv,eps;
      double rik2;
      double cutoff;

      cutoff = vdwcut*vdwcut;

      *evdw = 0.0;
      *e14  = 0.0;

      if (minim_values.iprint)
      {
          fprintf(pcmlogfile,"\nVDW Terms - Lennard-Jones Potential\n");
          fprintf(pcmlogfile,"     At1   At2    Rik    Radius     Eps       Evdw\n");
      }
      for (i=1; i < natom; i++)
      {
         for(j=i+1; j <= natom; j++)
         {
            if (skip[i][j] != i)
            {
              if (use[i] || use[j])
              {
                 ia = i;
                 ib = j;
                 ita = type[ia];
                 itb = type[ib];
                         
                 xi = x[ia];
                 yi = y[ia];
                 zi = z[ia];
                 xk = x[ib];
                 yk = y[ib];
                 zk = z[ib];

                 xr = xi - xk;
                 yr = yi - yk;
                 zr = zi - zk;
                 rik2 = xr*xr + yr*yr + zr*zr;
                 if ( rik2 < cutoff)
		   {
                     rv = vrad[ita][itb];
                     eps = veps[ita][itb];
                 
                     if (skip[i][j] == -i)
                       eps /= units.v14scale;
                     p2 = (rv*rv)/rik2;
                     p6 = p2*p2*p2;
                     p12 = p6*p6;
                     e = eps*(p12 - 2.0*p6);
                     if (skip[i][j] == -i)
                       *e14 += e;
                     else
                       *evdw += e;
                     if (minim_values.iprint)
                       fprintf(pcmlogfile,"VDW: %-4d - %-4d %-8.3f %-8.3f %-8.3f = %8.4f\n",ia, ib, sqrt(rik2),rv,eps,e);
                 }
              }
            }
         }
      }
}
// =====================================
void elj1(int natom,int *type,int *use,double *x,double *y,double *z,double vdwcut,int **skip,double **vrad,double **veps,double *evdw,double *e14,
           double **devdw,double **de14)
{
  int i,ia, ib, j, ita,itb;
      double xi,yi,zi,xr,yr,zr;
      double xk,yk,zk;
      double e,p2,p6,p12,rv,eps;
      double rik2, rik;
      double dedx,dedy,dedz,de;
      double cutoff;

      cutoff = vdwcut*vdwcut;
      *evdw = 0.0;
      *e14  = 0.0;

      for (i=1; i <= natom; i++)
      {
          devdw[i][0] = 0.0;
          devdw[i][1] = 0.0;
          devdw[i][2] = 0.0;
          de14[i][0] = 0.0;
          de14[i][1] = 0.0;
          de14[i][2] = 0.0;
      }


      for (i=1; i < natom; i++)
      {
         ia = i;
         ita = type[ia];
         xi = x[ia];
         yi = y[ia];
         zi = z[ia];
         for(j=i+1; j <= natom; j++)
         {
            if (use[i] || use[j])
            {
              if (skip[i][j] != i)
              {                
                 ib = j;
                 itb = type[j];
                 xk = x[ib];
                 yk = y[ib];
                 zk = z[ib];

                 xr = xi - xk;
                 yr = yi - yk;
                 zr = zi - zk;
                 rik2 = xr*xr + yr*yr + zr*zr;
                 if (rik2 < cutoff)
                 {
                     rv = vrad[ita][itb];
                     eps = veps[ita][itb];
                     if (skip[i][j] == -i)
                       eps /= units.v14scale;
                     rik = sqrt(rik2);
                     p2 = (rv*rv)/rik2;
                     p6 = p2*p2*p2;
                     p12 = p6*p6;
                     e = eps*(p12 - 2.0*p6);
                     de = eps * (p12 - p6) * (-12.0/rik);
                     de = de / rik;
                     dedx = de * xr;
                     dedy = de * yr;
                     dedz = de * zr;

                     if (skip[i][j] == -i)
                     {
                       *e14 += e;
                       de14[ia][0] += dedx;
                       de14[ia][1] += dedy;
                       de14[ia][2] += dedz;
                       de14[ib][0] -= dedx;
                       de14[ib][1] -= dedy;
                       de14[ib][2] -= dedz;
                     } else
                     {
                       *evdw += e;
                       devdw[ia][0] += dedx;
                       devdw[ia][1] += dedy;
                       devdw[ia][2] += dedz;
                       devdw[ib][0] -= dedx;
                       devdw[ib][1] -= dedy;
                       devdw[ib][2] -= dedz;
                     }
                 }
              }
            }
         }
      }
}
// ===============================
void elj2(int iatom,int natom,int *type,int *use,double *x,double *y,double *z,double vdwcut,int **skip,double **vrad,double **veps,
           float **hessx,float **hessy,float **hessz)
{
      int ia, ib,j, k, ita, itb;
      double xi,yi,zi,xr,yr,zr;
      double xk,yk,zk;
      double e,p2,p6,p12,rv,eps;
      double rik2, rik;
      double de,d2e;
      double d2edx,d2edy,d2edz,term[3][3];
      double cutoff;

      cutoff = vdwcut*vdwcut;      

    skip[iatom][iatom] = iatom;
    ia = iatom;
    ita = type[ia];
         
    xi = x[ia];
    yi = y[ia];
    zi = z[ia];
   for (j=1; j <= natom; j++)
   {
      ib = j;         
       if (skip[iatom][j] != iatom)
      {
          itb = type[ib];

         xk = x[ib];
         yk = y[ib];
         zk = z[ib];

         xr = xi - xk;
         yr = yi - yk;
         zr = zi - zk;
         rik2 = xr*xr + yr*yr + zr*zr;
         if (rik2 < cutoff)
         {
             rv = vrad[ita][itb];
             eps = veps[ita][itb];
             if (skip[iatom][j] == -iatom)
              eps /= units.v14scale;
             rik = sqrt(rik2);
             p2 = (rv*rv)/rik2;
             p6 = p2*p2*p2;
             p12 = p6*p6;
             e = eps*(p12 - 2.0*p6);
             de = eps * (p12 - p6) * (-12.0/rik);
             d2e = eps * (13.0*p12 - 7.0*p6) * (12.0/rik2);
             de = de / rik;
             d2e = (d2e-de) / rik2;
             d2edx = d2e * xr;
             d2edy = d2e * yr;
             d2edz = d2e * zr;
             term[0][0] = d2edx*xr + de;
             term[1][0] = d2edx*yr;
             term[2][0] = d2edx*zr;
             term[0][1] = term[1][0];
             term[1][1] = d2edy*yr + de;
             term[2][1] = d2edy*zr;
             term[0][2] = term[2][0];
             term[1][2] = term[2][1];
             term[2][2] = d2edz*zr + de;

           if (ia == iatom)
           {
                 for (k=0; k < 3; k++)
                 {
                     hessx[ia][k] += term[k][0];
                     hessy[ia][k] += term[k][1];
                     hessz[ia][k] += term[k][2];
                     hessx[ib][k] -= term[k][0];
                     hessy[ib][k] -= term[k][1];
                     hessz[ib][k] -= term[k][2];
                 }
           } else
           {
                 for (k=0; k < 3; k++)
                 {
                     hessx[ia][k] -= term[k][0];
                     hessy[ia][k] -= term[k][1];
                     hessz[ia][k] -= term[k][2];
                     hessx[ib][k] += term[k][0];
                     hessy[ib][k] += term[k][1];
                     hessz[ib][k] += term[k][2];
                 }
           }
         }
      }
   }
}
