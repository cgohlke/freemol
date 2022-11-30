#define EXTERN extern

#include "pcwin.h"

EXTERN struct t_minim_values {
        int iprint, ndc, nconst;
        float dielc;
        } minim_values;

void echarge(int natom,int *use, int **skip,double *x,double *y,double *z,double *charge,double chrgcut,double *eu);
void echarge1(int natom,int *use, int **skip,double *x,double *y,double *z,double *charge,double chrgcut,double *eu,double **deqq);
void echarge2(int i,int natom,int *use, int **skip,double *x,double *y,double *z,double *charge,double chrgcut,float **hessx,float **hessy,float **hessz);
       
void echarge(int natom,int *use, int **skip,double *x,double *y,double *z,double *charge,double chrgcut,double *eu)
{
   int i, j;
   double xi, yi, zi, xk,yk,zk;
   double xr, yr, zr;
   double e, rik, rik2,cc;
   double rdielc;
   double cutoff;

   cutoff = chrgcut*chrgcut;
   *eu = 0.0;
   cc = 4.80298*4.80298*14.39418;
   rdielc = 1.0/units.dielec;

   if (minim_values.iprint)
   {
       fprintf(pcmlogfile,"\nCharge-Charge Interactions : Diele = %7.3f\n",units.dielec);
       fprintf(pcmlogfile,"      At1      At2      q1        q2         Rik          Eqq\n");
   }
   for (i=1; i < natom; i++)
   {
      if (charge[i] != 0.0 )
      {
         xi = x[i];
         yi = y[i];
         zi = z[i];
         
         for(j=i+1; j <= natom; j++)
         {
           if (use[i] || use[j])
           {
            if ( skip[i][j] != i)
            {
             if (charge[j] != 0.0)
             {
                  
                xk = x[j];
                yk = y[j];
                zk = z[j];
                
                xr = xi - xk;
                yr = yi - yk;
                zr = zi - zk;
                rik2 =  xr*xr + yr*yr + zr*zr;
                if (rik2 < cutoff)
                {
                   rik = sqrt(rik2);
                   e = cc*rdielc*charge[i]*charge[j]/rik;
                   if (skip[i][j] == -i)
                     e /= units.chgscale;  
                   *eu += e;
                 
                  if (minim_values.iprint)
                     fprintf(pcmlogfile,"QQ: (%-3d)- (%-3d)   %-8.3f  %-8.3f   %-8.3f = %8.4f\n", i, j,charge[i], charge[j],rik,e);
                }
             }
            }
           }
         }
      }
   }
}
// ====================================================
void echarge1(int natom,int *use, int **skip,double *x,double *y,double *z,double *charge,double chrgcut,double *eu,double **deqq)
{
   int i,j;
   double xi, yi, zi;
   double xr, yr, zr;
   double xk,yk,zk;
   double e, rik, cc, rik2;
   double de, dedx, dedy, dedz;
   double rdielc;
   double cutoff;

   cutoff = chrgcut*chrgcut;
   cc = 4.80298*4.80298*14.39418;
   rdielc = 1.0/units.dielec;

   for (i=1; i < natom; i++)
   {
       if (charge[i] != 0.0 )
       {
         xi = x[i];
         yi = y[i];
         zi = z[i];
         
         for(j=i+1; j <= natom; j++)
         {
           if (use[i] || use[j])
           {
            if ( skip[i][j] != i)
            {
             if (charge[j] != 0.0)
             {
                  
                xk = x[j];
                yk = y[j];
                zk = z[j];
                xr = xi - xk;
                yr = yi - yk;
                zr = zi - zk;
                rik2 = xr*xr +yr*yr + zr*zr;
                if (rik2 < cutoff)
                {
                   rik = sqrt( rik2);
                   e = cc*rdielc*charge[i]*charge[j]/rik;
                   if (skip[i][j] == -i)
                     e /= units.chgscale;

                   de = -cc*rdielc*charge[i]*charge[j]/rik2;

                   de = de/rik;
                   dedx = de*xr;
                   dedy = de*yr;
                   dedz = de*zr;
                  
                   deqq[i][0] += dedx;
                   deqq[i][1] += dedy;
                   deqq[i][2] += dedz;
         
                   deqq[j][0] -= dedx;
                   deqq[j][1] -= dedy;
                   deqq[j][2] -= dedz;
                   *eu += e;
                }
               }
             }
           }
	 }
       }
   }
}
// ==========================================
void echarge2(int i,int natom,int *use, int **skip,double *x,double *y,double *z,double *charge,double chrgcut,float **hessx,float **hessy,float **hessz)
{
    int j,k;
    double fik,de,d2e,d2edx,d2edy,d2edz;
    double xi,yi,zi,xr,yr,zr,term[3][3];
    double xk,yk,zk;
    double r,r2;
    double cc, rdielc;
    double cutoff;
    
    if (charge[i] == 0.0)
       return;

    cutoff = chrgcut*chrgcut;
       
    cc = 4.80298*4.80298*14.39418;
    rdielc = 1.0/units.dielec;
      
    skip[i][i] = i;
    
    xi = x[i];
    yi = y[i];
    zi = z[i];

    
    for (k=1; k <= natom; k++)
    {
        if (i != k && charge[k] != 0.0  &&  (skip[i][k] != i) )
        {
            xk = x[k];
            yk = y[k];
            zk = z[k];

            xr = xi - xk;
            yr = yi - yk;
            zr = zi - zk;
            r2 = xr*xr + yr*yr + zr*zr;
            if (r2 < cutoff)
            {
              r = sqrt(r2);
              fik = cc*rdielc*charge[i]*charge[k];
              if (skip[i][k] == -i)
                 fik /= units.chgscale;
              de = -fik/r2;
              d2e = -2.0*de/r;
              de = de / r;
              d2e = (d2e-de) / r2;
              d2edx = d2e * xr;
              d2edy = d2e * yr;
              d2edz = d2e * zr;
              term[0][0]= d2edx*xr + de;
              term[1][0] = d2edx*yr;
              term[2][0] = d2edx*zr;
              term[0][1] = term[1][0];
              term[1][1] = d2edy*yr + de;
              term[2][1] = d2edy*zr;
              term[0][2] = term[2][0];
              term[1][2] = term[2][1];
              term[2][2] = d2edz*zr + de;

              for (j=0; j < 3; j++)
              {
                  hessx[i][j] += term[j][0];
                  hessy[i][j] += term[j][1];
                  hessz[i][j] += term[j][2];
                  hessx[k][j] -= term[j][0];
                  hessy[k][j] -= term[j][1];
                  hessz[k][j] -= term[j][2];
              }
            }
        }
    }
}
