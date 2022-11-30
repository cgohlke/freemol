#define EXTERN extern

#include "pcwin.h"
                
EXTERN struct t_minim_values {
        int iprint, ndc, nconst;
        float dielc;
        } minim_values;
        
int find_bond(int,int);
// ebond(nbnd,i12,type,tclass,use,x,y,z,bl,bk,*estr)
// ebond1(nbnd,i12,type,tclass,use,x,y,z,bl,bk,*estr,**destr)
void ebond(int nbnd,int **,int *use,double *x,double *y,double *z,double *bl,double *bk,double *estr);
void ebond1(int natom,int nbnd,int **,int *use,double *x,double *y,double *z,double *bl,double *bk,double *estr,double **deb);
void ebond2(int ia,int **,int **,int **,double *x,double *y,double *z,double *bl,double *bk,float **hessx,float **hessy,float **hessz);

//  =========================== 
void ebond(int nbnd,int **i12,int *use,double *x,double *y,double *z,double *bl,double *bk,double *estr)
{
/* compute stretching energy  */

   int i, it, kt;
   double xr, yr, zr, rik, rik2, bcorr;
   double dt, dt2, e;
 
   *estr = 0.0F;

   if (minim_values.iprint)
   {
       fprintf(pcmlogfile,"\nBond Terms \n");
       fprintf(pcmlogfile,"        At1       At2     R       BLen    Bconst      Eb\n");
   }
     
   for (i=0; i < nbnd; i++)
   {
      it = i12[i][0];
      kt = i12[i][1];
      bcorr = 0.0;
      if ( use[it] || use[kt] )
      {
         xr = x[it] - x[kt];
         yr = y[it] - y[kt];
         zr = z[it] - z[kt];
         rik2 = xr*xr + yr*yr + zr*zr;
         rik = sqrt(rik2);
         dt = rik - (bl[i] - bcorr);
         dt2 = dt*dt;
         e = units.bndunit *bk[i]*dt2*
            (1.0 + units.cbnd*dt + units.qbnd*dt2);
         *estr += e;
         if (minim_values.iprint)
           fprintf(pcmlogfile,"Bond: (%-3d) - (%-3d) %-8.3f %-8.3f %-8.3f = %-8.4f\n",it,kt,rik,bl[i]-bcorr,bk[i],e);
      }
   }
}
// =====================================
void ebond1(int natom,int nbnd,int **i12,int *use,double *x,double *y,double *z,double *bl,double *bk,double *estr,double **deb)
{
/* compute stretching energy and first derivatives */
   int i, it, kt;
   double xr, yr, zr, rik, rik2, bcorr;
   double dt, dt2, e, deddt;
   double de,dedx,dedy,dedz;
 
   *estr = 0.0F;
      for (i=0; i <= natom; i++)
      {
          deb[i][0] = 0.0;
          deb[i][1] = 0.0;
          deb[i][2] = 0.0;
      }
     
   for (i=0; i < nbnd; i++)
   {
      it = i12[i][0];
      kt = i12[i][1];
      bcorr = 0.0;
      if ( use[it] || use[kt] )
      {
         xr = x[it] - x[kt];
         yr = y[it] - y[kt];
         zr = z[it] - z[kt];
         rik2 = xr*xr + yr*yr + zr*zr;
         rik = sqrt(rik2);
         dt = rik - (bl[i] - bcorr);        
         dt2 = dt*dt;
         e = units.bndunit *bk[i]*dt2*
            (1.0 + units.cbnd*dt + units.qbnd*dt2);
         deddt = 2.0 * units.bndunit * bk[i] * dt
                 * (1.0+1.5*units.cbnd*dt+2.0*units.qbnd*dt2);
         if (rik == 0.0)
            de = 0.0;
         else
            de = deddt/rik;
         dedx = de*xr;
         dedy = de*yr;
         dedz = de*zr;
         *estr += e;
         deb[it][0] += dedx;
         deb[it][1] += dedy;
         deb[it][2] += dedz;
         deb[kt][0] -= dedx;
         deb[kt][1] -= dedy;
         deb[kt][2] -= dedz;         

      }
   }
}
// =============================================  
void ebond2(int ia,int **i12,int **iat,int **bo,double *x,double *y,double *z,double *bl,double *bk,float **hessx,float **hessy,float **hessz)
{
    int j,k,m,ibond;
    double  dt,dt2;
    double xr,yr,zr,rik,rik2,deddt,d2eddt2;
    double de,term,termx,termy,termz,d2e[3][3];
    
    for (m=0; m < MAXIAT; m++)
    {
        if (iat[ia][m]!= 0 && bo[ia][m] != 9)
        {
            ibond = find_bond(ia, iat[ia][m]);
            if (i12[ibond][0] == ia)
               k = i12[ibond][1];
            else
               k = i12[ibond][0];

            xr = x[ia] - x[k];
            yr = y[ia] - y[k];
            zr = z[ia] - z[k];
            rik2 = xr*xr + yr*yr + zr*zr;
            rik = sqrt(rik2);
            dt = rik - bl[ibond];
            dt2 = dt * dt;

            deddt = 2.0 * units.bndunit * bk[ibond] * dt
                      * (1.0+1.5*units.cbnd*dt+2.0*units.qbnd*dt2);
            d2eddt2 = 2.0 * units.bndunit * bk[ibond]
                        * (1.0+3.0*units.cbnd*dt+6.0*units.qbnd*dt2);
               
            if (rik2 == 0.0)
            {
              de = 0.0;
              term = 0.0;
            }else
            {
              de = deddt / rik;
              term = (d2eddt2-de) / rik2;
            }

            termx = term * xr;
            termy = term * yr;
            termz = term * zr;
            d2e[0][0] = termx*xr + de;
            d2e[1][0] = termx*yr;
            d2e[2][0] = termx*zr;
            d2e[0][1] = d2e[1][0];
            d2e[1][1] = termy*yr + de;
            d2e[2][1] = termy*zr;
            d2e[0][2] = d2e[2][0];
            d2e[1][2] = d2e[2][1];
            d2e[2][2] = termz*zr + de;
    
            for (j=0; j < 3; j++)
            {
               hessx[ia][j] += d2e[j][0];
               hessy[ia][j] += d2e[j][1];
               hessz[ia][j] += d2e[j][2];
               hessx[k][j] -= d2e[j][0];
               hessy[k][j] -= d2e[j][1];
               hessz[k][j] -= d2e[j][2];
            }
        }
    }  
}
