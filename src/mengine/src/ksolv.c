#define EXTERN extern

#include "pcwin.h"
#include "bonds_ff.h"
#include "angles.h"

#define STILL 1
#define HCT   2

EXTERN struct t_solvent {
    int type;
    double EPSin, EPSsolv;
    double doffset, p1,p2,p3,p4,p5;
    double *shct,*asolv,*rsolv,*vsolv,*gpol,*rborn;
    } solvent;


int find_bond(int,int);
void born(int natom,int **skip, double *x,double *y,double *z);
void ksolv(int natom,int *atomnum,int **iat,int **bo);
void message_alert(char *, char *);

// choices are gbsa models: Analytical Still and HCT
// allocate memory in get_mem

void ksolv(int natom,int *atomnum,int **iat,int **bo)
{
    int i,j,k, jji,it,kt;
    int ia,ib,ic,id;
    double ri,ri2,rk,cc;
    double r,r2,r4,rab,rbc;
    double cosine,factor;
    double h,ratio,term;

    solvent.doffset = -0.09;
    solvent.p1 = 0.073;
    solvent.p2 = 0.921;
    solvent.p3 = 6.211;
    solvent.p4 = 15.236;
    solvent.p5 = 1.254;

    if (solvent.type == HCT)
    {
        for (i=1; i <= natom; i++)
        {
            solvent.asolv[i] = 0.0054;
            solvent.shct[i] = 0.0;
            if (atomnum[i] == 1) solvent.shct[i] = 0.85;
            if (atomnum[i] == 6) solvent.shct[i] = 0.72;
            if (atomnum[i] == 7) solvent.shct[i] = 0.79;
            if (atomnum[i] == 8) solvent.shct[i] = 0.85;
            if (atomnum[i] == 9) solvent.shct[i] = 0.88;
            if (atomnum[i] == 15) solvent.shct[i] = 0.86;
            if (atomnum[i] == 16) solvent.shct[i] = 0.96;
            if (atomnum[i] == 26) solvent.shct[i] = 0.88;
        }
    } else if (solvent.type == STILL)
    {
        for (i=1; i <= natom; i++)
        {
            solvent.asolv[i] = 0.0049;
        }
    }
//  assign standard radii
    for (i=1; i <= natom; i++)
    {
        solvent.rsolv[i] = 0.0;
        if (atomnum[i] == 1)
        {
            solvent.rsolv[i] = 1.25;
            if (atomnum[iat[i][0]] == 7) solvent.rsolv[i] = 1.15;
            if (atomnum[iat[i][0]] == 8) solvent.rsolv[i] = 1.05;
        } else if (atomnum[i] == 3)
        {
            solvent.rsolv[i] = 1.432;
        } else if (atomnum[i] == 6)
        {
            solvent.rsolv[i] = 1.90;
            jji = 0;
            for (j=0; j < 4; j++)
            {
                if (iat[i][j] != 0)
                   jji++;
            }
            if (jji == 3) solvent.rsolv[i] = 1.875;
            if (jji == 2) solvent.rsolv[i] = 1.825;
        } else if (atomnum[i] == 7)
        {
            solvent.rsolv[i] = 1.7063;
            jji = 0;
            for (j=0; j < 4; j++)
            {
                if (iat[i][j] != 0)
                   jji++;
            }
            if (jji == 4) solvent.rsolv[i] = 1.625;
            if (jji == 1) solvent.rsolv[i] = 1.60;
        } else if (atomnum[i] == 8)
        {
            solvent.rsolv[i] = 1.535;
            jji = 0;
            for (j=0; j < 4; j++)
            {
                if (iat[i][j] != 0)
                   jji++;
            }
            if (jji == 1) solvent.rsolv[i] = 1.48;
        } else if (atomnum[i] == 9)
        {
            solvent.rsolv[i] = 1.47;
        } else if (atomnum[i] == 10)
        {
            solvent.rsolv[i] = 1.39;
        } else if (atomnum[i] == 11)
        {
            solvent.rsolv[i] = 1.992;
        } else if (atomnum[i] == 12)
        {
            solvent.rsolv[i] = 1.70;
        } else if (atomnum[i] == 14)
        {
            solvent.rsolv[i] = 1.80;
        } else if (atomnum[i] == 15)
        {
            solvent.rsolv[i] = 1.87;
        } else if (atomnum[i] == 16)
        {
            solvent.rsolv[i] = 1.775;
        } else if (atomnum[i] == 17)
        {
            solvent.rsolv[i] = 1.735;
        } else if (atomnum[i] == 18)
        {
            solvent.rsolv[i] = 1.70;
        } else if (atomnum[i] == 19)
        {
            solvent.rsolv[i] = 2.123;
        } else if (atomnum[i] == 20)
        {
            solvent.rsolv[i] = 1.817;
        } else if (atomnum[i] == 35)
        {
            solvent.rsolv[i] = 1.90;
        } else if (atomnum[i] == 36)
        {
            solvent.rsolv[i] = 1.812;
        } else if (atomnum[i] == 37)
        {
            solvent.rsolv[i] = 2.26;
        } else if (atomnum[i] == 53)
        {
            solvent.rsolv[i] = 2.10;
        } else if (atomnum[i] == 54)
        {
            solvent.rsolv[i] = 1.967;
        } else if (atomnum[i] == 55)
        {
            solvent.rsolv[i] = 2.507;
        } else if (atomnum[i] == 56)
        {
            solvent.rsolv[i] = 2.188;
        }else
        {
            message_alert("Error - Atom not available for GB/SA","Error");
            break;
        }
    }
//  atomic volumes for Still
    if (solvent.type == STILL)
    {
        for (i=1; i <= natom; i++)
        {
               solvent.vsolv[i] = (4.0*PI/3.0)*solvent.rsolv[i]*solvent.rsolv[i]*solvent.rsolv[i];
               ri = solvent.rsolv[i];
               ri2 = ri*ri;
               for (j=0; j < MAXIAT; j++)
               {
                  if (iat[i][j] != 0 && bo[i][j] != 9)
                  {
                    k = find_bond(i,iat[i][j]);
                    rk = solvent.rsolv[iat[i][j]];
                    r = 1.01 * bonds_ff.bl[k];
                    ratio = (rk*rk - ri2 - r*r)/(2.0*ri*r);
                    h = ri*(1.0+ratio);
                    term = (PI/3.0)*h*h*(3.0*ri-h);
                    solvent.vsolv[i] -= term;
                  }
               }
 
        }

// 1,2 and 1,3 polarization
        cc = 4.80298*4.80298*14.39418;

        for (i=1; i <= natom; i++)
           solvent.gpol[i] = -0.5*cc/(solvent.rsolv[i]+solvent.doffset+solvent.p1);
        
        
       for (i=0; i < bonds_ff.nbnd; i++)
       {
          it = bonds_ff.i12[i][0];
          kt = bonds_ff.i12[i][1];
          r = bonds_ff.bl[i];
          r4 = r*r*r*r;
          solvent.gpol[it] += solvent.p2*solvent.vsolv[kt]/r4;
          solvent.gpol[kt] += solvent.p2*solvent.vsolv[it]/r4;
  
       }

       for (i=0; i < angles.nang; i++)
       {
           ia = angles.i13[i][0];
           ib = angles.i13[i][1];
           ic = angles.i13[i][2];
             factor = 1.0;
             for (j=0; j < MAXIAT; j++)
             {
               if (iat[ia][j] != 0)
               {
                   id = iat[ia][j];
                   if (id == ic)
                      factor = 0.0;
                   else if (id != ib)
                   {
                       for (k=0; k < MAXIAT; k++)
                       {
                           if (iat[ic][k] != 0)
                           {
                               if (iat[ic][k] == id)
                                   factor = 0.5;
                           }
                       }
                   }
               }
             }
     
             id = find_bond(ia,ib);
             rab = bonds_ff.bl[id];
             id = find_bond(ib,ic);
             rbc = bonds_ff.bl[id];
             cosine = cos(angles.anat[i]/radian);
             r2 = rab*rab + rbc*rbc - 2.0*rab*rbc*cosine;
             r4 = r2*r2;
             solvent.gpol[ia] += factor*solvent.p3*solvent.vsolv[ic]/r4;
             solvent.gpol[ic] += factor*solvent.p3*solvent.vsolv[ia]/r4;
           
       }
    }
}
/* =================================================== */
void born(int natom,int **skip, double *x,double *y,double *z)
{
    int i,k;
    double cc;
    double ratio;
    double xi,yi,zi,ri;
    double rk,sk,sk2,sum;
    double lik,lik2,uik,uik2;
    double xr,yr,zr,rvdw;
    double r,r2,r4;
    double gpi,pip5,p5inv;
    double theta,term,ccf;

    cc = 4.80298*4.80298*14.39418;
    if (solvent.type == STILL)
    {
        p5inv = 1.0/solvent.p5;
        pip5 = PI*solvent.p5;
        for (i=1; i <= natom; i++)
        {
               xi = x[i];
               yi = y[i];
               zi = z[i];
               gpi = solvent.gpol[i];
               skip[i][i] = i;
               for (k = 1; k <= natom; k++)
               {
                  if (skip[i][k] != i)
                  {
                    xr = x[k] - xi;
                    yr = y[k] - yi;
                    zr = z[k] - zi;
                    r2 = xr*xr + yr*yr + zr*zr;
                    r4 = r2*r2;
                    rvdw = solvent.rsolv[i] + solvent.rsolv[k];
                    ratio = r2/(rvdw*rvdw);
                    if (ratio > p5inv)
                    {
                        ccf = 1.0;
                    } else
                    {
                        theta = ratio*pip5;
                        term = 0.5*(1.0-cos(theta));
                        ccf = term*term;
                    }
                    gpi += solvent.p4*ccf*solvent.vsolv[k]/r4;
                  }
               }
               solvent.rborn[i] = -0.5*cc/gpi;
            
        }
    }else if (solvent.type == HCT)
    {
        for (i=1; i <= natom; i++)
        {
               xi = x[i];
               yi = y[i];
               zi = z[i];
               ri = solvent.rsolv[i] + solvent.doffset;
               sum = 1.0/ri;
               for (k=1; k <= natom; k++)
               {
                   if (i != k)
                   {
                      xr = x[k] - xi;
                      yr = y[k] - yi;
                      zr = z[k] - zi;
                      r2 = xr*xr + yr*yr + zr*zr;
                      r = sqrt(r2);
                      rk = solvent.rsolv[k]+solvent.doffset;
                      sk = rk * solvent.shct[k];
                      sk2 = sk*sk;
                      if (ri < (r+sk))
                      {
                        lik = 1.0/ MaxFun(ri,r-sk);
                        uik = 1.0/(r+sk);
                        lik2 = lik*lik;
                        uik2 = uik*uik;
                        term = lik - uik + 0.25*r*(uik2-lik2) + (0.5/r)*log(uik/lik)
                              + (0.25*sk2/r)*(lik2-uik2);
                        sum -= 0.5*term;
                      }
                   }
               }
               theta = 1.0/sum;
               solvent.rborn[i] = MaxFun(ri,theta);
 
        }
    }
}
                   
            

              
      
