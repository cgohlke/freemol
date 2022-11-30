#define EXTERN extern

#include "pcwin.h"

#define STILL 1
#define HCT   2

EXTERN struct t_solvent {
    int type;
    double EPSin, EPSsolv;
    double doffset, p1,p2,p3,p4,p5;
    double *shct,*asolv,*rsolv,*vsolv,*gpol,*rborn;
    } solvent;

void esolv(int natom,double chrgcut,int **skip,int *use,double *charge,double *x,double *y,double *z,double *esolv);
void esolv1(int natom,double chrgcut,int **skip,int *use,double *charge,double *x,double *y,double *z,double *esolv,double **desolv,double *drb);
void esolv2(int ia,int natom,double chrgcut,double *charge,double *x,double *y,double *z,float **hessx,float **hessy,float **hessz);
void born1(int natom,int **skip,double *x,double *y,double *z,double **desolv,double *drb);
void born(int natom,int **skip, double *x,double *y,double *z);
// ===============================
void esolv(int natom,double chrgcut,int **skip,int *use, double *charge,double *x,double *y,double *z,double *esolv)
{
    int i,j;
    double e,ai,ri,rb;
    double term,probe,term1;
    double dwater, cc,f,cutoff;
    double xi,yi,zi,xr,yr,zr,rik2;
    double rb2,fgb;

    probe = 1.40;
    dwater = 78.3;
    cc = 4.80298*4.80298*14.39418;
    cutoff = chrgcut*chrgcut;

    born(natom,skip,x,y,z);
    term = 4.0*PI;
    for (i=1; i <= natom; i++)
    {
        ai = solvent.asolv[i];
        ri = solvent.rsolv[i];
        rb = solvent.rborn[i];
        if (rb > 0.001)
        {
            term1 = (ri/rb)*(ri/rb);
            e = ai*term*(ri+probe)*(ri+probe)*term1*term1*term1;
            *esolv += e;
        }

    }
// polarization term
    f = -cc*(1.0/solvent.EPSin - 1.0/solvent.EPSsolv);

    for (i=1; i <= natom; i++)
    {
      if (charge[i] != 0.0)
      {
         for(j=i; j <= natom; j++)
         {
           if (use[i] || use[j])
           {
              if (charge[j] != 0.0)
              {
                 xi = x[i];
                 yi = y[i];
                 zi = z[i];

                 xr = xi - x[j];
                 yr = yi - y[j];
                 zr = zi - z[j];
                 rik2 =  xr*xr + yr*yr + zr*zr;
                 if (rik2 < cutoff)
                 {
                    rb2 = solvent.rborn[i]*solvent.rborn[j];
                    fgb = sqrt(rik2 + rb2*exp(-0.25*rik2/rb2));
                    e = f*charge[i]*charge[j]/fgb;
                    if (i == j) e *= 0.5;
                    *esolv += e;
                }
              }
           }
         }
      }
    }    
}
// =======================
void esolv1(int natom,double chrgcut,int **skip,int *use,double *charge,double *x,double *y,double *z,double *esolv,double **desolv,double *drb)
{
    int i,j;
    double e,de,derb,ai,ri,rb;
    double term,probe,term1;
    double dwater, cc,f,cutoff;
    double xi,yi,zi,xr,yr,zr,rik2,r;
    double rb2,fgb,fgb2,expterm;
    double rbi,rbk,dedx,dedy,dedz,drbi,drbk;

    probe = 1.40;
    dwater = 78.3;
    cc = 4.80298*4.80298*14.39418;
    cutoff = chrgcut*chrgcut;
   for (i=1; i <= natom; i++)
   {
       drb[i] = 0.0;
       desolv[i][0] = 0.0;
       desolv[i][1] = 0.0;
       desolv[i][2] = 0.0;
   }

    born(natom,skip,x,y,z);
    term = 4.0*PI;
    for (i=1; i <= natom; i++)
    {
        ai = solvent.asolv[i];
        ri = solvent.rsolv[i];
        rb = solvent.rborn[i];
        if (rb > 0.01)
        {
            term1 = (ri/rb)*(ri/rb);
            e = ai*term*(ri+probe)*(ri+probe)*term1*term1*term1;
            *esolv += e;
            drb[i] -= 6.0*e/rb;
        }
    }
// polarization term
    f = -cc*(1.0/solvent.EPSin - 1.0/solvent.EPSsolv);

    for (i=1; i <= natom; i++)
    {
      if (fabs(charge[i]) > 0.001)
      {
         xi = x[i];
         yi = y[i];
         zi = z[i];
         rbi = solvent.rborn[i];
         for(j=i; j <= natom; j++)
         {
           if (use[i] || use[j])
           {
              if (fabs(charge[j]) > 0.001)
              {
                 xr = xi - x[j];
                 yr = yi - y[j];
                 zr = zi - z[j];
                 rik2 =  xr*xr + yr*yr + zr*zr;
                 if (rik2 < cutoff)
                 {
                    r = sqrt(rik2);
                    rbk = solvent.rborn[j];
                    rb2 = rbi*rbk;
                    expterm = exp(-0.25*rik2/rb2);
                    fgb2 = rik2 + rb2*expterm;
                    fgb = sqrt(fgb2);
                    e = f*charge[i]*charge[j]/fgb;
                    de = -e*(r-0.25*r*expterm)/fgb2;
                    derb = -e*expterm*(0.5+0.125*rik2/rb2)/fgb2;
                    if (i == j)
                    {
                        e *= 0.5;
                        *esolv += e;
                        drbi = derb*rbk;
                        drb[i] += drbi;
                    } else
                    {
                        *esolv += e;
                        de /= r;
                        dedx = de*xr;
                        dedy = de*yr;
                        dedz = de*zr;
                        desolv[i][0] += dedx;
                        desolv[i][1] += dedy;
                        desolv[i][2] += dedz;
                        desolv[j][0] -= dedx;
                        desolv[j][1] -= dedy;
                        desolv[j][2] -= dedz;
                        
                        drbi = derb*rbk;
                        drbk = derb*rbi;
                        drb[i] += drbi;
                        drb[j] += drbk;
                    }
                }
             }
           }
         }
      }
    }
//
    born1(natom,skip,x,y,z,desolv,drb);
}
/* ===================================================== */
void born1(int natom,int **skip,double *x,double *y,double *z,double **desolv,double *drb)
{
    int i,k;
    double xi,yi,zi,cc,rvdw;
    double xr,yr,zr;
    double de;
    double r,r2,r6;
    double p5inv,pip5;
    double gpi,vk,ratio;
    double ccf,cosq,dccf;
    double sinq,term,theta;
    double rb2,ri,rk,sk,sk2;
    double lik,lik2,lik3;
    double uik,uik2,uik3;
    double dlik,duik;
    double t1,t2,t3;
    double dedx,dedy,dedz;

    cc = 4.80298*4.80298*14.39418;
    if (solvent.type == STILL)
    {
        p5inv = 1.0/solvent.p5;
        pip5 = PI*solvent.p5;
        for (i=1; i <= natom; i++)
        {
	    skip[i][i] = i;
            xi = x[i];
            yi = y[i];
            zi = z[i];
            gpi = 2.0*solvent.rborn[i]*solvent.rborn[i]/cc;
            for (k = 1; k <= natom; k++)
            {
                if (skip[i][k] != i)
                {
                    xr = x[k] - xi;
                    yr = y[k] - yi;
                    zr = z[k] - zi;
                    vk = solvent.vsolv[k];
                    r2 = xr*xr + yr*yr + zr*zr;
                    r = sqrt(r2);
                    r6 = r2*r2*r2;
                    rvdw = solvent.rsolv[i] + solvent.rsolv[k];
                    ratio = r2/(rvdw*rvdw);
                    if (ratio > p5inv)
                    {
                        ccf = 1.0;
                        dccf = 0.0;
                    } else
                    {
                        theta = ratio*pip5;
                        cosq = cos(theta);
                        term = 0.5*(1.0-cosq);
                        ccf = term*term;
                        sinq = sin(theta);
                        dccf = 2.0*term*sinq*pip5*ratio;
                    }
                    de = drb[i]*solvent.p4*gpi*vk*(4.0*ccf-dccf)/r6;
                    dedx = de*xr;
                    dedy = de*yr;
                    dedz = de*zr;
                    desolv[i][0] += dedx;
                    desolv[i][1] += dedy;
                    desolv[i][2] += dedz;
                    desolv[k][0] -= dedx;
                    desolv[k][1] -= dedy;
                    desolv[k][2] -= dedz;
                    
                }
            }
          
        }
    } else if (solvent.type == HCT)
    {
        for (i=1; i <= natom; i++)
        {
            xi = x[i];
            yi = y[i];
            zi = z[i];
            ri = solvent.rsolv[i] + solvent.doffset;
            rb2 = solvent.rborn[i]*solvent.rborn[i];
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
                    if (ri < (r+sk) )
                    {
                        lik = 1.0/ MaxFun(ri,r-sk);
                        uik = 1.0/(r+sk);
                        lik2 = lik*lik;
                        uik2 = uik*uik;
                        lik3 = lik * lik2;
                        uik3 = uik * uik2;
                        dlik = 1.0;
                        if (ri >= (r-sk)) dlik = 0.0;
                        duik = 1.0;
                        t1 = 0.5*lik2 + 0.25*sk2*lik3/r - 0.25*((lik/r)+(lik3*r));
                        t2 = -0.5*uik2 - 0.25*sk2*(uik3/r) + 0.25*((uik/r)+uik3*r);
                        t3 = 0.125*(1.0+sk2/r2)*(lik2-uik2) + 0.25*log(uik/lik)/r2;
                        de = drb[i] * rb2 * (dlik*t1+duik*t2+t3) / r;
                        dedx = de*xr;
                        dedy = de*yr;
                        dedz = de*zr;
                        desolv[i][0] += dedx;
                        desolv[i][1] += dedy;
                        desolv[i][2] += dedz;
                        desolv[k][0] -= dedx;
                        desolv[k][1] -= dedy;
                        desolv[k][2] -= dedz;
                    }
                }
            }

        }
    }
}
/*    ============================================== */
void esolv2(int ia,int natom,double chrgcut,double *charge,double *x,double *y,double *z,float **hessx,float **hessy,float **hessz)
{
    int jj,k;
    double probe,cc,cutoff;
    double fi,fik,de,d2e;
    double d2edx,d2edy,d2edz;
    double xi,yi,zi,xr,yr,zr;
    double r,r2;
    double dwater,rb2;
    double expterm;
    double fgb,fgb2,dfgb,dfgb2,d2fgb;
    double term[3][3];

    probe = 1.40;
    dwater = 78.3;
    cc = 4.80298*4.80298*14.39418;

    if (fabs(charge[ia]) > 0.001)
       return;
        
    cutoff = chrgcut*chrgcut;
    fi = -cc*(1.0/solvent.EPSin - 1.0/solvent.EPSsolv)*charge[ia];
    xi = x[ia];
    yi = y[ia];
    zi = z[ia];

    for (k = 1; k <= natom; k++)
    {
        if (ia != k && charge[k] != 0.0)
        {
            xr = xi - x[k];
            yr = yi - y[k];
            zr = zi - z[k];
            r2 = xr*xr + yr*yr + zr*zr;
            if (r2 < cutoff)
            {
                r = sqrt(r2);
                fik = fi*charge[k];
                rb2 = solvent.rborn[ia]*solvent.rborn[k];
                expterm = exp(-0.25*r2/rb2);
                fgb2 = r2 + rb2*expterm;
                fgb = sqrt(fgb2);
                dfgb = (1.0-0.25*expterm)*r/fgb;
                dfgb2 = dfgb*dfgb;
                d2fgb = -dfgb2/fgb + dfgb/r + 0.125*(r2/rb2)*expterm/fgb;
                de = -fik*dfgb/fgb2;
                d2e = -fik*(d2fgb-2.0*dfgb2/fgb)/fgb2;
                de /= r;
                d2e = (d2e-de)/r;
                d2edx = d2e*xr;
                d2edy = d2e*yr;
                d2edz = d2e*zr;
                term[0][0] = d2edx*xr + de;
                term[1][0] = d2edx*yr;
                term[2][0] = d2edx*zr;
                term[0][1] = term[1][0];
                term[1][1] = d2edy*yr + de;
                term[2][1] = d2edy*zr;
                term[0][2] = term[2][0];
                term[1][2] = term[2][1];
                term[2][2] = d2edz*zr + de;
                for (jj=0; jj < 3; jj++)
                {
                    hessx[ia][jj] += term[jj][0];
                    hessy[ia][jj] += term[jj][1];
                    hessz[ia][jj] += term[jj][2];
                    hessx[k][jj] -= term[jj][0];
                    hessy[k][jj] -= term[jj][1];
                    hessz[k][jj] -= term[jj][2];
                }
            }
        }
    }
}

                
                
                
            
    
