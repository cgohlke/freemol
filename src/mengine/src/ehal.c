#define EXTERN extern

#ifdef KEG_OPENMP
#include <omp.h>
#endif

#include "pcwin.h"
#include "utility.h"
        
EXTERN struct t_minim_values {
        int iprint, ndc, nconst;
        float dielc;
        } minim_values;
        
//void ehal(natom,cutoffs.vdwcut,atom.x,atom.y,atom.z, atom.type,atom.use,skip,nonbond.vrad,nonbond.veps,double *energies.evdw,double *energies.e14)
//void ehal1(natom,cutoffs.vdwcut,atom.x,atom.y,atom.z, atom.type,atom.use,skip,nonbond.vrad,nonbond.veps,double *energies.evdw,double *energies.e14,
//           derivs.devdw,derivs.de14)
//void ehal2(natom,cutoffs.vdwcut,atom.x,atom.y,atom.z, atom.type,atom.use,skip,nonbond.vrad,nonbond.veps,hessian.hessx,hessian.hessy,hessian.hessz)

void ehal(int natom,int *type,int *use,double *x,double *y,double *z,double vdwcut,int **skip,double **vrad,double **veps,double *evdw,double *e14);
void ehal1(int natom,int *type,int *use,double *x,double *y,double *z,double vdwcut,int **skip,double **vrad,double **veps,double *evdw,double *e14,
           double **devdw,double **de14);
void ehal2(int iatom,int natom,int *type,int *use,double *x,double *y,double *z,double vdwcut,int **skip,double **vrad,double **veps,
           float **hessx,float **hessy,float **hessz);

static int icount = 0;

// =============================      
void ehal(int natom,int *type, int *use, double *x, double *y, double *z,double vdwcut,int **skip,double **vrad, double **veps,double *evdw,double *e14)
{
      int i,ia, ib, ita, itb, j;
      double xi,yi,zi,xr,yr,zr;
      double xk,yk,zk;
      double e,rv,eps;
      double rik,rik2,rik7;
      double sigma, kappa, kappa7;
      double rv7,rho, tau;
      double cutoff;
      
      cutoff = vdwcut*vdwcut;
      *evdw = 0.0;
      *e14  = 0.0;
      sigma = 1.12;
      kappa = 1.07;
      kappa7 = pow(kappa,7.0);

      if (minim_values.iprint)
      {
          fprintf(pcmlogfile,"\nVDW Terms - MMFF Buffered 14-7 Potential\n");
          fprintf(pcmlogfile,"     At1   At2    Rik    Radius     Eps       Evdw\n");
      }
      for (i=1; i < natom; i++)
      {

         for(j=i+1; j <= natom; j++)
         {
             if (skip[i][j] != i)
             {
               ia = i;
               ib = j;
               if (use[ia] || use[ib])
               {
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
                 if (rik2 < cutoff)
                 {
                     rv = vrad[ita][itb];
                     eps = veps[ita][itb];
                     if (skip[i][j] == -i)
                       eps /= units.v14scale;
                    rik = sqrt(rik2);
                    rv7 = pow(rv,7.0);
                    rik7 = pow(rik,7.0);
                    rho = rik7 + (sigma-1.00)*rv7;
                    tau = rik + (kappa-1.00)*rv;

                    e = eps*kappa7*(rv7/pow(tau,7.0))*(sigma*rv7/rho-2.00);
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
// ==============================
void ehal1(int natom,int *type, int *use, double *x, double *y, double *z,double vdwcut,int **skip,double **vrad, double **veps,double *evdw,double *e14,
	   double **devdw,double **de14)
{
  int i,ia, ib, j, ita, itb,skipij;
      double xi,yi,zi,xr,yr,zr;
      double xk,yk,zk;
      double e,rv,eps;
      double rik,rik2;
      double rik6,rik7;
      double sigma, kappa, kappa7;
      double rv7,rho, tau,tau7,tau8, rv14;
      double dedx,dedy,dedz,de;
      double cutoff,v14scale;
      double sum_evdw,sum_e14;
      double **ldevdw,**lde14;
      
      icount++;
      cutoff = vdwcut*vdwcut;
      v14scale = units.v14scale;
      
      *evdw = 0.0;
      *e14  = 0.0;
      sum_evdw = 0.0;
      sum_e14 = 0.0;
      ldevdw = dmatrix(0,natom+1,0,3);
      lde14 = dmatrix(0,natom+1,0,3);
      
      sigma = 1.12;
      kappa = 1.07;
      kappa7 = kappa*kappa*kappa*kappa*kappa*kappa*kappa;
      
      for (i=1; i <= natom; i++)
      {
          ldevdw[i][0] = 0.0;
          ldevdw[i][1] = 0.0;
          ldevdw[i][2] = 0.0;
          lde14[i][0] = 0.0;
          lde14[i][1] = 0.0;
          lde14[i][2] = 0.0;
      }

      /*#pragma omp parallel for private(i,j,ia,ib,ita,itb,xi,yi,zi,xk,yk,zk,xr,yr,zr,rik,rik2,rv,eps,rv7,rv14,rik6,rik7,rho,tau,tau7,tau8,e,de,dedx,dedy,dedz,skipij) \
      shared(x,y,z,use,type,skip,vrad,veps,sigma,kappa7,cutoff,natom,v14scale) reduction(+:sum_evdw,sum_e14) firstprivate(ldevdw,lde14) lastprivate(ldevdw,lde14)
      */
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
		skipij = skip[i][j];
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
                    if (skipij == -i)
                        eps /= v14scale;
                    rik = sqrt(rik2);
		    rv7 = rv*rv*rv*rv*rv*rv*rv;  //   pow(rv,7.0);
                    rv14 = rv7*rv7;
                    rik6 = rik2*rik2*rik2;
                    rik7 = rik6*rik;
                    rho = rik7 + (sigma-1.00)*rv7;
                    tau = rik + (kappa-1.00)*rv;
		    tau7 = tau*tau*tau*tau*tau*tau*tau;
       		    tau8 = tau7*tau;
                    e = eps*kappa7*(rv7/tau7)*(sigma*rv7/rho-2.00);
         
                    de = -7.00*eps*kappa7*(rv7/tau8)*(sigma*rv7/rho-2.00)
                         -7.00*eps*kappa7*sigma*rv14*rik6/(rho*rho*tau7);

                    de /= rik;
                    dedx = de*xr;
                    dedy = de*yr;
                    dedz = de*zr;
                    if (skipij == -i)
                    {
                      sum_e14 += e;
                      lde14[ia][0] += dedx;
                      lde14[ia][1] += dedy;
                      lde14[ia][2] += dedz;
                      lde14[ib][0] -= dedx;
                      lde14[ib][1] -= dedy;
                      lde14[ib][2] -= dedz;
                    } else
                    {
                      sum_evdw += e;
                      ldevdw[ia][0] += dedx;
                      ldevdw[ia][1] += dedy;
                      ldevdw[ia][2] += dedz;
                      ldevdw[ib][0] -= dedx;
                      ldevdw[ib][1] -= dedy;
                      ldevdw[ib][2] -= dedz;
                    }
		 }
              }
            }
         }
      }
      *evdw = sum_evdw;
      *e14 = sum_e14;
      for (i=1;i<= natom; i++)
	{
	  devdw[i][0] = ldevdw[i][0];
	  devdw[i][1] = ldevdw[i][1];
	  devdw[i][2] = ldevdw[i][2];
	  de14[i][0] = lde14[i][0];
	  de14[i][1] = lde14[i][1];
	  de14[i][2] = lde14[i][2];
	}
      free_dmatrix(ldevdw, 0,natom+1,0,3);
      free_dmatrix(lde14,0,natom+1,0,3);
}
// =============================
void ehal2(int iatom,int natom,int *type,int *use,double *x,double *y,double *z,double vdwcut,int **skip,double **vrad,double **veps,
          float **hessx,float **hessy,float **hessz)
{
      int j,ia, ib, k, ita, itb;
      double xi,yi,zi,xr,yr,zr;
      double xk,yk,zk;
      double rv,eps;
      double rik,rik2, rik5;
      double rik6,rik7, rik12;
      double sigma, kappa, kappa7;
      double rv7,rho, tau, rv14;
      double de,d2e;
      double d2edx,d2edy,d2edz,term[3][3];
      double cutoff;
      
      cutoff = vdwcut*vdwcut;
     
      sigma = 1.12;
      kappa = 1.07;
      kappa7 = pow(kappa,7.0);

    ita = type[iatom];

    ia = iatom;
    ita = type[ia];    
    xi = x[ia];
    yi = y[ia];
    zi = z[ia];
    skip[iatom][iatom] = iatom;
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
            rv7 = pow(rv,7.0);
            rv14 = rv7*rv7;
            rik6 = rik2*rik2*rik2;
            rik5 = rik6/rik;
            rik7 = rik6*rik;
            rik12 = rik6*rik6;
            rho = rik7 + (sigma-1.00)*rv7;
            tau = rik + (kappa-1.00)*rv;

            de = -7.00*eps*kappa7*(rv7/pow(tau,8.00))*(sigma*rv7/rho-2.00)
                 -7.00*eps*kappa7*sigma*rv14*rik6/(rho*rho*pow(tau,7.00));

            d2e = 56.0*eps*kappa7*(rv7/pow(tau,9.00))*(sigma*rv7/rho-2.00)
                 + 98.0*eps*kappa7*sigma*rv14*rik6/(rho*rho*pow(tau,8.00))
                 + 98.0*eps*kappa7*sigma*rv14*rik12/(rho*rho*rho*pow(tau,7.00))
                 - 42.0*eps*kappa7*sigma*rv14*rik5/(rho*rho*pow(tau,7.00));

           de /= rik;
           d2e = (d2e-de)/rik2;
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

