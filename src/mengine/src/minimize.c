#define EXTERN extern

#include "pcwin.h"
#include "pcmod.h"
#include "utility.h"
#include "job_control.h"
#include <sys/timeb.h>
#include <sys/time.h>
#include <time.h>

// mode
#define NONE      0
#define  Failure     6

double energy(void);
void tncg(double,int,int,int *,double *,double *, double, double (*)(double *, double *),
   void (*)(int, double *,double *,int *, int *, int *, double *));
void mqn(double,int , int, int *,double *,double *, double *, double (*)(double *,double *) );
void search(int,double *,double *,double *,double *,double,double *,int *,
 double (*)(double *,double *),int *);
void newton2(int, double *,double *,int *, int *, int *, double *);
void hessian(int, double *, int *, int *, int *,double *);
void gradient(void);
double minimiz1(double *, double *);
void minimize(int natom,int *use,double *x,double *y,double *z);
double get_total_energy(void);
double get_total_deriv_x(int i);
double get_total_deriv_y(int i);
double get_total_deriv_z(int i);
void write_sdf(int);

struct t_minvar{
     double cappa, stpmin, stpmax, angmax;
     int   intmax;
     }  minvar;
     
EXTERN struct t_minim_control {
        int type, method, field, added_const;
        char added_path[256],added_name[256];
        } minim_control;
        
EXTERN struct t_minim_values {
        int iprint, ndc, nconst;
        float dielc;
        } minim_values;
    
double scale2;

// =====================================
void minimize(int natom,int *use,double *x,double *y,double *z)
{
   int i,nvar, iter, icount;
   double minimum,grdmin;
   //   double minimiz1();
   double etot;
   double *xx;  //  xx[maxvar];
   int method=0, maxvar;
   double start_time;
   
   maxvar = 3*natom;
   xx = dvector(0,maxvar);
      
   grdmin = 1.0;
   scale2 = 12.0; // bfgs  12.0

   if (job_control.monitor)
     {
#ifdef WIN32
   struct _timeb timebuffer;
   _ftime( &timebuffer );
   start_time = (timebuffer.time+(timebuffer.millitm/((double)1000.0)));
#else
   struct timeval tv;
   gettimeofday(&tv,NULL);
   start_time = tv.tv_sec + (tv.tv_usec/(double)1000000.0);
#endif
     }
   else 
     start_time = 0;
   if (minim_control.method == 1 || minim_control.method == 3 || minim_control.method == 4)
   {
     nvar = 0;
     icount = 0;
     for (i=1; i <= natom; i++)
     {
       if (use[i])
       {
         xx[nvar] = x[i]*scale2;
         nvar++;
         xx[nvar] = y[i]*scale2;
         nvar++;
         xx[nvar] = z[i]*scale2;
         nvar++;
       }
     }

     method = 1;
     grdmin = 0.5;
     if (minim_control.method == 1)
        grdmin = 0.1;
     else
        grdmin = 0.5;
       
     mqn(start_time,nvar,method, &iter,xx, &minimum, &grdmin, minimiz1 );    
   }

   if (grdmin > 1.00)
   {
       free_dvector(xx, 0, maxvar);
       return;
   }
      if (minim_control.method == 2 || minim_control.method == 3 || minim_control.method == 4)
      {
        scale2 = 1.0;   // tcng
        grdmin = 0.0001;
        nvar = 0;
        icount = 0;
        for (i=1; i <= natom; i++)
        {
          if (use[i])
          {
            xx[nvar] = x[i]*scale2;
            nvar++;
            xx[nvar] = y[i]*scale2;
            nvar++;
            xx[nvar] = z[i]*scale2;
            nvar++;
          }
        }
        tncg(start_time, nvar,method,&iter, xx, &minimum, grdmin,minimiz1, newton2); 
       
   
        nvar = 0;
        for (i=1; i <= natom; i++)
        {
          if (use[i])
          {
            x[i] =  xx[nvar]/scale2;
            nvar++;
            y[i] =  xx[nvar]/scale2;
            nvar++;
            z[i] =  xx[nvar]/scale2;
            nvar++;
          }
        }
      }
   
   if (minimum < -1000.0)
   {
    
       etot = energy();    
       free_dvector(xx, 0, maxvar);
       return;      
   }   
   free_dvector(xx, 0, maxvar);
   if (job_control.monitor)
     job_control.monitor = FALSE;
}
// ======================================
void mqn(double start_time,int nvar,int method,int *iter, double *x,  double *minimum, double *grdmin, double (*fgvalue) (double *,double *))
{
      int i,ncalls,nerror;
      int niter,period,nstart;
      double fast,slow,epsln,d1temp,d2temp;
      double f,f_old,f_new,f_move;
      double rms,beta,x_move,g_norm,g_rms;
      double gg,gg_old;
      double sg,dg,sd,dd,angle;
      double *g, *p;                         //   g[maxvar];
      double *x_old, *g_old;             //   x_old[maxvar],g_old[maxvar];
      double  *s, *d;                 //   p[maxvar],s[maxvar],d[maxvar];
      double fctmin;
      int restart, terminate;
      int maxiter, nextiter,status, maxvar;
      double end_time,old_time,ttime;

      old_time = start_time;
      maxvar = 3*natom;
      x_old = dvector(0,maxvar);
      g_old = dvector(0,maxvar);
      g = dvector(0,maxvar);
      p = dvector(0,maxvar);
      s = dvector(0,maxvar);
      d = dvector(0,maxvar);
      
      ncalls = 0;
      rms = sqrt((float)nvar)/ sqrt(3.0);
      restart = TRUE;
      terminate = FALSE;
      status = 0;
      nerror = 0;

      fctmin = -10000.0;
      maxiter = 1000;
      nextiter = 1;
      fast = 0.5;
      slow = 0.0;
      epsln = 1.0e-16;
      if (nvar > 200)
         period = nvar;
      else
         period = 200;
      minvar.cappa = .1;
      minvar.stpmin = 1.0e-20;
      minvar.stpmax = 5.0;
      minvar.angmax = 100.0;
      minvar.intmax = 5;

      niter = nextiter -1;
      maxiter = niter + maxiter;
      ncalls = ncalls + 1;
      f = fgvalue(x, g); // get function and first deriv at original point
      g_norm = 0.0;
      for (i=0; i < nvar; i++)
      {
          x_old[i] = x[i];
          g_old[i] = g[i];
          g_norm += g[i]*g[i];
      }
      g_norm = sqrt(g_norm);
      f_move = 0.5*minvar.stpmax*g_norm;
      g_rms = g_norm*scale2/rms;


     if (niter > maxiter)
        terminate = TRUE;
     if (f < fctmin)
        terminate = TRUE;
     if (g_rms < *grdmin)
         terminate = TRUE;

     while ( ! terminate)
     {
         niter++;
         status = 0;   


         if (restart || method == 0)
         {
             for (i=0; i < nvar; i++)
                 p[i] = -g[i];
             nstart = niter;
             restart = FALSE;
         } else if (method == 1) // BFGS method
         {
             sg = 0.0;
             dg = 0.0;
             dd = 0.0;
             sd = 0.0;
             for (i=0; i < nvar; i++)
             {
               sg += s[i]*g[i];
               dg += d[i]*g[i];
               dd += d[i]*d[i];
               sd += s[i]*d[i];
             }
             for (i=0; i < nvar; i++)
             {
                 d1temp = (d[i]*sg + s[i]*dg)/sd;
                 d2temp = (1.0+dd/sd)*(s[i]*sg/sd);
                 p[i] = -g[i] + d1temp - d2temp;
             }
         } else if (method == 2)  // Fletcher Reeves
         {
             gg = 0.0;
             gg_old = 0.0;
             for (i=0; i < nvar; i++)
             {
                 gg += g[i]*g[i];
                 gg_old += g_old[i]*g_old[i];
             }
             beta = gg/gg_old;
             for (i=0; i < nvar; i++)
                p[i] = -g[i] + beta*p[i];
         } else if (method == 3)  // Polak Ribere
         {
             dg = 0.0;
             gg_old = 0.0;
             for (i=0; i < nvar; i++)
             {
                 dg += d[i]*g[i];
                 gg_old += g_old[i]*g_old[i];
             }
             beta = dg/gg_old;
             for (i=0; i < nvar; i++)
                p[i] = -g[i] + beta*p[i];
        }

         // do a line search
         f_old = f;
         search(nvar,&f,g,x,p,f_move,&angle,&ncalls,fgvalue,&status);
         if (status == Failure)
         {
             g_rms = 1000.0;
             terminate = TRUE;
             goto L_DONE;
         }
         f_new = f;
         
         f_move = f_old - f_new;
         x_move = 0.0;
         g_norm = 0.0;
         for (i=0; i < nvar; i++)
         {
            s[i] = x[i] - x_old[i];
            d[i] = g[i] - g_old[i];
            x_move += s[i]*s[i];
            g_norm += g[i]*g[i];
            x_old[i] = x[i];
            g_old[i] = g[i];
         }
         x_move = sqrt(x_move) / (scale2 * rms);
         g_norm = sqrt(g_norm);
         g_rms = g_norm * scale2/rms;
         
// function increase
         if (f_move <= 0.0)
         {
 //           status = Increase;
            nerror = nerror + 1;
            if (nerror == 3)
               terminate = TRUE;
            else
               restart = TRUE;

            for(i=0; i < nvar; i++)
            {
               x[i] = x_old[i];
               g[i] = g_old[i];
            }
         }
         if (x_move < epsln)
         {
             nerror++;
             if (nerror > 3)
               terminate = TRUE;
             else
               restart = TRUE;
         }
// normal termination
         if (f < fctmin)
         {
 //           status = SmallFct;
            terminate = TRUE;
         }
         if (g_rms < *grdmin)
         {
//            status = SmallGrad;
            nerror++;
            if (nerror > 1)
             terminate = TRUE;
            else
             restart = TRUE;
         }

	 //	 if (niter%5)
	 //  printf("mqn: Energy %f grad %f\n",f_new, g_rms);
	 // check monitor
	 if (job_control.monitor)
	   {
#ifdef WIN32
        struct _timeb timebuffer;
        _ftime( &timebuffer );
        end_time = (timebuffer.time+(timebuffer.millitm/((double)1000.0)));
#else
        struct timeval tv;
        gettimeofday(&tv,NULL);
        end_time = tv.tv_sec + (tv.tv_usec/(double)1000000.0);
#endif
	     ttime = (end_time-old_time);
	     if ( ttime > 0 &&  (ttime > job_control.interval) )
	       {
		 //	        fprintf(pcmlogfile,"Times: %f %f %f\n",start_time,end_time,ttime);
		 old_time = end_time;
	         write_sdf(0);
	       }
	   }
	 // check abort signal

     }
L_DONE:
     *minimum = f;
     *grdmin = g_rms;
     *iter = niter;
     free_dvector(x_old ,0,maxvar);
     free_dvector(g_old,0,maxvar);
     free_dvector( g ,0,maxvar);
     free_dvector( p ,0,maxvar);
     free_dvector( s ,0,maxvar);
     free_dvector( d ,0,maxvar);             
}
// ==============================
double minimiz1(double *xx, double *g)
{
    int i,nvar;
    double e_min;

    nvar = 0;
    for (i = 1; i <= natom; i++)
    {
          if (atom.use[i])
          {
            atom.x[i] = xx[nvar]/scale2;
            nvar++;
            atom.y[i] = xx[nvar]/scale2;
            nvar++;
            atom.z[i] = xx[nvar]/scale2;
            nvar++;
          }
    }

    gradient();
    e_min = get_total_energy();

    nvar = 0;
    for (i=1; i <= natom; i++)
    {
        if (atom.use[i])
        {
            xx[nvar] = atom.x[i]*scale2;
            g[nvar] = get_total_deriv_x(i)/scale2;
            nvar++;
            xx[nvar] = atom.y[i]*scale2;
            g[nvar] = get_total_deriv_y(i)/scale2;
            nvar++;
            xx[nvar] = atom.z[i]*scale2;
            g[nvar] = get_total_deriv_z(i)/scale2;
            nvar++;
        }
    }
    return(e_min);
}
// =====================
void newton2(int mode, double *xx,double *h,int *hinit,
               int *hstop, int *hindex, double *hdiag)
{
    int i,j,k,nvar, maxvar, nuse, maxhess;
    int *hvar, *huse;   //  hvar[maxvar],huse[maxvar];

    if (mode == NONE)
       return;

    maxvar = 3*natom;
    hvar = ivector(0,maxvar);
    huse = ivector(0,maxvar);
    
    nvar = 0;
    nuse = TRUE;
    for (i=1; i <= natom; i++)
    {
        if (atom.use[i])
        {
            atom.x[i] = xx[nvar];
            nvar++;
            atom.y[i] = xx[nvar];
            nvar++;
            atom.z[i] = xx[nvar];
            nvar++;
        } else
           nuse = FALSE;
    }

    if (natom < 300)
      maxhess = (3*natom*(3*natom-1))/2;
    else if (natom < 800)
      maxhess = (3*natom*(3*natom-1))/3;
    else
      maxhess = (3*natom*(3*natom-1))/20;

    hessian(maxhess, h,hinit,hstop,hindex,hdiag);

    nvar = 0;
    if (nuse == FALSE)
    {
       for (i=1; i <= natom; i++)
       {
           k = 3*(i-1);
           if (atom.use[i])
           {
              for (j=0; j < 3; j++)
              {
                hvar[nvar] = j+k;
                huse[j+k] = nvar;
                nvar++;
              }
           } else
           {
              for (j=0; j < 3; j++)
                 huse[j+k] = 0;
           }
       }
       for (i=0; i < nvar; i++)
       {
          k = hvar[i];
          hinit[i] = hinit[k];
          hstop[i] = hstop[k];
          hdiag[i] = hdiag[k];
          for (j=hinit[i]; j < hstop[i]; j++)
            hindex[j] = huse[hindex[j]];
       }
    }
//
    nvar = 0;
    for (i=1; i <= natom; i++)
    {
        if (atom.use[i])
        {
            xx[nvar] = atom.x[i];
            nvar++;
            xx[nvar] = atom.y[i];
            nvar++;
            xx[nvar] = atom.z[i];
            nvar++;
        }
    }         
    free_ivector(hvar ,0,maxvar);
    free_ivector(huse ,0,maxvar);

}
            
            
            
