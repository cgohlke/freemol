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
#define NEWTON    1
#define TNCG      2
#define DTNCG     3
#define AUTO      4
// method
#define DIAG      5
#define SSOR      6
#define ICCG      7
#define BLOCK     8
#define FULL      9
// search results
#define  Success     0
#define  ReSearch    1
#define  WideAngle   2
#define  BadIntpln   3
#define  IntplnErr   4
#define  blank       5

void outminstat(int , double ,double );
void tncg(double,int,int,int *,double *,double *, double, double (*)(double *,double *),
          double (*)(int,double *,double *,int *, int *, int *, double *));
double fgvalue(double *, double *);
void search(int,double *,double *,double *,double *,double,double *,int *,double (*)(double *,double *),int *);
void solve(int,int,int,int, double *,double *,double *,double *,int *,
              int *,int *,double *,int *,int *,int *, double (*)(double *, double *));
void hmatrix(int,double *,double *,int *, int *, int *, double *);
void write_sdf(int);

EXTERN struct t_minvar{
     double cappa, stpmin, stpmax, angmax;
     int   intmax;
     }  minvar;
     
double hesscut;

void tncg(double start_time,int nvar,int imethod,int *iter,double *x,double *minimum,
          double grdmin, double (*fgvalue) (double *,double *),
          double (*hmatrix)(int,double *,double *,int *, int *, int *, double *))
{
    int i,  iter_tn, iter_cg, fg_call, newhess;
    int *h_init, *h_stop, *h_index;         // h_init[maxvar], h_stop[maxvar], h_index[MAXHESS];
    double *x_old, *g, *p, *h_diag;          // x_old[maxvar], g[maxvar], p[maxvar], h_diag[maxvar];
    double  *h;                   // h[MAXHESS];
    double f, angle, rms, x_move, f_move, f_old, g_norm, g_rms;
    int   done, negtest, nextiter, maxiter;
    double fctmin;
    int mode, hmode, status, maxvar, maxhess, method;
    int maxerr,nerr;
    double end_time,old_time,ttime;

    old_time = start_time;
    maxerr = 3;
    nerr = 0;
    maxvar = 3*natom;
   if (natom < 500)
      maxhess = (3*natom*(3*natom-1))/2;
   else if (natom < 800)
      maxhess = (3*natom*(3*natom-1))/3;
   else
      maxhess = (3*natom*(3*natom-1))/20;


    h_init = ivector(0,maxvar);
    h_stop = ivector(0,maxvar);
    h_index = ivector(0,maxhess);
    x_old = dvector(0,maxvar);
    g = dvector(0,maxvar);
    p = dvector(0,maxvar);
    h_diag = dvector(0,maxvar);
    h = dvector(0,maxhess);
    
    rms = sqrt((float) nvar)/ sqrt(3.0);

    fctmin = -1000000.0;
    nextiter = 1;
    newhess = 1;
    done = FALSE;
    minvar.cappa = 0.1;
    minvar.stpmin = 1.0e-20;
    minvar.stpmax = 1.5;
    minvar.angmax = 180.0;
    minvar.intmax = 8;
    negtest = TRUE;
    maxiter = 1000;
    
    iter_tn = nextiter - 1;
    maxiter += iter_tn;


    iter_cg = 0;
    fg_call = 1;
    f = fgvalue(x,g);
    f_old = f;
    g_norm = 0.0;
    for (i=0; i < nvar; i++)
    {
         x_old[i] = x[i];
         g_norm += g[i]*g[i];
    }
    g_norm = sqrt(g_norm);
    f_move = 0.5 * minvar.stpmax * g_norm;
    g_rms = g_norm / rms;

    done = FALSE;

    if (g_rms <= grdmin)
    {
        done = TRUE;
        *minimum = f;
    } else if ( f <= fctmin)
    {
        done = TRUE;
        *minimum = f;
    } else if (iter_tn >= maxiter)
    {
        done = TRUE;
        *minimum = f;
    }
    while ( ! done)
    {
        iter_tn++;


//   mode = tcng   method = diag
       if ( g_rms > 3.0)
        mode = TNCG;
       else
        mode = DTNCG;
                   
        hesscut = 0.0;
        if (g_rms > 10.0)
        {
            method = DIAG;
            hesscut = 1.0;
        } else if (nvar < 10)
        {
            method = SSOR;
            hesscut = 1.0;
        }else if (g_rms < 1.0)
        {
            method = ICCG;
            hesscut = .001*nvar;
            if (hesscut > 0.1) hesscut = 0.1;
        } else
        {
            method = ICCG;
            hesscut = 0.001*nvar;
            if (hesscut > 1.0) hesscut = 1.0;
        }
        hmode = FULL;
        
        if ( (iter_tn%newhess) != 0) hmode = NONE;
        if (mode == DTNCG && method == NONE) hmode = NONE;
        if (mode == DTNCG && method == DIAG) hmode = DIAG;

        hmatrix(hmode, x, h, h_init, h_stop, h_index, h_diag);

        solve(mode,method,negtest,nvar, p,x,g,h,h_init,
              h_stop,h_index,h_diag,&iter_tn,&iter_cg,&fg_call,fgvalue);
        search(nvar,&f,g,x,p,f_move,&angle,&fg_call,fgvalue, &status);    
        
         if (status != Success)
         {
                nerr++;
                if (nerr > maxerr)
                {
		  //                    fprintf(pcmlogfile,"Error in tncg :: bad search\n");
                    done = TRUE;
                }
         }
         
         f_move = f_old - f;
         f_old = f;
         x_move = 0.0;
         g_norm = 0.0;
         for(i=0; i < nvar; i++)
         {
            x_move += (x[i]-x_old[i])*(x[i]-x_old[i]);
            x_old[i] = x[i];
            g_norm += g[i]*g[i];
         }
  
         x_move = sqrt(x_move);
         x_move = x_move / rms;
         g_norm = sqrt(g_norm);
         g_rms = g_norm / rms;

         if (iter_tn >= maxiter)
            done = TRUE;

         if (f_move <= 0.0000001)
         {
            done = TRUE;
         }

         if (g_rms <= grdmin)
            done = TRUE;
         else if (f <= fctmin)
            done = TRUE;
            
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
		 // fprintf(pcmlogfile,"Times: %f %f %f\n",start_time,end_time,ttime);
		 old_time = end_time;
	         write_sdf(0);
	       }
	   }
              
         if (done)
         {
             *minimum = f;
             *iter = iter_tn;
             free_ivector(h_init ,0,maxvar);
             free_ivector(h_stop ,0,maxvar);
             free_ivector(h_index ,0,maxhess);
             free_dvector(x_old ,0,maxvar);
             free_dvector(g ,0,maxvar);
             free_dvector(p ,0,maxvar);
             free_dvector(h_diag ,0,maxvar);
             free_dvector(h ,0,maxhess);

             return;

         }       
    }
    free_ivector(h_init ,0,maxvar);
    free_ivector(h_stop ,0,maxvar);
    free_ivector(h_index ,0,maxhess);
    free_dvector(x_old ,0,maxvar);
    free_dvector(g ,0,maxvar);
    free_dvector(p ,0,maxvar);
    free_dvector(h_diag ,0,maxvar);
    free_dvector(h ,0,maxhess);

    return;   
}

