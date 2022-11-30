#define EXTERN extern

#include "pcwin.h"


EXTERN struct t_minim_values {
        int iprint, ndc, nconst;
        float dielc;
        } minim_values;

void estrbnd(int nstrbnd,int **,int **,int *use,double *x,double *y,double *z,float *anat,double *bl,double *bk,float *ksb1,
             float *ksb2,double *estrbnd);
void estrbnd1(int natom,int nstrbnd,int **,int **,int *use,double *x,double *y,double *z,float *anat,double *bl,double *bk,float *ksb1,
	      float *ksb2,double *estrbnd,double **destbn);
void estrbnd2(int iatom,int nstrbnd,int **,int **,int *use,double *x,double *y,double *z,float *anat,double *bl,double *bk,float *ksb1,
	      float *ksb2,float **hessx,float **hessy,float **hessz);
// ===================================       
void estrbnd(int nstrbnd,int **isb,int **i13,int *use,double *x,double *y,double *z,float *anat,double *bl,double *bk,
          float *ksb1,float *ksb2,double *estrbnd)
{
      int i,ij,j,k,ia,ib,ic;
      double e,dr,dt,angle,dot,cosine,dr1;
      double xia,yia,zia,xib,yib,zib;
      double xic,yic,zic;
      double rab,xab,yab,zab;
      double rcb,xcb,ycb,zcb;

      *estrbnd = 0;
      if (minim_values.iprint && nstrbnd > 0)
      {
          fprintf(pcmlogfile,"\nStretch-Bend Terms\n");
          fprintf(pcmlogfile,"          At1     At2    At3      Angle     Anat   StbnCst1   StbnCst2  Estb\n");
      }
      
      for (i=0; i < nstrbnd; i++)
      {
          ij = isb[i][0];
          ia = i13[ij][0];
          ib = i13[ij][1];
          ic = i13[ij][2];
          if (use[ia] || use[ib] || use[ic])
          {
            xia = x[ia];
            yia = y[ia];
            zia = z[ia];
            xib = x[ib];
            yib = y[ib];
            zib = z[ib];
            xic = x[ic];
            yic = y[ic];
            zic = z[ic];

            xab = xia - xib;
            yab = yia - yib;
            zab = zia - zib;
            rab = sqrt(xab*xab + yab*yab + zab*zab);
            xcb = xic - xib;
            ycb = yic - yib;
            zcb = zic - zib;
            rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb);
            if (rab*rcb != 0.0)
            {
               dot = xab*xcb + yab*ycb + zab*zcb;
               cosine = dot / (rab*rcb);
               if (cosine < -1.0)
                 cosine = -1.0;
               if (cosine > 1.0)
                 cosine = 1.0;
               angle = radian * acos(cosine);
               dt = angle - anat[ij];
               j = isb[i][1];
               k = isb[i][2];
               dr = 0.0;
               dr1 = 0.0;
               if (j > -1) dr = (rab - bl[j])*ksb1[i];
               if (k > -1) dr1 = (rcb - bl[k])*ksb2[i];
               e = units.stbnunit*dt*(dr+dr1);
               *estrbnd += e;

               if (minim_values.iprint)
		 {
		   fprintf(pcmlogfile,"StrBnd: (%-3d)-(%-3d)-(%-3d) : %-8.3f %-8.3f %-8.3f %-8.3f = %-8.4f\n",
			   ia, ib, ic,angle,dt, rab-bl[j],ksb1[i], units.stbnunit*dt*dr);
		   fprintf(pcmlogfile,"StrBnd: (%-3d)-(%-3d)-(%-3d) : %-8.3f %-8.3f %-8.3f %-8.3f = %-8.4f\n",
			   ia, ib, ic,angle,dt, rcb-bl[k],ksb2[i], units.stbnunit*dt*dr1);
		 }
            }
          }
      }
}
// =================
void estrbnd1(int natom,int nstrbnd,int **isb,int **i13,int *use,double *x,double *y,double *z,float *anat,double *bl,double *bk,float *ksb1,
	      float *ksb2,double *estrbnd,double **destb)
{
      int i,ij,j,k,ia,ib,ic;
      double e,dr,dr1,dt,angle,dot,cosine;
      double xia,yia,zia,xib,yib,zib;
      double xic,yic,zic;
      double rab,xab,yab,zab;
      double rcb,xcb,ycb,zcb;
      double xp,yp,zp,rp;
      double dedxia,dedyia,dedzia;
      double dedxib,dedyib,dedzib;
      double dedxic,dedyic,dedzic;
      double term,terma,termc, rab2, rcb2;
      double ksb11, ksb21;
      double dtemp;

      *estrbnd = 0;
      for (i=0; i <= natom; i++)
      {
          destb[i][0] = 0.0;
          destb[i][1] = 0.0;
          destb[i][2] = 0.0;
      }

      for (i=0; i < nstrbnd; i++)
      {
          ij = isb[i][0];
          ia = i13[ij][0];
          ib = i13[ij][1];
          ic = i13[ij][2];
          if (use[ia] || use[ib] || use[ic])
          {
            xia = x[ia];
            yia = y[ia];
            zia = z[ia];
            xib = x[ib];
            yib = y[ib];
            zib = z[ib];
            xic = x[ic];
            yic = y[ic];
            zic = z[ic];

            xab = xia - xib;
            yab = yia - yib;
            zab = zia - zib;
            rab2 = xab*xab + yab*yab + zab*zab;
            rab = sqrt(xab*xab + yab*yab + zab*zab);
            xcb = xic - xib;
            ycb = yic - yib;
            zcb = zic - zib;
            rcb2 = xcb*xcb + ycb*ycb + zcb*zcb;
            rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb);
            xp = ycb*zab - zcb*yab;
            yp = zcb*xab - xcb*zab;
            zp = xcb*yab - ycb*xab;
            rp = sqrt(xp*xp + yp*yp + zp*zp);
            if (rp != 0.0)
            {
               dot = xab*xcb + yab*ycb + zab*zcb;
               cosine = dot / (rab*rcb);
               if (cosine < -1.0)
                 cosine = -1.0;
               if (cosine > 1.0)
                 cosine = 1.0;
               angle = radian * acos(cosine);
               dt = angle - anat[ij];

               dtemp = dot*dot/(rab2*rcb2);
               if (fabs(dtemp) >= 1.0)
                  dtemp = 0.9999*dtemp/fabs(dtemp);
               termc = sqrt(1.00 - dtemp);
               terma = rab*rcb;
               j = isb[i][1];
               k = isb[i][2];
               term = -units.stbnunit*radian; 
               if (j < 0)
               {
                  dr = 0.0;
                  ksb11 = 0.0;
               } else
               {
                  dr = ksb1[i]*(rab - bl[j]);
                  ksb11 = ksb1[i];
               }
               if (k < 0)
               {
                  dr1 = 0.0;
                  ksb21 = 0.0;
               }  else
               {
                  dr1 = ksb2[i]*(rcb - bl[k]);
                  ksb21 = ksb2[i];
               }
               e = units.stbnunit*dt*(dr+dr1);
               dedxia = (term *(xcb/terma - (xab*dot)/(rab2*terma))*(dr+dr1))/ termc + (units.stbnunit*dt*ksb11*xab)/rab;
               dedyia = (term *(ycb/terma - (yab*dot)/(rab2*terma))*(dr+dr1))/ termc + (units.stbnunit*dt*ksb11*yab)/rab;
               dedzia = (term *(zcb/terma - (zab*dot)/(rab2*terma))*(dr+dr1))/ termc + (units.stbnunit*dt*ksb11*zab)/rab;

               dedxic = (term *(xab/terma - (xcb*dot)/(rcb2*terma))*(dr+dr1))/ termc + (units.stbnunit*dt*ksb21*xcb)/rcb;
               dedyic = (term *(yab/terma - (ycb*dot)/(rcb2*terma))*(dr+dr1))/ termc + (units.stbnunit*dt*ksb21*ycb)/rcb;
               dedzic = (term *(zab/terma - (zcb*dot)/(rcb2*terma))*(dr+dr1))/ termc + (units.stbnunit*dt*ksb21*zcb)/rcb;

               dedxib = ( (term /(rab*rcb))*( (xcb*dot)/rcb2 + (xab*dot)/rab2 + (2.0*xib-xia-xic))*(dr+dr1))/termc +
                             units.stbnunit*dt*( -(ksb11*xab)/rab - ksb21*xcb/rcb);
               dedyib = ( (term /(rab*rcb))*( (ycb*dot)/rcb2 + (yab*dot)/rab2 + (2.0*yib-yia-yic))*(dr+dr1))/termc +
                             units.stbnunit*dt*( -(ksb11*yab)/rab - ksb21*ycb/rcb);
               dedzib = ( (term /(rab*rcb))*( (zcb*dot)/rcb2 + (zab*dot)/rab2 + (2.0*zib-zia-zic))*(dr+dr1))/termc +
                             units.stbnunit*dt*( -(ksb11*zab)/rab - ksb21*zcb/rcb);
              
               *estrbnd += e;
               destb[ia][0] += dedxia;
               destb[ia][1] += dedyia;
               destb[ia][2] += dedzia;

               destb[ib][0] += dedxib;
               destb[ib][1] += dedyib;
               destb[ib][2] += dedzib;

               destb[ic][0] += dedxic;
               destb[ic][1] += dedyic;
               destb[ic][2] += dedzic;
               
            }
          }
      }
}
// =============================
void estrbnd2(int iatom,int nstrbnd,int **isb,int **i13,int *use,double *x,double *y,double *z,float *anat,double *bl,double *bk,float *ksb1,
	      float *ksb2,float **hessx,float **hessy,float **hessz)
{
    int i,j,k,ij;
    int ia,ib,ic;
    double dt,dr,angle,dot,cosine;
    double xia,yia,zia;
    double xib,yib,zib;
    double xic,yic,zic;
    double xab,yab,zab,rab;
    double xcb,ycb,zcb,rcb;
    double xp,yp,zp,rp,rp2;
    double terma,termc;
    double xrab,yrab,zrab,rab2;
    double xrcb,yrcb,zrcb,rcb2;
    double xabp,yabp,zabp;
    double xcbp,ycbp,zcbp;
    double ddtdxia,ddtdyia,ddtdzia;
    double ddtdxib,ddtdyib,ddtdzib;
    double ddtdxic,ddtdyic,ddtdzic;
    double ddrdxia,ddrdyia,ddrdzia;
    double ddrdxib,ddrdyib,ddrdzib;
    double ddrdxic,ddrdyic,ddrdzic;
    double dtxiaxia,dtxiayia,dtxiazia;
    double dtxibxib,dtxibyib,dtxibzib;
    double dtxicxic,dtxicyic,dtxiczic;
    double dtyiayia,dtyiazia,dtziazia;
    double dtyibyib,dtyibzib,dtzibzib;
    double dtyicyic,dtyiczic,dtziczic;
    double dtxibxia,dtxibyia,dtxibzia;
    double dtyibxia,dtyibyia,dtyibzia;
    double dtzibxia,dtzibyia,dtzibzia;
    double dtxibxic,dtxibyic,dtxibzic;
    double dtyibxic,dtyibyic,dtyibzic;
    double dtzibxic,dtzibyic,dtzibzic;
    double dtxiaxic,dtxiayic,dtxiazic;
    double dtyiaxic,dtyiayic,dtyiazic;
    double dtziaxic,dtziayic,dtziazic;
    double drxiaxia,drxiayia,drxiazia;
    double drxibxib,drxibyib,drxibzib;
    double drxicxic,drxicyic,drxiczic;
    double dryiayia,dryiazia,drziazia;
    double dryibyib,dryibzib,drzibzib;
    double dryicyic,dryiczic,drziczic;
    
      for (i=0; i < nstrbnd; i++)
      {
          ij = isb[i][0];
          ia = i13[ij][0];
          ib = i13[ij][1];
          ic = i13[ij][2];
          if (iatom == ia || iatom == ib || iatom == ic)
          {
            xia = x[ia];
            yia = y[ia];
            zia = z[ia];
            xib = x[ib];
            yib = y[ib];
            zib = z[ib];
            xic = x[ic];
            yic = y[ic];
            zic = z[ic];
            xab = xia - xib;
            yab = yia - yib;
            zab = zia - zib;
            rab = sqrt(xab*xab + yab*yab + zab*zab);
            xcb = xic - xib;
            ycb = yic - yib;
            zcb = zic - zib;
            rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb);
            xp = ycb*zab - zcb*yab;
            yp = zcb*xab - xcb*zab;
            zp = xcb*yab - ycb*xab;
            rp = sqrt(xp*xp + yp*yp + zp*zp);
            if (rp != 0.0)
            {
              dot = xab*xcb + yab*ycb + zab*zcb;
              cosine = dot / (rab*rcb);
               if (cosine < -1.0)
                 cosine = -1.0;
               if (cosine > 1.0)
                 cosine = 1.0;
               angle = radian * acos(cosine);

               dt = angle - anat[ij];
               terma = -radian / (rab*rab*rp);
               termc = radian / (rcb*rcb*rp);
               ddtdxia = terma * (yab*zp-zab*yp);
               ddtdyia = terma * (zab*xp-xab*zp);
               ddtdzia = terma * (xab*yp-yab*xp);
               ddtdxic = termc * (ycb*zp-zcb*yp);
               ddtdyic = termc * (zcb*xp-xcb*zp);
               ddtdzic = termc * (xcb*yp-ycb*xp);
               ddtdxib = -ddtdxia - ddtdxic;
               ddtdyib = -ddtdyia - ddtdyic;
               ddtdzib = -ddtdzia - ddtdzic;

               rab2 = 2.0 / (rab*rab);
               xrab = xab * rab2;
               yrab = yab * rab2;
               zrab = zab * rab2;
               rcb2 = 2.0 / (rcb*rcb);
               xrcb = xcb * rcb2;
               yrcb = ycb * rcb2;
               zrcb = zcb * rcb2;
               rp2 = 1.0 / (rp*rp);
               xabp = (yab*zp-zab*yp) * rp2;
               yabp = (zab*xp-xab*zp) * rp2;
               zabp = (xab*yp-yab*xp) * rp2;
               xcbp = (ycb*zp-zcb*yp) * rp2;
               ycbp = (zcb*xp-xcb*zp) * rp2;
               zcbp = (xcb*yp-ycb*xp) * rp2;

               dtxiaxia = terma*(xab*xcb-dot) + ddtdxia*(xcbp-xrab);
               dtxiayia = terma*(zp+yab*xcb) + ddtdxia*(ycbp-yrab);
               dtxiazia = terma*(zab*xcb-yp) + ddtdxia*(zcbp-zrab);
               dtyiayia = terma*(yab*ycb-dot) + ddtdyia*(ycbp-yrab);
               dtyiazia = terma*(xp+zab*ycb) + ddtdyia*(zcbp-zrab);
               dtziazia = terma*(zab*zcb-dot) + ddtdzia*(zcbp-zrab);
               dtxicxic = termc*(dot-xab*xcb) - ddtdxic*(xabp+xrcb);
               dtxicyic = termc*(zp-ycb*xab) - ddtdxic*(yabp+yrcb);
               dtxiczic = -termc*(yp+zcb*xab) - ddtdxic*(zabp+zrcb);
               dtyicyic = termc*(dot-yab*ycb) - ddtdyic*(yabp+yrcb);
               dtyiczic = termc*(xp-zcb*yab) - ddtdyic*(zabp+zrcb);
               dtziczic = termc*(dot-zab*zcb) - ddtdzic*(zabp+zrcb);
               dtxiaxic = terma*(yab*yab+zab*zab) - ddtdxia*xabp;
               dtxiayic = -terma*xab*yab - ddtdxia*yabp;
               dtxiazic = -terma*xab*zab - ddtdxia*zabp;
               dtyiaxic = -terma*xab*yab - ddtdyia*xabp;
               dtyiayic = terma*(xab*xab+zab*zab) - ddtdyia*yabp;
               dtyiazic = -terma*yab*zab - ddtdyia*zabp;
               dtziaxic = -terma*xab*zab - ddtdzia*xabp;
               dtziayic = -terma*yab*zab - ddtdzia*yabp;
               dtziazic = terma*(xab*xab+yab*yab) - ddtdzia*zabp;

               dtxibxia = -dtxiaxia - dtxiaxic;
               dtxibyia = -dtxiayia - dtyiaxic;
               dtxibzia = -dtxiazia - dtziaxic;
               dtyibxia = -dtxiayia - dtxiayic;
               dtyibyia = -dtyiayia - dtyiayic;
               dtyibzia = -dtyiazia - dtziayic;
               dtzibxia = -dtxiazia - dtxiazic;
               dtzibyia = -dtyiazia - dtyiazic;
               dtzibzia = -dtziazia - dtziazic;
               dtxibxic = -dtxicxic - dtxiaxic;
               dtxibyic = -dtxicyic - dtxiayic;
               dtxibzic = -dtxiczic - dtxiazic;
               dtyibxic = -dtxicyic - dtyiaxic;
               dtyibyic = -dtyicyic - dtyiayic;
               dtyibzic = -dtyiczic - dtyiazic;
               dtzibxic = -dtxiczic - dtziaxic;
               dtzibyic = -dtyiczic - dtziayic;
               dtzibzic = -dtziczic - dtziazic;
               dtxibxib = -dtxibxia - dtxibxic;
               dtxibyib = -dtxibyia - dtxibyic;
               dtxibzib = -dtxibzia - dtxibzic;
               dtyibyib = -dtyibyia - dtyibyic;
               dtyibzib = -dtyibzia - dtyibzic;
               dtzibzib = -dtzibzia - dtzibzic;

               j = isb[i][1];
               k = isb[i][2];
               dr = 0.0;
               terma = 0.0;
               termc = 0.0;
               if (j > -1)
               {
                  terma = units.stbnunit* ksb1[i] ;
                  dr = terma * (rab-bl[j]);
                  terma = terma / rab;
               }
               if (k > -1)
               {
                  termc = units.stbnunit* ksb2[i];
                  dr = dr + termc*(rcb-bl[k]);
                  termc = termc / rcb;
               }
               ddrdxia = terma * xab;
               ddrdyia = terma * yab;
               ddrdzia = terma * zab;
               ddrdxic = termc * xcb;
               ddrdyic = termc * ycb;
               ddrdzic = termc * zcb;
               ddrdxib = -ddrdxia - ddrdxic;
               ddrdyib = -ddrdyia - ddrdyic;
               ddrdzib = -ddrdzia - ddrdzic;

               xab = xab / rab;
               yab = yab / rab;
               zab = zab / rab;
               xcb = xcb / rcb;
               ycb = ycb / rcb;
               zcb = zcb / rcb;

               drxiaxia = terma * (1.0-xab*xab);
               drxiayia = -terma * xab*yab;
               drxiazia = -terma * xab*zab;
               dryiayia = terma * (1.0-yab*yab);
               dryiazia = -terma * yab*zab;
               drziazia = terma * (1.0-zab*zab);
               drxicxic = termc * (1.0-xcb*xcb);
               drxicyic = -termc * xcb*ycb;
               drxiczic = -termc * xcb*zcb;
               dryicyic = termc * (1.0-ycb*ycb);
               dryiczic = -termc * ycb*zcb;
               drziczic = termc * (1.0-zcb*zcb);
               drxibxib = drxiaxia + drxicxic;
               drxibyib = drxiayia + drxicyic;
               drxibzib = drxiazia + drxiczic;
               dryibyib = dryiayia + dryicyic;
               dryibzib = dryiazia + dryiczic;
               drzibzib = drziazia + drziczic;

            if (ia == iatom)
            {
               hessx[ia][0] += dt*drxiaxia + dr*dtxiaxia + 2.0*ddtdxia*ddrdxia;
               hessx[ia][1] += dt*drxiayia + dr*dtxiayia + ddtdxia*ddrdyia + ddtdyia*ddrdxia;
               hessx[ia][2] += dt*drxiazia + dr*dtxiazia + ddtdxia*ddrdzia + ddtdzia*ddrdxia;
               hessy[ia][0] += dt*drxiayia + dr*dtxiayia + ddtdyia*ddrdxia + ddtdxia*ddrdyia;
               hessy[ia][1] += dt*dryiayia + dr*dtyiayia + 2.0*ddtdyia*ddrdyia;
               hessy[ia][2] += dt*dryiazia + dr*dtyiazia + ddtdyia*ddrdzia + ddtdzia*ddrdyia;
               hessz[ia][0] += dt*drxiazia + dr*dtxiazia + ddtdzia*ddrdxia + ddtdxia*ddrdzia;
               hessz[ia][1] += dt*dryiazia + dr*dtyiazia + ddtdzia*ddrdyia + ddtdyia*ddrdzia;
               hessz[ia][2] += dt*drziazia + dr*dtziazia + 2.0*ddtdzia*ddrdzia;
               hessx[ib][0] += -dt*drxiaxia + dr*dtxibxia + ddtdxia*ddrdxib + ddtdxib*ddrdxia;
               hessx[ib][1] += -dt*drxiayia + dr*dtxibyia + ddtdxia*ddrdyib + ddtdyib*ddrdxia;
               hessx[ib][2] += -dt*drxiazia + dr*dtxibzia + ddtdxia*ddrdzib + ddtdzib*ddrdxia;
               hessy[ib][0] += -dt*drxiayia + dr*dtyibxia + ddtdyia*ddrdxib + ddtdxib*ddrdyia;
               hessy[ib][1] += -dt*dryiayia + dr*dtyibyia + ddtdyia*ddrdyib + ddtdyib*ddrdyia;
               hessy[ib][2] += -dt*dryiazia + dr*dtyibzia + ddtdyia*ddrdzib + ddtdzib*ddrdyia;
               hessz[ib][0] += -dt*drxiazia + dr*dtzibxia + ddtdzia*ddrdxib + ddtdxib*ddrdzia;
               hessz[ib][1] += -dt*dryiazia + dr*dtzibyia + ddtdzia*ddrdyib + ddtdyib*ddrdzia;
               hessz[ib][2] += -dt*drziazia + dr*dtzibzia + ddtdzia*ddrdzib + ddtdzib*ddrdzia;
               hessx[ic][0] += dr*dtxiaxic + ddtdxia*ddrdxic + ddtdxic*ddrdxia;
               hessx[ic][1] += dr*dtxiayic + ddtdxia*ddrdyic + ddtdyic*ddrdxia;
               hessx[ic][2] += dr*dtxiazic + ddtdxia*ddrdzic + ddtdzic*ddrdxia;
               hessy[ic][0] += dr*dtyiaxic + ddtdyia*ddrdxic + ddtdxic*ddrdyia;
               hessy[ic][1] += dr*dtyiayic + ddtdyia*ddrdyic + ddtdyic*ddrdyia;
               hessy[ic][2] += dr*dtyiazic + ddtdyia*ddrdzic + ddtdzic*ddrdyia;
               hessz[ic][0] += dr*dtziaxic + ddtdzia*ddrdxic + ddtdxic*ddrdzia;
               hessz[ic][1] += dr*dtziayic + ddtdzia*ddrdyic + ddtdyic*ddrdzia;
               hessz[ic][2] += dr*dtziazic + ddtdzia*ddrdzic + ddtdzic*ddrdzia;
            } else if (ib == iatom)
            {
               hessx[ib][0] += dt*drxibxib + dr*dtxibxib + 2.0*ddtdxib*ddrdxib;
               hessx[ib][1] += dt*drxibyib + dr*dtxibyib + ddtdxib*ddrdyib + ddtdyib*ddrdxib;
               hessx[ib][2] += dt*drxibzib + dr*dtxibzib + ddtdxib*ddrdzib + ddtdzib*ddrdxib;
               hessy[ib][0] += dt*drxibyib + dr*dtxibyib + ddtdyib*ddrdxib + ddtdxib*ddrdyib;
               hessy[ib][1] += dt*dryibyib + dr*dtyibyib + 2.0*ddtdyib*ddrdyib;
               hessy[ib][2] += dt*dryibzib + dr*dtyibzib + ddtdyib*ddrdzib + ddtdzib*ddrdyib;
               hessz[ib][0] += dt*drxibzib + dr*dtxibzib + ddtdzib*ddrdxib + ddtdxib*ddrdzib;
               hessz[ib][1] += dt*dryibzib + dr*dtyibzib + ddtdzib*ddrdyib + ddtdyib*ddrdzib;
               hessz[ib][2] += dt*drzibzib + dr*dtzibzib + 2.0*ddtdzib*ddrdzib;
               hessx[ia][0] += -dt*drxiaxia + dr*dtxibxia + ddtdxib*ddrdxia + ddtdxia*ddrdxib;
               hessx[ia][1] += -dt*drxiayia + dr*dtxibyia + ddtdxib*ddrdyia + ddtdyia*ddrdxib;
               hessx[ia][2] += -dt*drxiazia + dr*dtxibzia + ddtdxib*ddrdzia + ddtdzia*ddrdxib;
               hessy[ia][0] += -dt*drxiayia + dr*dtyibxia + ddtdyib*ddrdxia + ddtdxia*ddrdyib;
               hessy[ia][1] += -dt*dryiayia + dr*dtyibyia + ddtdyib*ddrdyia + ddtdyia*ddrdyib;
               hessy[ia][2] += -dt*dryiazia + dr*dtyibzia + ddtdyib*ddrdzia + ddtdzia*ddrdyib;
               hessz[ia][0] += -dt*drxiazia + dr*dtzibxia + ddtdzib*ddrdxia + ddtdxia*ddrdzib;
               hessz[ia][1] += -dt*dryiazia + dr*dtzibyia + ddtdzib*ddrdyia + ddtdyia*ddrdzib;
               hessz[ia][2] += -dt*drziazia + dr*dtzibzia + ddtdzib*ddrdzia + ddtdzia*ddrdzib;
               hessx[ic][0] += -dt*drxicxic + dr*dtxibxic + ddtdxib*ddrdxic + ddtdxic*ddrdxib;
               hessx[ic][1] += -dt*drxicyic + dr*dtxibyic + ddtdxib*ddrdyic + ddtdyic*ddrdxib;
               hessx[ic][2] += -dt*drxiczic + dr*dtxibzic + ddtdxib*ddrdzic + ddtdzic*ddrdxib;
               hessy[ic][0] += -dt*drxicyic + dr*dtyibxic + ddtdyib*ddrdxic + ddtdxic*ddrdyib;
               hessy[ic][1] += -dt*dryicyic + dr*dtyibyic + ddtdyib*ddrdyic + ddtdyic*ddrdyib;
               hessy[ic][2] += -dt*dryiczic + dr*dtyibzic + ddtdyib*ddrdzic + ddtdzic*ddrdyib;
               hessz[ic][0] += -dt*drxiczic + dr*dtzibxic + ddtdzib*ddrdxic + ddtdxic*ddrdzib;
               hessz[ic][1] += -dt*dryiczic + dr*dtzibyic + ddtdzib*ddrdyic + ddtdyic*ddrdzib;
               hessz[ic][2] += -dt*drziczic + dr*dtzibzic + ddtdzib*ddrdzic + ddtdzic*ddrdzib;
            }else if (ic == iatom)
            {
               hessx[ic][0] += dt*drxicxic + dr*dtxicxic + 2.0*ddtdxic*ddrdxic;
               hessx[ic][1] += dt*drxicyic + dr*dtxicyic + ddtdxic*ddrdyic + ddtdyic*ddrdxic;
               hessx[ic][2] += dt*drxiczic + dr*dtxiczic + ddtdxic*ddrdzic + ddtdzic*ddrdxic;
               hessy[ic][0] += dt*drxicyic + dr*dtxicyic + ddtdyic*ddrdxic + ddtdxic*ddrdyic;
               hessy[ic][1] += dt*dryicyic + dr*dtyicyic + 2.0*ddtdyic*ddrdyic;
               hessy[ic][2] += dt*dryiczic + dr*dtyiczic + ddtdyic*ddrdzic + ddtdzic*ddrdyic;
               hessz[ic][0] += dt*drxiczic + dr*dtxiczic + ddtdzic*ddrdxic + ddtdxic*ddrdzic;
               hessz[ic][1] += dt*dryiczic + dr*dtyiczic + ddtdzic*ddrdyic + ddtdyic*ddrdzic;
               hessz[ic][2] += dt*drziczic + dr*dtziczic + 2.0*ddtdzic*ddrdzic;
               hessx[ib][0] += -dt*drxicxic + dr*dtxibxic + ddtdxic*ddrdxib + ddtdxib*ddrdxic;
               hessx[ib][1] += -dt*drxicyic + dr*dtxibyic + ddtdxic*ddrdyib + ddtdyib*ddrdxic;
               hessx[ib][2] += -dt*drxiczic + dr*dtxibzic + ddtdxic*ddrdzib + ddtdzib*ddrdxic;
               hessy[ib][0] += -dt*drxicyic + dr*dtyibxic + ddtdyic*ddrdxib + ddtdxib*ddrdyic;
               hessy[ib][1] += -dt*dryicyic + dr*dtyibyic + ddtdyic*ddrdyib + ddtdyib*ddrdyic;
               hessy[ib][2] += -dt*dryiczic + dr*dtyibzic + ddtdyic*ddrdzib + ddtdzib*ddrdyic;
               hessz[ib][0] += -dt*drxiczic + dr*dtzibxic + ddtdzic*ddrdxib + ddtdxib*ddrdzic;
               hessz[ib][1] += -dt*dryiczic + dr*dtzibyic + ddtdzic*ddrdyib + ddtdyib*ddrdzic;
               hessz[ib][2] += -dt*drziczic + dr*dtzibzic + ddtdzic*ddrdzib + ddtdzib*ddrdzic;
               hessx[ia][0] += dr*dtxiaxic + ddtdxic*ddrdxia + ddtdxia*ddrdxic;
               hessx[ia][1] += dr*dtyiaxic + ddtdxic*ddrdyia + ddtdyia*ddrdxic;
               hessx[ia][2] += dr*dtziaxic + ddtdxic*ddrdzia + ddtdzia*ddrdxic;
               hessy[ia][0] += dr*dtxiayic + ddtdyic*ddrdxia + ddtdxia*ddrdyic;
               hessy[ia][1] += dr*dtyiayic + ddtdyic*ddrdyia + ddtdyia*ddrdyic;
               hessy[ia][2] += dr*dtziayic + ddtdyic*ddrdzia + ddtdzia*ddrdyic;
               hessz[ia][0] += dr*dtxiazic + ddtdzic*ddrdxia + ddtdxia*ddrdzic;
               hessz[ia][1] += dr*dtyiazic + ddtdzic*ddrdyia + ddtdyia*ddrdzic;
               hessz[ia][2] += dr*dtziazic + ddtdzic*ddrdzia + ddtdzia*ddrdzic;
            }
            }
          }
      }
}
                  
