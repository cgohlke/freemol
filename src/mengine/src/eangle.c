#define EXTERN extern

#include "pcwin.h"

//double acos(double);


EXTERN struct t_minim_values {
        int iprint, ndc, nconst;
        float dielc;
        } minim_values;

void eangle(int nang,int **,int *use,double *x,double *y,double *z,int *angin,int *angtype,float *anat,float *acon,double *ebend);
void eangle1(int nang,int natom,int **,int *use,double *x,double *y,double *z,int *angin,int *angtype,float *anat,float *acon,double *ebend,double **dea);
void eangle2(int i,int nang,int **,int *angin, double *x,double *y,double *z,float *acon,float *anat,double **dea,float **hessx,float **hessy,float **hessz);
void eangle2a(int iatom,int nang,int **,int *angin,double *x,double *y,double *z,float *acon,float *anat,float **hessx,float **hessy,float **hessz);
void eangle2b(int i,int nang,int **,int *angin,double *x,double *y,double *z,float *acon,float *anat,double **dea);
// ======================        
void eangle(int nang,int **i13,int *use,double *x,double *y,double *z,int *angin,int *angtype,float *anat,float *acon,double *ebend)
{
    int i, ia, ib, ic, id;
    double e, angle, dot, cosine;
    double dt, dt2, dt3, dt4;
    double xia, yia, zia;
    double xib, yib, zib;
    double xic, yic, zic;
    double xab, yab, zab, rab2;
    double xcb, ycb, zcb, rcb2;


    *ebend = 0.0;
    if (minim_values.iprint)
    {
        fprintf(pcmlogfile,"\nAngle Terms \n");
        fprintf(pcmlogfile,"             At1      At2      At3     Angle     Thet0   Tconst      Ebend\n");
    }

    for (i=0; i < nang; i++)
    {
        ia = i13[i][0];
        ib = i13[i][1];
        ic = i13[i][2];
        id = i13[i][3];
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
                xcb = xic - xib;
                ycb = yic - yib;
                zcb = zic - zib;
                rcb2 = xcb*xcb + ycb*ycb + zcb*zcb;
                if (rab2 != 0.0 && rcb2 != 0.00)
                {
                    dot = xab*xcb + yab*ycb + zab*zcb;
                    cosine = dot / sqrt(rab2*rcb2);
                    if (cosine > 1.0)
                       cosine = 1.0;
                    if (cosine < -1.0)
                       cosine = -1.0;
                    angle = radian*acos(cosine);

                        dt = angle - anat[i];
                        dt2 = dt*dt;
                        dt3 = dt2*dt;
                        dt4 = dt2*dt2;
                        e = units.angunit*acon[i]*dt2
                          *(1.0 + units.cang*dt + units.qang*dt2 + units.pang*dt3 + units.sang*dt4);
                    
                    *ebend += e;
                    if (minim_values.iprint)
                    {
                          fprintf(pcmlogfile,"Angle 1-3: (%-3d)- (%-3d)- (%-3d)  %-8.3f  %-8.3f %-8.4f = %-8.4f\n",
                             ia, ib, ic, angle, anat[i],acon[i], e);
                    }
                }   
        }         
    }
}
// ===============================================
void eangle1(int nang,int natom,int **i13,int *use,double *x,double *y,double *z,int *angin,int *angtype,float *anat,float *acon,double *ebend,double **dea)
{
  int i, ia, ib, ic, id;
    double e, angle, dot, cosine;
    double dt, dt2, dt3, dt4;
    double deddt,terma,termc;
    double xia, yia, zia;
    double xib, yib, zib;
    double xic, yic, zic;
    double xab, yab, zab, rab2;
    double xcb, ycb, zcb, rcb2;
    double xp,yp,zp,rp;
    double dedxia,dedyia,dedzia;
    double dedxib,dedyib,dedzib;
    double dedxic,dedyic,dedzic;
    double dtemp;

    *ebend = 0.0F;
      for (i=0; i <= natom; i++)
      {
          dea[i][0] = 0.0;
          dea[i][1] = 0.0;
          dea[i][2] = 0.0;
      }
    for (i=0; i < nang; i++)
    {
        ia = i13[i][0];
        ib = i13[i][1];
        ic = i13[i][2];
        id = i13[i][3];
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
                xcb = xic - xib;
                ycb = yic - yib;
                zcb = zic - zib;
                rcb2 = xcb*xcb + ycb*ycb + zcb*zcb;
                xp = ycb*zab - zcb*yab;
                yp = zcb*xab - xcb*zab;
                zp = xcb*yab - ycb*xab;
                rp = sqrt(xp*xp + yp*yp + zp*zp);
                if (rp != 0.0 )
                {
                    dot = xab*xcb + yab*ycb + zab*zcb;
                    cosine = dot / sqrt(rab2*rcb2);
                    if (cosine > 1.0)
                       cosine = 1.0;
                    if (cosine < -1.0)
                       cosine = -1.0;
                    angle = radian*acos(cosine);
                        dt = angle - anat[i];
                        dt2 = dt*dt;
                        dt3 = dt2*dt;
                        dt4 = dt2*dt2;

                        e = units.angunit*acon[i]*dt2
                          *(1.0 + units.cang*dt + units.qang*dt2 + units.pang*dt3 + units.sang*dt4);
                        deddt = units.angunit*acon[i]* dt * radian
                            * (2.0 + 3.0*units.cang*dt + 4.0*units.qang*dt2
                                 + 5.0*units.pang*dt3 + 6.0*units.sang*dt4);
                   

                   dtemp = dot*dot/(rab2*rcb2);
                   if (fabs(dtemp) >= 1.0)
                     dtemp = 0.9999*dtemp/(fabs(dtemp));
                   terma = sqrt(1.0-dtemp);
                   termc = sqrt(rab2*rcb2);
                   dedxia = -deddt*( xcb/termc - (xab*cosine)/rab2)/terma;
                   dedyia = -deddt*( ycb/termc - (yab*cosine)/rab2)/terma;
                   dedzia = -deddt*( zcb/termc - (zab*cosine)/rab2)/terma;

                   dedxic = -deddt*( xab/termc - (xcb*cosine)/rcb2)/terma;
                   dedyic = -deddt*( yab/termc - (ycb*cosine)/rcb2)/terma;
                   dedzic = -deddt*( zab/termc - (zcb*cosine)/rcb2)/terma;

                   dedxib = -deddt*( (-xab-xcb)/termc + (dot*(xcb*rab2 + xab*rcb2))/(rab2*rcb2*termc))/terma;
                   dedyib = -deddt*( (-yab-ycb)/termc + (dot*(ycb*rab2 + yab*rcb2))/(rab2*rcb2*termc))/terma;
                   dedzib = -deddt*( (-zab-zcb)/termc + (dot*(zcb*rab2 + zab*rcb2))/(rab2*rcb2*termc))/terma;  
                  
                              
                   dea[ia][0] += dedxia;
                   dea[ia][1] += dedyia;
                   dea[ia][2] += dedzia;

                   dea[ib][0] += dedxib;
                   dea[ib][1] += dedyib;
                   dea[ib][2] += dedzib;

                   dea[ic][0] += dedxic;
                   dea[ic][1] += dedyic;
                   dea[ic][2] += dedzic;
                   
                   *ebend += e;
                }  
        }         
    }
}
// =========================
void eangle2(int i,int nang,int **i13,int *angin, double *x,double *y,double *z,float *acon,float *anat,double **dea,float **hessx,float **hessy,float **hessz)
{
   eangle2a(i,nang,i13,angin,x,y,z,acon,anat,hessx,hessy,hessz);
}
// =================================================                    
void eangle2a(int iatom,int nang,int **i13,int *angin,double *x,double *y,double *z,float *acon,float *anat,float **hessx,float **hessy,float **hessz)
{
    int i,ia,ib,ic;
    double angle,dot,cosine;
    double dt,dt2,dt3,dt4;
    double xia,yia,zia;
    double xib,yib,zib;
    double xic,yic,zic;
    double xab,yab,zab,rab;
    double xcb,ycb,zcb,rcb;
    double xp,yp,zp,rp,rp2;
    double xrab,yrab,zrab,rab2;
    double xrcb,yrcb,zrcb,rcb2;
    double xabp,yabp,zabp;
    double xcbp,ycbp,zcbp;
    double deddt,d2eddt2,terma,termc;
    double ddtdxia,ddtdyia,ddtdzia;
    double ddtdxib,ddtdyib,ddtdzib;
    double ddtdxic,ddtdyic,ddtdzic;
    double dxiaxia,dxiayia,dxiazia;
    double dxibxib,dxibyib,dxibzib;
    double dxicxic,dxicyic,dxiczic;
    double dyiayia,dyiazia,dziazia;
    double dyibyib,dyibzib,dzibzib;
    double dyicyic,dyiczic,dziczic;
    double dxibxia,dxibyia,dxibzia;
    double dyibxia,dyibyia,dyibzia;
    double dzibxia,dzibyia,dzibzia;
    double dxibxic,dxibyic,dxibzic;
    double dyibxic,dyibyic,dyibzic;
    double dzibxic,dzibyic,dzibzic;
    double dxiaxic,dxiayic,dxiazic;
    double dyiaxic,dyiayic,dyiazic;
    double dziaxic,dziayic,dziazic;

    for (i=0; i < nang; i++)
    {
        ia = i13[i][0];
        ib = i13[i][1];
        ic = i13[i][2];
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

                      dt = angle - anat[i];
                      dt2 = dt * dt;
                      dt3 = dt2 * dt;
                      dt4 = dt3 * dt;
                      deddt = units.angunit * acon[i] * dt 
                           * (2.0 + 3.0*units.cang*dt + 4.0*units.qang*dt2
                               + 5.0*units.pang*dt3 + 6.0*units.sang*dt4);
                      d2eddt2 = units.angunit * acon[i]
                            * (2.0 + 6.0*units.cang*dt + 12.0*units.qang*dt2
                                + 20.0*units.pang*dt3 + 30.0*units.sang*dt4);
                        terma = -radian/ (rab*rab*rp);
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

                  dxiaxia = terma*(xab*xcb-dot) + ddtdxia*(xcbp-xrab);
                  dxiayia = terma*(zp+yab*xcb) + ddtdxia*(ycbp-yrab);
                  dxiazia = terma*(zab*xcb-yp) + ddtdxia*(zcbp-zrab);
                  dyiayia = terma*(yab*ycb-dot) + ddtdyia*(ycbp-yrab);
                  dyiazia = terma*(xp+zab*ycb) + ddtdyia*(zcbp-zrab);
                  dziazia = terma*(zab*zcb-dot) + ddtdzia*(zcbp-zrab);
                  dxicxic = termc*(dot-xab*xcb) - ddtdxic*(xabp+xrcb);
                  dxicyic = termc*(zp-ycb*xab) - ddtdxic*(yabp+yrcb);
                  dxiczic = -termc*(yp+zcb*xab) - ddtdxic*(zabp+zrcb);
                  dyicyic = termc*(dot-yab*ycb) - ddtdyic*(yabp+yrcb);
                  dyiczic = termc*(xp-zcb*yab) - ddtdyic*(zabp+zrcb);
                  dziczic = termc*(dot-zab*zcb) - ddtdzic*(zabp+zrcb);
                  dxiaxic = terma*(yab*yab+zab*zab) - ddtdxia*xabp;
                  dxiayic = -terma*xab*yab - ddtdxia*yabp;
                  dxiazic = -terma*xab*zab - ddtdxia*zabp;
                  dyiaxic = -terma*xab*yab - ddtdyia*xabp;
                  dyiayic = terma*(xab*xab+zab*zab) - ddtdyia*yabp;
                  dyiazic = -terma*yab*zab - ddtdyia*zabp;
                  dziaxic = -terma*xab*zab - ddtdzia*xabp;
                  dziayic = -terma*yab*zab - ddtdzia*yabp;
                  dziazic = terma*(xab*xab+yab*yab) - ddtdzia*zabp;

                  dxibxia = -dxiaxia - dxiaxic;
                  dxibyia = -dxiayia - dyiaxic;
                  dxibzia = -dxiazia - dziaxic;
                  dyibxia = -dxiayia - dxiayic;
                  dyibyia = -dyiayia - dyiayic;
                  dyibzia = -dyiazia - dziayic;
                  dzibxia = -dxiazia - dxiazic;
                  dzibyia = -dyiazia - dyiazic;
                  dzibzia = -dziazia - dziazic;
                  dxibxic = -dxicxic - dxiaxic;
                  dxibyic = -dxicyic - dxiayic;
                  dxibzic = -dxiczic - dxiazic;
                  dyibxic = -dxicyic - dyiaxic;
                  dyibyic = -dyicyic - dyiayic;
                  dyibzic = -dyiczic - dyiazic;
                  dzibxic = -dxiczic - dziaxic;
                  dzibyic = -dyiczic - dziayic;
                  dzibzic = -dziczic - dziazic;
                  dxibxib = -dxibxia - dxibxic;
                  dxibyib = -dxibyia - dxibyic;
                  dxibzib = -dxibzia - dxibzic;
                  dyibyib = -dyibyia - dyibyic;
                  dyibzib = -dyibzia - dyibzic;
                  dzibzib = -dzibzia - dzibzic;

                  if (ia == iatom)
                  {
                     hessx[ia][0] += deddt*dxiaxia + d2eddt2*ddtdxia*ddtdxia;
                     hessx[ia][1] += deddt*dxiayia + d2eddt2*ddtdxia*ddtdyia;
                     hessx[ia][2] += deddt*dxiazia + d2eddt2*ddtdxia*ddtdzia;
                     hessy[ia][0] += deddt*dxiayia + d2eddt2*ddtdyia*ddtdxia;
                     hessy[ia][1] += deddt*dyiayia + d2eddt2*ddtdyia*ddtdyia;
                     hessy[ia][2] += deddt*dyiazia + d2eddt2*ddtdyia*ddtdzia;
                     hessz[ia][0] += deddt*dxiazia + d2eddt2*ddtdzia*ddtdxia;
                     hessz[ia][1] += deddt*dyiazia + d2eddt2*ddtdzia*ddtdyia;
                     hessz[ia][2] += deddt*dziazia + d2eddt2*ddtdzia*ddtdzia;
                     hessx[ib][0] += deddt*dxibxia + d2eddt2*ddtdxia*ddtdxib;
                     hessx[ib][1] += deddt*dyibxia + d2eddt2*ddtdxia*ddtdyib;
                     hessx[ib][2] += deddt*dzibxia + d2eddt2*ddtdxia*ddtdzib;
                     hessy[ib][0] += deddt*dxibyia + d2eddt2*ddtdyia*ddtdxib;
                     hessy[ib][1] += deddt*dyibyia + d2eddt2*ddtdyia*ddtdyib;
                     hessy[ib][2] += deddt*dzibyia + d2eddt2*ddtdyia*ddtdzib;
                     hessz[ib][0] += deddt*dxibzia + d2eddt2*ddtdzia*ddtdxib;
                     hessz[ib][1] += deddt*dyibzia + d2eddt2*ddtdzia*ddtdyib;
                     hessz[ib][2] += deddt*dzibzia + d2eddt2*ddtdzia*ddtdzib;
                     hessx[ic][0] += deddt*dxiaxic + d2eddt2*ddtdxia*ddtdxic;
                     hessx[ic][1] += deddt*dxiayic + d2eddt2*ddtdxia*ddtdyic;
                     hessx[ic][2] += deddt*dxiazic + d2eddt2*ddtdxia*ddtdzic;
                     hessy[ic][0] += deddt*dyiaxic + d2eddt2*ddtdyia*ddtdxic;
                     hessy[ic][1] += deddt*dyiayic + d2eddt2*ddtdyia*ddtdyic;
                     hessy[ic][2] += deddt*dyiazic + d2eddt2*ddtdyia*ddtdzic;
                     hessz[ic][0] += deddt*dziaxic + d2eddt2*ddtdzia*ddtdxic;
                     hessz[ic][1] += deddt*dziayic + d2eddt2*ddtdzia*ddtdyic;
                     hessz[ic][2] += deddt*dziazic + d2eddt2*ddtdzia*ddtdzic;
                  } else if (ib == iatom)
                  {
                     hessx[ib][0] += deddt*dxibxib + d2eddt2*ddtdxib*ddtdxib;
                     hessx[ib][1] += deddt*dxibyib + d2eddt2*ddtdxib*ddtdyib;
                     hessx[ib][2] += deddt*dxibzib + d2eddt2*ddtdxib*ddtdzib;
                     hessy[ib][0] += deddt*dxibyib + d2eddt2*ddtdyib*ddtdxib;
                     hessy[ib][1] += deddt*dyibyib + d2eddt2*ddtdyib*ddtdyib;
                     hessy[ib][2] += deddt*dyibzib + d2eddt2*ddtdyib*ddtdzib;
                     hessz[ib][0] += deddt*dxibzib + d2eddt2*ddtdzib*ddtdxib;
                     hessz[ib][1] += deddt*dyibzib + d2eddt2*ddtdzib*ddtdyib;
                     hessz[ib][2] += deddt*dzibzib + d2eddt2*ddtdzib*ddtdzib;
                     hessx[ia][0] += deddt*dxibxia + d2eddt2*ddtdxib*ddtdxia;
                     hessx[ia][1] += deddt*dxibyia + d2eddt2*ddtdxib*ddtdyia;
                     hessx[ia][2] += deddt*dxibzia + d2eddt2*ddtdxib*ddtdzia;
                     hessy[ia][0] += deddt*dyibxia + d2eddt2*ddtdyib*ddtdxia;
                     hessy[ia][1] += deddt*dyibyia + d2eddt2*ddtdyib*ddtdyia;
                     hessy[ia][2] += deddt*dyibzia + d2eddt2*ddtdyib*ddtdzia;
                     hessz[ia][0] += deddt*dzibxia + d2eddt2*ddtdzib*ddtdxia;
                     hessz[ia][1] += deddt*dzibyia + d2eddt2*ddtdzib*ddtdyia;
                     hessz[ia][2] += deddt*dzibzia + d2eddt2*ddtdzib*ddtdzia;
                     hessx[ic][0] += deddt*dxibxic + d2eddt2*ddtdxib*ddtdxic;
                     hessx[ic][1] += deddt*dxibyic + d2eddt2*ddtdxib*ddtdyic;
                     hessx[ic][2] += deddt*dxibzic + d2eddt2*ddtdxib*ddtdzic;
                     hessy[ic][0] += deddt*dyibxic + d2eddt2*ddtdyib*ddtdxic;
                     hessy[ic][1] += deddt*dyibyic + d2eddt2*ddtdyib*ddtdyic;
                     hessy[ic][2] += deddt*dyibzic + d2eddt2*ddtdyib*ddtdzic;
                     hessz[ic][0] += deddt*dzibxic + d2eddt2*ddtdzib*ddtdxic;
                     hessz[ic][1] += deddt*dzibyic + d2eddt2*ddtdzib*ddtdyic;
                     hessz[ic][2] += deddt*dzibzic + d2eddt2*ddtdzib*ddtdzic;
                  }else if (ic == iatom)
                  {
                     hessx[ic][0] += deddt*dxicxic + d2eddt2*ddtdxic*ddtdxic;
                     hessx[ic][1] += deddt*dxicyic + d2eddt2*ddtdxic*ddtdyic;
                     hessx[ic][2] += deddt*dxiczic + d2eddt2*ddtdxic*ddtdzic;
                     hessy[ic][0] += deddt*dxicyic + d2eddt2*ddtdyic*ddtdxic;
                     hessy[ic][1] += deddt*dyicyic + d2eddt2*ddtdyic*ddtdyic;
                     hessy[ic][2] += deddt*dyiczic + d2eddt2*ddtdyic*ddtdzic;
                     hessz[ic][0] += deddt*dxiczic + d2eddt2*ddtdzic*ddtdxic;
                     hessz[ic][1] += deddt*dyiczic + d2eddt2*ddtdzic*ddtdyic;
                     hessz[ic][2] += deddt*dziczic + d2eddt2*ddtdzic*ddtdzic;
                     hessx[ib][0] += deddt*dxibxic + d2eddt2*ddtdxic*ddtdxib;
                     hessx[ib][1] += deddt*dyibxic + d2eddt2*ddtdxic*ddtdyib;
                     hessx[ib][2] += deddt*dzibxic + d2eddt2*ddtdxic*ddtdzib;
                     hessy[ib][0] += deddt*dxibyic + d2eddt2*ddtdyic*ddtdxib;
                     hessy[ib][1] += deddt*dyibyic + d2eddt2*ddtdyic*ddtdyib;
                     hessy[ib][2] += deddt*dzibyic + d2eddt2*ddtdyic*ddtdzib;
                     hessz[ib][0] += deddt*dxibzic + d2eddt2*ddtdzic*ddtdxib;
                     hessz[ib][1] += deddt*dyibzic + d2eddt2*ddtdzic*ddtdyib;
                     hessz[ib][2] += deddt*dzibzic + d2eddt2*ddtdzic*ddtdzib;
                     hessx[ia][0] += deddt*dxiaxic + d2eddt2*ddtdxic*ddtdxia;
                     hessx[ia][1] += deddt*dyiaxic + d2eddt2*ddtdxic*ddtdyia;
                     hessx[ia][2] += deddt*dziaxic + d2eddt2*ddtdxic*ddtdzia;
                     hessy[ia][0] += deddt*dxiayic + d2eddt2*ddtdyic*ddtdxia;
                     hessy[ia][1] += deddt*dyiayic + d2eddt2*ddtdyic*ddtdyia;
                     hessy[ia][2] += deddt*dziayic + d2eddt2*ddtdyic*ddtdzia;
                     hessz[ia][0] += deddt*dxiazic + d2eddt2*ddtdzic*ddtdxia;
                     hessz[ia][1] += deddt*dyiazic + d2eddt2*ddtdzic*ddtdyia;
                     hessz[ia][2] += deddt*dziazic + d2eddt2*ddtdzic*ddtdzia;
                  }
               }
            }
    }
}
// ================================================================
void eangle2b(int i,int nang,int **i13,int *angin,double *x,double *y,double *z,float *acon,float *anat,double **dea)
{
    int ia,ib,ic,id;
    double angle,dot,cosine;
    double dt,dt2,dt3,dt4;
    double deddt,terma,termc,term;
    double xia,yia,zia;
    double xib,yib,zib;
    double xic,yic,zic;
    double xid,yid,zid;
    double xad,yad,zad;
    double xbd,ybd,zbd;
    double xcd,ycd,zcd;
    double xip,yip,zip;
    double xap,yap,zap,rap2;
    double xcp,ycp,zcp,rcp2;
    double xt,yt,zt,rt2,ptrt2;
    double xm,ym,zm,rm,delta,delta2;
    double dedxia,dedyia,dedzia;
    double dedxib,dedyib,dedzib;
    double dedxic,dedyic,dedzic;
    double dedxid,dedyid,dedzid;
    double dedxip,dedyip,dedzip;
    double dpdxia,dpdyia,dpdzia;
    double dpdxic,dpdyic,dpdzic;

    ia = i13[i][0];
    ib = i13[i][1];
    ic = i13[i][2];
    id = i13[i][3];
    xia = x[ia];
    yia = y[ia];
    zia = z[ia];
    xib = x[ib];
    yib = y[ib];
    zib = z[ib];
    xic = x[ic];
    yic = y[ic];
    zic = z[ic];
    xid = x[id];
    yid = y[id];
    zid = z[id];

    dea[ia][0] = 0.0;
    dea[ia][1] = 0.0;
    dea[ia][2] = 0.0;
    dea[ib][0] = 0.0;
    dea[ib][1] = 0.0;
    dea[ib][2] = 0.0;
    dea[ic][0] = 0.0;
    dea[ic][1] = 0.0;
    dea[ic][2] = 0.0;
    dea[id][0] = 0.0;
    dea[id][1] = 0.0;
    dea[id][2] = 0.0;

      xad = xia - xid;
      yad = yia - yid;
      zad = zia - zid;
      xbd = xib - xid;
      ybd = yib - yid;
      zbd = zib - zid;
      xcd = xic - xid;
      ycd = yic - yid;
      zcd = zic - zid;
      xt = yad*zcd - zad*ycd;
      yt = zad*xcd - xad*zcd;
      zt = xad*ycd - yad*xcd;
      rt2 = xt*xt + yt*yt + zt*zt;
      delta = -(xt*xbd + yt*ybd + zt*zbd) / rt2;
      xip = xib + xt*delta;
      yip = yib + yt*delta;
      zip = zib + zt*delta;
      xap = xia - xip;
      yap = yia - yip;
      zap = zia - zip;
      rap2 = xap*xap + yap*yap + zap*zap;
      xcp = xic - xip;
      ycp = yic - yip;
      zcp = zic - zip;
      rcp2 = xcp*xcp + ycp*ycp + zcp*zcp;
      xm = ycp*zap - zcp*yap;
      ym = zcp*xap - xcp*zap;
      zm = xcp*yap - ycp*xap;
      rm = sqrt(xm*xm + ym*ym + zm*zm);
      if (rm != 0.0)
      {
         dot = xap*xcp + yap*ycp + zap*zcp;
         cosine = dot / sqrt(rap2*rcp2);
         if (cosine < -1.0)
           cosine = -1.0;
         if (cosine > 1.0)
           cosine = 1.0;
         angle = radian * acos(cosine);

         dt = angle - anat[i];
         dt2 = dt * dt;
         dt3 = dt2 * dt;
         dt4 = dt2 * dt2;
         deddt = units.angunit * acon[i] * dt
                 * (2.0 + 3.0*units.cang*dt + 4.0*units.qang*dt2
                    + 5.0*units.pang*dt3 + 6.0*units.sang*dt4);

         deddt = deddt * radian;
         terma = -deddt / (rap2*rm);
         termc = deddt / (rcp2*rm);
         dedxia = terma * (yap*zm-zap*ym);
         dedyia = terma * (zap*xm-xap*zm);
         dedzia = terma * (xap*ym-yap*xm);
         dedxic = termc * (ycp*zm-zcp*ym);
         dedyic = termc * (zcp*xm-xcp*zm);
         dedzic = termc * (xcp*ym-ycp*xm);
         dedxip = -dedxia - dedxic;
         dedyip = -dedyia - dedyic;
         dedzip = -dedzia - dedzic;

         delta2 = 2.0 * delta;
         ptrt2 = (dedxip*xt + dedyip*yt + dedzip*zt) / rt2;
         term = (zcd*ybd-ycd*zbd) + delta2*(yt*zcd-zt*ycd);
         dpdxia = delta*(ycd*dedzip-zcd*dedyip) + term*ptrt2;
         term = (xcd*zbd-zcd*xbd) + delta2*(zt*xcd-xt*zcd);
         dpdyia = delta*(zcd*dedxip-xcd*dedzip) + term*ptrt2;
         term = (ycd*xbd-xcd*ybd) + delta2*(xt*ycd-yt*xcd);
         dpdzia = delta*(xcd*dedyip-ycd*dedxip) + term*ptrt2;
         term = (yad*zbd-zad*ybd) + delta2*(zt*yad-yt*zad);
         dpdxic = delta*(zad*dedyip-yad*dedzip) + term*ptrt2;
         term = (zad*xbd-xad*zbd) + delta2*(xt*zad-zt*xad);
         dpdyic = delta*(xad*dedzip-zad*dedxip) + term*ptrt2;
         term = (xad*ybd-yad*xbd) + delta2*(yt*xad-xt*yad);
         dpdzic = delta*(yad*dedxip-xad*dedyip) + term*ptrt2;
         dedxia = dedxia + dpdxia;
         dedyia = dedyia + dpdyia;
         dedzia = dedzia + dpdzia;
         dedxib = dedxip;
         dedyib = dedyip;
         dedzib = dedzip;
         dedxic = dedxic + dpdxic;
         dedyic = dedyic + dpdyic;
         dedzic = dedzic + dpdzic;
         dedxid = -dedxia - dedxib - dedxic;
         dedyid = -dedyia - dedyib - dedyic;
         dedzid = -dedzia - dedzib - dedzic;

         dea[ia][0] = dedxia;
         dea[ia][1] = dedyia;
         dea[ia][2] = dedzia;
         dea[ib][0] = dedxib;
         dea[ib][1] = dedyib;
         dea[ib][2] = dedzib;
         dea[ic][0] = dedxic;
         dea[ic][1] = dedyic;
         dea[ic][2] = dedzic;
         dea[id][0] = dedxid;
         dea[id][1] = dedyid;
         dea[id][2] = dedzid;
      }
}

