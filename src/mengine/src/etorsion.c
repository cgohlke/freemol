#define EXTERN extern

#include "pcwin.h"

double dihdrl(int,int,int,int);

EXTERN struct t_minim_values {
        int iprint, ndc, nconst;
        float dielc;
        } minim_values;

void etorsion(int ntor,int **i14,int *use,int *type,double *x,double *y,double *z,int *phin1,int *phin2,int *phin3,int *phin4,int *phin5,int *phin6,float *vin1,float *vin2,
              float *vin3,float *vin4,float *vin5,float *vin6,double *etor);
void etorsion1(int natom,int ntor,int **i14,int *use,double *x,double *y,double *z,int *phin1,int *phin2,int *phin3,int *phin4,int *phin5,int *phin6,float *vin1,float *vin2,
             float *vin3,float *vin4,float *vin5,float *vin6,double *etor,double **detor);
void etorsion2(int iatom,int ntor,int **i14,int *use,double *x,double *y,double *z,int *phin1,int *phin2,int *phin3,int *phin4,int *phin5,int *phin6,float *vin1,float *vin2,        
             float *vin3,float *vin4,float *vin5,float *vin6,float **hessx,float **hessy,float **hessz);
// ==============================
void etorsion(int ntor,int **i14,int *use,int *type,double *x,double *y,double *z,int *phin1,int *phin2,int *phin3,int *phin4,int *phin5,int *phin6,float *vin1,float *vin2,float *vin3,float *vin4,
	      float *vin5,float *vin6,double *etor)
{
    int i,ia, ib, ic, id;
    double e, rt2, ru2, rtru;
    double xt,yt,zt,xu,yu,zu;
    double v1,v2,v3,v4,v5,v6;
    double s1,s2,s3,s4,s5,s6;
    double cosine,cosine2,cosine3;
    double cosine4,cosine5,cosine6;
    double phi1,phi2,phi3;
    double phi4,phi5,phi6;
    double xia,yia,zia,xib,yib,zib;
    double xic,yic,zic,xid,yid,zid;
    double xba,yba,zba,xcb,ycb,zcb;
    double xdc,ydc,zdc;
    double rdihed;
  
    
    *etor = 0.0;
    if (minim_values.iprint)
    {
      fprintf(pcmlogfile,"\nTorsion Terms\n");
      fprintf(pcmlogfile,"        At1      At2      At3      At4       Types          Angle     V1      V2        V3         Etor\n");
    }

    
    for (i=0; i < ntor; i++)
    {
        ia = i14[i][0];
        ib = i14[i][1];
        ic = i14[i][2];
        id = i14[i][3];
          if (use[ia] || use[ib] || use[ic] || use[id] )
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
            xid = x[id];
            yid = y[id];
            zid = z[id];
            xba = xib - xia;
            yba = yib - yia;
            zba = zib - zia;
            xcb = xic - xib;
            ycb = yic - yib;
            zcb = zic - zib;
            xdc = xid - xic;
            ydc = yid - yic;
            zdc = zid - zic;
            xt = yba*zcb - ycb*zba;
            yt = zba*xcb - zcb*xba;
            zt = xba*ycb - xcb*yba;
            xu = ycb*zdc - ydc*zcb;
            yu = zcb*xdc - zdc*xcb;
            zu = xcb*ydc - xdc*ycb;
            rt2 = xt*xt + yt*yt + zt*zt;
            ru2 = xu*xu + yu*yu + zu*zu;
            rtru = sqrt(rt2 * ru2);
            if (rtru != 0.0)
            { 
               cosine = (xt*xu + yt*yu + zt*zu) / rtru;

               v1 = vin1[i];
               s1 = phin1[i];
               v2 = vin2[i];
               s2 = phin2[i];
               v3 = vin3[i];
               s3 = phin3[i];
               v4 = vin4[i];
               s4 = phin4[i];
               v5 = vin5[i];
               s5 = phin5[i];
               v6 = vin6[i];
               s6 = phin6[i];

//     compute the powers of the cosine of this angle

               cosine2 = cosine * cosine;
               cosine3 = cosine2 * cosine;
               cosine4 = cosine2 * cosine2;
               cosine5 = cosine3 * cosine2;
               cosine6 = cosine3 * cosine3;

               phi1 = 1.0 + s1*cosine;
               phi2 = 1.0 + s2*(2.0*cosine2 - 1.0);
               phi3 = 1.0 + s3*(4.0*cosine3 - 3.0*cosine);
               phi4 = 1.0 + s4*(8.0*(cosine4-cosine2) + 1.0);
               phi5 = 1.0 + s5*(16.0*cosine5 - 20.0*cosine3 + 5.0*cosine);
               phi6 = 1.0 + s6*(32.0*cosine6 - 48.0*cosine4 + 18.0*cosine2 - 1.0);

               e = units.torsunit * (v1*phi1 + v2*phi2 + v3*phi3
                                   + v4*phi4 + v5*phi5 + v6*phi6);
               *etor += e;
               if(minim_values.iprint)
               {
                   rdihed = dihdrl(ia,ib,ic,id);
                   fprintf(pcmlogfile,"Tor:  (%-3d)- (%-3d)- (%-3d)- (%-3d)    %d %d %d %d     %8.2f    %-8.3f %-8.3f %-8.3f = %-8.4f\n",
                   ia,ib,ic,id, type[ia], type[ib],
                    type[ic], type[id], rdihed,v1,v2,v3,e);
               }
            }
        }
     }
}
// ===========================================================
void etorsion1(int natom,int ntor,int **i14,int *use,double *x,double *y,double *z,int *phin1,int *phin2,int *phin3,int *phin4,int *phin5,int *phin6,float *vin1,float *vin2,float *vin3,float *vin4,
	       float *vin5,float *vin6,double *etor,double **detor)
{
    int i,ia, ib, ic, id;
    double e, rt2, ru2, rtru,dedphi,sine;
    double xt,yt,zt,xu,yu,zu;
    double v1,v2,v3,v4,v5,v6;
    double s1,s2,s3,s4,s5,s6;
    double cosine,cosine2,cosine3;
    double cosine4,cosine5,cosine6;
    double phi1,phi2,phi3;
    double phi4,phi5,phi6;
    double xia,yia,zia,xib,yib,zib;
    double xic,yic,zic,xid,yid,zid;
    double xba,yba,zba,xcb,ycb,zcb;
    double xdc,ydc,zdc;
    double xca,yca,zca,xdb,ydb,zdb;
    double xtu,ytu,ztu,rcb;
    double dphi1,dphi2,dphi3;
    double dphi4,dphi5,dphi6;
    double dphidxt,dphidyt,dphidzt;
    double dphidxu,dphidyu,dphidzu;
    double dphidxia,dphidyia,dphidzia;
    double dphidxib,dphidyib,dphidzib;
    double dphidxic,dphidyic,dphidzic;
    double dphidxid,dphidyid,dphidzid;
    double dedxia,dedyia,dedzia;
    double dedxib,dedyib,dedzib;
    double dedxic,dedyic,dedzic;
    double dedxid,dedyid,dedzid;
    
    *etor = 0.0;
      for (i=0; i <= natom; i++)
      {
          detor[i][0] = 0.0;
          detor[i][1] = 0.0;
          detor[i][2] = 0.0;
      }
    
    for (i=0; i < ntor; i++)
    {
        ia = i14[i][0];
        ib = i14[i][1];
        ic = i14[i][2];
        id = i14[i][3];
          if (use[ia] || use[ib] || use[ic] || use[id] )
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
            xid = x[id];
            yid = y[id];
            zid = z[id];
            xba = xib - xia;
            yba = yib - yia;
            zba = zib - zia;
            xcb = xic - xib;
            ycb = yic - yib;
            zcb = zic - zib;
            xdc = xid - xic;
            ydc = yid - yic;
            zdc = zid - zic;
            xt = yba*zcb - ycb*zba;
            yt = zba*xcb - zcb*xba;
            zt = xba*ycb - xcb*yba;
            xu = ycb*zdc - ydc*zcb;
            yu = zcb*xdc - zdc*xcb;
            zu = xcb*ydc - xdc*ycb;
            xtu = yt*zu - yu*zt;
            ytu = zt*xu - zu*xt;
            ztu = xt*yu - xu*yt;
            rt2 = xt*xt + yt*yt + zt*zt;
            ru2 = xu*xu + yu*yu + zu*zu;
            rtru = sqrt(rt2 * ru2);
            if (rtru != 0.0)
            {
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb);
               cosine = (xt*xu + yt*yu + zt*zu) / rtru;  
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru);

               v1 = vin1[i];
               s1 = phin1[i];
               v2 = vin2[i];
               s2 = phin2[i];
               v3 = vin3[i];
               s3 = phin3[i];
               v4 = vin4[i];
               s4 = phin4[i];
               v5 = vin5[i];
               s5 = phin5[i];
               v6 = vin6[i];
               s6 = phin6[i];
//     compute the powers of the cosine of this angle

               cosine2 = cosine * cosine;
               cosine3 = cosine2 * cosine;
               cosine4 = cosine2 * cosine2;
               cosine5 = cosine3 * cosine2;
               cosine6 = cosine3 * cosine3;

               phi1 = 1.0 + s1*cosine;
               phi2 = 1.0 + s2*(2.0*cosine2 - 1.0);
               phi3 = 1.0 + s3*(4.0*cosine3 - 3.0*cosine);
               phi4 = 1.0 + s4*(8.0*(cosine4-cosine2) + 1.0);
               phi5 = 1.0 + s5*(16.0*cosine5 - 20.0*cosine3 + 5.0*cosine);
               phi6 = 1.0 + s6*(32.0*cosine6 - 48.0*cosine4 + 18.0*cosine2 - 1.0);

               dphi1 = s1;
               dphi2 = s2 * (4.0*cosine);
               dphi3 = s3 * (12.0*cosine2 - 3.0);
               dphi4 = s4 * (32.0*cosine3 - 16.0*cosine);
               dphi5 = s5 * (80.0*cosine4 - 60.0*cosine2 + 5.0);
               dphi6 = s6 * (192.0*(cosine5-cosine3) + 36.0*cosine);

               e = units.torsunit * (v1*phi1 + v2*phi2 + v3*phi3
                                   + v4*phi4 + v5*phi5 + v6*phi6);
               dedphi = -sine * units.torsunit
                       * (v1*dphi1+v2*dphi2+v3*dphi3+v4*dphi4+v5*dphi5+v6*dphi6);

               xca = xic - xia;
               yca = yic - yia;
               zca = zic - zia;
               xdb = xid - xib;
               ydb = yid - yib;
               zdb = zid - zib;
               dphidxt = (yt*zcb - ycb*zt) / (rt2*rcb);
               dphidyt = (zt*xcb - zcb*xt) / (rt2*rcb);
               dphidzt = (xt*ycb - xcb*yt) / (rt2*rcb);
               dphidxu = -(yu*zcb - ycb*zu) / (ru2*rcb);
               dphidyu = -(zu*xcb - zcb*xu) / (ru2*rcb);
               dphidzu = -(xu*ycb - xcb*yu) / (ru2*rcb);

               dphidxia = zcb*dphidyt - ycb*dphidzt;
               dphidyia = xcb*dphidzt - zcb*dphidxt;
               dphidzia = ycb*dphidxt - xcb*dphidyt;
               dphidxib = yca*dphidzt - zca*dphidyt + zdc*dphidyu - ydc*dphidzu;
               dphidyib = zca*dphidxt - xca*dphidzt + xdc*dphidzu - zdc*dphidxu;
               dphidzib = xca*dphidyt - yca*dphidxt + ydc*dphidxu - xdc*dphidyu;
               dphidxic = zba*dphidyt - yba*dphidzt + ydb*dphidzu - zdb*dphidyu;
               dphidyic = xba*dphidzt - zba*dphidxt + zdb*dphidxu - xdb*dphidzu;
               dphidzic = yba*dphidxt - xba*dphidyt + xdb*dphidyu - ydb*dphidxu;
               dphidxid = zcb*dphidyu - ycb*dphidzu;
               dphidyid = xcb*dphidzu - zcb*dphidxu;
               dphidzid = ycb*dphidxu - xcb*dphidyu;
               dedxia = dedphi * dphidxia;
               dedyia = dedphi * dphidyia;
               dedzia = dedphi * dphidzia;
               dedxib = dedphi * dphidxib;
               dedyib = dedphi * dphidyib;
               dedzib = dedphi * dphidzib;
               dedxic = dedphi * dphidxic;
               dedyic = dedphi * dphidyic;
               dedzic = dedphi * dphidzic;
               dedxid = dedphi * dphidxid;
               dedyid = dedphi * dphidyid;
               dedzid = dedphi * dphidzid;

               *etor += e;
               detor[ia][0] += dedxia;
               detor[ia][1] += dedyia;
               detor[ia][2] += dedzia;
               
               detor[ib][0] += dedxib;
               detor[ib][1] += dedyib;
               detor[ib][2] += dedzib;

               detor[ic][0] += dedxic;
               detor[ic][1] += dedyic;
               detor[ic][2] += dedzic;

               detor[id][0] += dedxid;
               detor[id][1] += dedyid;
               detor[id][2] += dedzid;
               
            }
        }
     }
}
// ===========================================================
void etorsion2(int iatom,int ntor,int **i14,int *use,double *x,double *y,double *z,int *phin1,int *phin2,int *phin3,int *phin4,int *phin5,int *phin6,float *vin1,float *vin2,
           float *vin3,float *vin4,float *vin5,float *vin6,float **hessx,float **hessy,float **hessz)
{
    int i,ia,ib,ic,id;
    double dedphi,d2edphi2,sine;
    double v1,v2,v3,v4,v5,v6;
    double s1,s2,s3,s4,s5,s6;
    double cosine,cosine2,cosine3;
    double cosine4,cosine5,cosine6;
    double xia,yia,zia,xib,yib,zib;
    double xic,yic,zic,xid,yid,zid;
    double xba,yba,zba,xcb,ycb,zcb;
    double xdc,ydc,zdc;
    double xca,yca,zca,xdb,ydb,zdb;
    double xt,yt,zt,xu,yu,zu,xtu,ytu,ztu;
    double rt2,ru2,rtru,rcb;
    double phi1,phi2,phi3;
    double phi4,phi5,phi6;
    double dphi1,dphi2,dphi3;
    double dphi4,dphi5,dphi6;
    double d2phi1,d2phi2,d2phi3;
    double d2phi4,d2phi5,d2phi6;
    double dphidxt,dphidyt,dphidzt;
    double dphidxu,dphidyu,dphidzu;
    double dphidxia,dphidyia,dphidzia;
    double dphidxib,dphidyib,dphidzib;
    double dphidxic,dphidyic,dphidzic;
    double dphidxid,dphidyid,dphidzid;
    double xycb2,xzcb2,yzcb2;
    double rcbxt,rcbyt,rcbzt,rcbt2;
    double rcbxu,rcbyu,rcbzu,rcbu2;
    double dphidxibt,dphidyibt,dphidzibt;
    double dphidxibu,dphidyibu,dphidzibu;
    double dphidxict,dphidyict,dphidzict;
    double dphidxicu,dphidyicu,dphidzicu;
    double dxiaxia,dyiayia,dziazia,dxibxib,dyibyib,dzibzib;
    double dxicxic,dyicyic,dziczic,dxidxid,dyidyid,dzidzid;
    double dxiayia,dxiazia,dyiazia,dxibyib,dxibzib,dyibzib;
    double dxicyic,dxiczic,dyiczic,dxidyid,dxidzid,dyidzid;
    double dxiaxib,dxiayib,dxiazib,dyiaxib,dyiayib,dyiazib;
    double dziaxib,dziayib,dziazib,dxiaxic,dxiayic,dxiazic;
    double dyiaxic,dyiayic,dyiazic,dziaxic,dziayic,dziazic;
    double dxiaxid,dxiayid,dxiazid,dyiaxid,dyiayid,dyiazid;
    double dziaxid,dziayid,dziazid,dxibxic,dxibyic,dxibzic;
    double dyibxic,dyibyic,dyibzic,dzibxic,dzibyic,dzibzic;
    double dxibxid,dxibyid,dxibzid,dyibxid,dyibyid,dyibzid;
    double dzibxid,dzibyid,dzibzid,dxicxid,dxicyid,dxiczid;
    double dyicxid,dyicyid,dyiczid,dzicxid,dzicyid,dziczid;

    for (i=0; i < ntor; i++)
    {
        ia = i14[i][0];
        ib = i14[i][1];
        ic = i14[i][2];
        id = i14[i][3];
        if (ia == iatom || ib == iatom || ic == iatom || id == iatom)
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
            xid = x[id];
            yid = y[id];
            zid = z[id];
            xba = xib - xia;
            yba = yib - yia;
            zba = zib - zia;
            xcb = xic - xib;
            ycb = yic - yib;
            zcb = zic - zib;
            xdc = xid - xic;
            ydc = yid - yic;
            zdc = zid - zic;
            xt = yba*zcb - ycb*zba;
            yt = zba*xcb - zcb*xba;
            zt = xba*ycb - xcb*yba;
            xu = ycb*zdc - ydc*zcb;
            yu = zcb*xdc - zdc*xcb;
            zu = xcb*ydc - xdc*ycb;
            xtu = yt*zu - yu*zt;
            ytu = zt*xu - zu*xt;
            ztu = xt*yu - xu*yt;
            rt2 = xt*xt + yt*yt + zt*zt;
            ru2 = xu*xu + yu*yu + zu*zu;
            rtru = sqrt(rt2 * ru2);
            if (rtru != 0.0)
            {
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb);
               cosine = (xt*xu + yt*yu + zt*zu) / rtru;
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru);

               v1 = vin1[i];
               s1 = phin1[i];
               v2 = vin2[i];
               s2 = phin2[i];
               v3 = vin3[i];
               s3 = phin3[i];
               v4 = vin4[i];
               s4 = phin4[i];
               v5 = vin5[i];
               s5 = phin5[i];
               v6 = vin6[i];
               s6 = phin6[i];

//     compute the powers of the cosine of this angle

               cosine2 = cosine * cosine;
               cosine3 = cosine2 * cosine;
               cosine4 = cosine2 * cosine2;
               cosine5 = cosine3 * cosine2;
               cosine6 = cosine3 * cosine3;

               phi1 = 1.0 + s1*cosine;
               phi2 = 1.0 + s2*(2.0*cosine2 - 1.0);
               phi3 = 1.0 + s3*(4.0*cosine3 - 3.0*cosine);
               phi4 = 1.0 + s4*(8.0*(cosine4-cosine2) + 1.0);
               phi5 = 1.0 + s5*(16.0*cosine5 - 20.0*cosine3 + 5.0*cosine);
               phi6 = 1.0 + s6*(32.0*cosine6 - 48.0*cosine4 + 18.0*cosine2 - 1.0);

               dphi1 = s1;
               dphi2 = s2 * (4.0*cosine);
               dphi3 = s3 * (12.0*cosine2 - 3.0);
               dphi4 = s4 * (32.0*cosine3 - 16.0*cosine);
               dphi5 = s5 * (80.0*cosine4 - 60.0*cosine2 + 5.0);
               dphi6 = s6 * (192.0*(cosine5-cosine3) + 36.0*cosine);
            
               d2phi1 = -s1 * cosine;
               d2phi2 = -s2 * (8.0*cosine2 - 4.0);
               d2phi3 = -s3 * (36.0*cosine3 - 27.0*cosine);
               d2phi4 = -s4 * (128.0*(cosine4-cosine2) + 16.0);
               d2phi5 = -s5 * (400.0*cosine5 - 500.0*cosine + 125.0*cosine);
               d2phi6 = -s6 * (1152.0*cosine6 - 1728.0*cosine2 + 648.0*cosine2);

               dedphi = -sine * units.torsunit * (v1*dphi1 + v2*dphi2 + v3*dphi3
                             + v4*dphi4 + v5*dphi5 + v6*dphi6);
               d2edphi2 = units.torsunit * (v1*d2phi1 + v2*d2phi2 + v3*d2phi3
                             + v4*d2phi4 + v5*d2phi5 + v6*d2phi6);

               xca = xic - xia;
               yca = yic - yia;
               zca = zic - zia;
               xdb = xid - xib;
               ydb = yid - yib;
               zdb = zid - zib;
               dphidxt = (yt*zcb - ycb*zt) / (rt2*rcb);
               dphidyt = (zt*xcb - zcb*xt) / (rt2*rcb);
               dphidzt = (xt*ycb - xcb*yt) / (rt2*rcb);
               dphidxu = -(yu*zcb - ycb*zu) / (ru2*rcb);
               dphidyu = -(zu*xcb - zcb*xu) / (ru2*rcb);
               dphidzu = -(xu*ycb - xcb*yu) / (ru2*rcb);

               xycb2 = xcb*xcb + ycb*ycb;
               xzcb2 = xcb*xcb + zcb*zcb;
               yzcb2 = ycb*ycb + zcb*zcb;
               rcbxt = -2.0 * rcb * dphidxt;
               rcbyt = -2.0 * rcb * dphidyt;
               rcbzt = -2.0 * rcb * dphidzt;
               rcbt2 = rcb * rt2;
               rcbxu = 2.0 * rcb * dphidxu;
               rcbyu = 2.0 * rcb * dphidyu;
               rcbzu = 2.0 * rcb * dphidzu;
               rcbu2 = rcb * ru2;
               dphidxibt = yca*dphidzt - zca*dphidyt;
               dphidxibu = zdc*dphidyu - ydc*dphidzu;
               dphidyibt = zca*dphidxt - xca*dphidzt;
               dphidyibu = xdc*dphidzu - zdc*dphidxu;
               dphidzibt = xca*dphidyt - yca*dphidxt;
               dphidzibu = ydc*dphidxu - xdc*dphidyu;
               dphidxict = zba*dphidyt - yba*dphidzt;
               dphidxicu = ydb*dphidzu - zdb*dphidyu;
               dphidyict = xba*dphidzt - zba*dphidxt;
               dphidyicu = zdb*dphidxu - xdb*dphidzu;
               dphidzict = yba*dphidxt - xba*dphidyt;
               dphidzicu = xdb*dphidyu - ydb*dphidxu;
               
               dphidxia = zcb*dphidyt - ycb*dphidzt;
               dphidyia = xcb*dphidzt - zcb*dphidxt;
               dphidzia = ycb*dphidxt - xcb*dphidyt;
               dphidxib = dphidxibt + dphidxibu;
               dphidyib = dphidyibt + dphidyibu;
               dphidzib = dphidzibt + dphidzibu;
               dphidxic = dphidxict + dphidxicu;
               dphidyic = dphidyict + dphidyicu;
               dphidzic = dphidzict + dphidzicu;
               dphidxid = zcb*dphidyu - ycb*dphidzu;
               dphidyid = xcb*dphidzu - zcb*dphidxu;
               dphidzid = ycb*dphidxu - xcb*dphidyu;

               dxiaxia = rcbxt*dphidxia;
               dxiayia = rcbxt*dphidyia - zcb*rcb/rt2;
               dxiazia = rcbxt*dphidzia + ycb*rcb/rt2;
               dxiaxib = rcbxt*dphidxibt + xcb*(zca*ycb-yca*zcb)/rcbt2;
               dxiayib = rcbxt*dphidyibt + dphidzt + (xca*zcb*xcb+zca*yzcb2)/rcbt2;
               dxiazib = rcbxt*dphidzibt - dphidyt - (xca*ycb*xcb+yca*yzcb2)/rcbt2;
               dxiaxic = rcbxt*dphidxict + xcb*xt/rcbt2;
               dxiayic = rcbxt*dphidyict - dphidzt - (xba*zcb*xcb+zba*yzcb2)/rcbt2;
               dxiazic = rcbxt*dphidzict + dphidyt + (xba*ycb*xcb+yba*yzcb2)/rcbt2;
               dxiaxid = 0.0;
               dxiayid = 0.0;
               dxiazid = 0.0;
               dyiayia = rcbyt*dphidyia;
               dyiazia = rcbyt*dphidzia - xcb*rcb/rt2;
               dyiaxib = rcbyt*dphidxibt - dphidzt - (yca*zcb*ycb+zca*xzcb2)/rcbt2;
               dyiayib = rcbyt*dphidyibt + ycb*(xca*zcb-zca*xcb)/rcbt2;
               dyiazib = rcbyt*dphidzibt + dphidxt + (yca*xcb*ycb+xca*xzcb2)/rcbt2;
               dyiaxic = rcbyt*dphidxict + dphidzt + (yba*zcb*ycb+zba*xzcb2)/rcbt2;
               dyiayic = rcbyt*dphidyict + ycb*yt/rcbt2;
               dyiazic = rcbyt*dphidzict - dphidxt - (yba*xcb*ycb+xba*xzcb2)/rcbt2;
               dyiaxid = 0.0;
               dyiayid = 0.0;
               dyiazid = 0.0;
               dziazia = rcbzt*dphidzia;
               dziaxib = rcbzt*dphidxibt + dphidyt + (zca*ycb*zcb+yca*xycb2)/rcbt2;
               dziayib = rcbzt*dphidyibt - dphidxt - (zca*xcb*zcb+xca*xycb2)/rcbt2;
               dziazib = rcbzt*dphidzibt + zcb*(yca*xcb-xca*ycb)/rcbt2;
               dziaxic = rcbzt*dphidxict - dphidyt - (zba*ycb*zcb+yba*xycb2)/rcbt2;
               dziayic = rcbzt*dphidyict + dphidxt + (zba*xcb*zcb+xba*xycb2)/rcbt2;
               dziazic = rcbzt*dphidzict + zcb*zt/rcbt2;
               dziaxid = 0.0;
               dziayid = 0.0;
               dziazid = 0.0;
               dxibxic = -xcb*dphidxib/(rcb*rcb) - (yca*(zba*xcb+yt)-zca*(yba*xcb-zt))/rcbt2
                  - 2.0*(yt*zba-yba*zt)*dphidxibt/rt2 - (zdc*(ydb*xcb+zu)-ydc*(zdb*xcb-yu))/rcbu2
                  + 2.0*(yu*zdb-ydb*zu)*dphidxibu/ru2;
               dxibyic = -ycb*dphidxib/(rcb*rcb) + dphidzt + dphidzu - (yca*(zba*ycb-xt)+zca*(xba*xcb+zcb*zba))/rcbt2
                  - 2.0*(zt*xba-zba*xt)*dphidxibt/rt2 + (zdc*(xdb*xcb+zcb*zdb)+ydc*(zdb*ycb+xu))/rcbu2
                  + 2.0*(zu*xdb-zdb*xu)*dphidxibu/ru2;
               dxibxid = rcbxu*dphidxibu + xcb*xu/rcbu2;
               dxibyid = rcbyu*dphidxibu - dphidzu - (ydc*zcb*ycb+zdc*xzcb2)/rcbu2;
               dxibzid = rcbzu*dphidxibu + dphidyu + (zdc*ycb*zcb+ydc*xycb2)/rcbu2;
               dyibzib = ycb*dphidzib/(rcb*rcb) - (xca*(xca*xcb+zcb*zca)+yca*(ycb*xca+zt))/rcbt2
                  - 2.0*(xt*zca-xca*zt)*dphidzibt/rt2 + (ydc*(xdc*ycb-zu)+xdc*(xdc*xcb+zcb*zdc))/rcbu2
                  + 2.0*(xu*zdc-xdc*zu)*dphidzibu/ru2;
               dyibxic = -xcb*dphidyib/(rcb*rcb) - dphidzt - dphidzu + (xca*(zba*xcb+yt)+zca*(zba*zcb+ycb*yba))/rcbt2
                  - 2.0*(yt*zba-yba*zt)*dphidyibt/rt2 - (zdc*(zdb*zcb+ycb*ydb)+xdc*(zdb*xcb-yu))/rcbu2
                 + 2.0*(yu*zdb-ydb*zu)*dphidyibu/ru2;
               dyibyic = -ycb*dphidyib/(rcb*rcb)- (zca*(xba*ycb+zt)-xca*(zba*ycb-xt))/rcbt2
                  - 2.0*(zt*xba-zba*xt)*dphidyibt/rt2 - (xdc*(zdb*ycb+xu)-zdc*(xdb*ycb-zu))/rcbu2
                  + 2.0*(zu*xdb-zdb*xu)*dphidyibu/ru2;
               dyibxid = rcbxu*dphidyibu + dphidzu + (xdc*zcb*xcb+zdc*yzcb2)/rcbu2;
               dyibyid = rcbyu*dphidyibu + ycb*yu/rcbu2;
               dyibzid = rcbzu*dphidyibu - dphidxu - (zdc*xcb*zcb+xdc*xycb2)/rcbu2;
               dzibxic = -xcb*dphidzib/(rcb*rcb) + dphidyt + dphidyu - (xca*(yba*xcb-zt)+yca*(zba*zcb+ycb*yba))/rcbt2
                  - 2.0*(yt*zba-yba*zt)*dphidzibt/rt2 + (ydc*(zdb*zcb+ycb*ydb)+xdc*(ydb*xcb+zu))/rcbu2
                  + 2.0*(yu*zdb-ydb*zu)*dphidzibu/ru2;
               dzibzic = -zcb*dphidzib/(rcb*rcb) - (xca*(yba*zcb+xt)-yca*(xba*zcb-yt))/rcbt2
                  - 2.0*(xt*yba-xba*yt)*dphidzibt/rt2 - (ydc*(xdb*zcb+yu)-xdc*(ydb*zcb-xu))/rcbu2
                  + 2.0*(xu*ydb-xdb*yu)*dphidzibu/ru2;
               dzibxid = rcbxu*dphidzibu - dphidyu - (xdc*ycb*xcb+ydc*yzcb2)/rcbu2;
               dzibyid = rcbyu*dphidzibu + dphidxu + (ydc*xcb*ycb+xdc*xzcb2)/rcbu2;
               dzibzid = rcbzu*dphidzibu + zcb*zu/rcbu2;
               dxicxid = rcbxu*dphidxicu - xcb*(zdb*ycb-ydb*zcb)/rcbu2;
               dxicyid = rcbyu*dphidxicu + dphidzu + (ydb*zcb*ycb+zdb*xzcb2)/rcbu2;
               dxiczid = rcbzu*dphidxicu - dphidyu - (zdb*ycb*zcb+ydb*xycb2)/rcbu2;
               dyicxid = rcbxu*dphidyicu - dphidzu - (xdb*zcb*xcb+zdb*yzcb2)/rcbu2;
               dyicyid = rcbyu*dphidyicu - ycb*(xdb*zcb-zdb*xcb)/rcbu2;
               dyiczid = rcbzu*dphidyicu + dphidxu + (zdb*xcb*zcb+xdb*xycb2)/rcbu2;
               dzicxid = rcbxu*dphidzicu + dphidyu + (xdb*ycb*xcb+ydb*yzcb2)/rcbu2;
               dzicyid = rcbyu*dphidzicu - dphidxu - (ydb*xcb*ycb+xdb*xzcb2)/rcbu2;
               dziczid = rcbzu*dphidzicu - zcb*(ydb*xcb-xdb*ycb)/rcbu2;
               dxidxid = rcbxu*dphidxid;
               dxidyid = rcbxu*dphidyid + zcb*rcb/ru2;
               dxidzid = rcbxu*dphidzid - ycb*rcb/ru2;
               dyidyid = rcbyu*dphidyid;
               dyidzid = rcbyu*dphidzid + xcb*rcb/ru2;
               dzidzid = rcbzu*dphidzid;

               dxibxib = -dxiaxib - dxibxic - dxibxid;
               dxibyib = -dyiaxib - dxibyic - dxibyid;
               dxibzib = -dxiazib - dzibxic - dzibxid;
               dxibzic = -dziaxib - dxibzib - dxibzid;
               dyibyib = -dyiayib - dyibyic - dyibyid;
               dyibzic = -dziayib - dyibzib - dyibzid;
               dzibzib = -dziazib - dzibzic - dzibzid;
               dzibyic = -dyiazib - dyibzib - dzibyid;
               dxicxic = -dxiaxic - dxibxic - dxicxid;
               dxicyic = -dyiaxic - dyibxic - dxicyid;
               dxiczic = -dziaxic - dzibxic - dxiczid;
               dyicyic = -dyiayic - dyibyic - dyicyid;
               dyiczic = -dziayic - dzibyic - dyiczid;
               dziczic = -dziazic - dzibzic - dziczid;

         if (iatom == ia)
         {
            hessx[ia][0] += dedphi*dxiaxia + d2edphi2*dphidxia*dphidxia;
            hessy[ia][0] += dedphi*dxiayia + d2edphi2*dphidxia*dphidyia;
            hessz[ia][0] += dedphi*dxiazia + d2edphi2*dphidxia*dphidzia;
            hessx[ia][1] += dedphi*dxiayia + d2edphi2*dphidxia*dphidyia;
            hessy[ia][1] += dedphi*dyiayia + d2edphi2*dphidyia*dphidyia;
            hessz[ia][1] += dedphi*dyiazia + d2edphi2*dphidyia*dphidzia;
            hessx[ia][2] += dedphi*dxiazia + d2edphi2*dphidxia*dphidzia;
            hessy[ia][2] += dedphi*dyiazia + d2edphi2*dphidyia*dphidzia;
            hessz[ia][2] += dedphi*dziazia + d2edphi2*dphidzia*dphidzia;
            hessx[ib][0] += dedphi*dxiaxib + d2edphi2*dphidxia*dphidxib;
            hessy[ib][0] += dedphi*dyiaxib + d2edphi2*dphidyia*dphidxib;
            hessz[ib][0] += dedphi*dziaxib + d2edphi2*dphidzia*dphidxib;
            hessx[ib][1] += dedphi*dxiayib + d2edphi2*dphidxia*dphidyib;
            hessy[ib][1] += dedphi*dyiayib + d2edphi2*dphidyia*dphidyib;
            hessz[ib][1] += dedphi*dziayib + d2edphi2*dphidzia*dphidyib;
            hessx[ib][2] += dedphi*dxiazib + d2edphi2*dphidxia*dphidzib;
            hessy[ib][2] += dedphi*dyiazib + d2edphi2*dphidyia*dphidzib;
            hessz[ib][2] += dedphi*dziazib + d2edphi2*dphidzia*dphidzib;
            hessx[ic][0] += dedphi*dxiaxic + d2edphi2*dphidxia*dphidxic;
            hessy[ic][0] += dedphi*dyiaxic + d2edphi2*dphidyia*dphidxic;
            hessz[ic][0] += dedphi*dziaxic + d2edphi2*dphidzia*dphidxic;
            hessx[ic][1] += dedphi*dxiayic + d2edphi2*dphidxia*dphidyic;
            hessy[ic][1] += dedphi*dyiayic + d2edphi2*dphidyia*dphidyic;
            hessz[ic][1] += dedphi*dziayic + d2edphi2*dphidzia*dphidyic;
            hessx[ic][2] += dedphi*dxiazic + d2edphi2*dphidxia*dphidzic;
            hessy[ic][2] += dedphi*dyiazic + d2edphi2*dphidyia*dphidzic;
            hessz[ic][2] += dedphi*dziazic + d2edphi2*dphidzia*dphidzic;
            hessx[id][0] += dedphi*dxiaxid + d2edphi2*dphidxia*dphidxid;
            hessy[id][0] += dedphi*dyiaxid + d2edphi2*dphidyia*dphidxid;
            hessz[id][0] += dedphi*dziaxid + d2edphi2*dphidzia*dphidxid;
            hessx[id][1] += dedphi*dxiayid + d2edphi2*dphidxia*dphidyid;
            hessy[id][1] += dedphi*dyiayid + d2edphi2*dphidyia*dphidyid;
            hessz[id][1] += dedphi*dziayid + d2edphi2*dphidzia*dphidyid;
            hessx[id][2] += dedphi*dxiazid + d2edphi2*dphidxia*dphidzid;
            hessy[id][2] += dedphi*dyiazid + d2edphi2*dphidyia*dphidzid;
            hessz[id][2] += dedphi*dziazid + d2edphi2*dphidzia*dphidzid;
         }else if (iatom == ib)
         {
            hessx[ib][0] += dedphi*dxibxib + d2edphi2*dphidxib*dphidxib;
            hessy[ib][0] += dedphi*dxibyib + d2edphi2*dphidxib*dphidyib;
            hessz[ib][0] += dedphi*dxibzib + d2edphi2*dphidxib*dphidzib;
            hessx[ib][1] += dedphi*dxibyib + d2edphi2*dphidxib*dphidyib;
            hessy[ib][1] += dedphi*dyibyib + d2edphi2*dphidyib*dphidyib;
            hessz[ib][1] += dedphi*dyibzib + d2edphi2*dphidyib*dphidzib;
            hessx[ib][2] += dedphi*dxibzib + d2edphi2*dphidxib*dphidzib;
            hessy[ib][2] += dedphi*dyibzib + d2edphi2*dphidyib*dphidzib;
            hessz[ib][2] += dedphi*dzibzib + d2edphi2*dphidzib*dphidzib;
            hessx[ia][0] += dedphi*dxiaxib + d2edphi2*dphidxib*dphidxia;
            hessy[ia][0] += dedphi*dxiayib + d2edphi2*dphidyib*dphidxia;
            hessz[ia][0] += dedphi*dxiazib + d2edphi2*dphidzib*dphidxia;
            hessx[ia][1] += dedphi*dyiaxib + d2edphi2*dphidxib*dphidyia;
            hessy[ia][1] += dedphi*dyiayib + d2edphi2*dphidyib*dphidyia;
            hessz[ia][1] += dedphi*dyiazib + d2edphi2*dphidzib*dphidyia;
            hessx[ia][2] += dedphi*dziaxib + d2edphi2*dphidxib*dphidzia;
            hessy[ia][2] += dedphi*dziayib + d2edphi2*dphidyib*dphidzia;
            hessz[ia][2] += dedphi*dziazib + d2edphi2*dphidzib*dphidzia;
            hessx[ic][0] += dedphi*dxibxic + d2edphi2*dphidxib*dphidxic;
            hessy[ic][0] += dedphi*dyibxic + d2edphi2*dphidyib*dphidxic;
            hessz[ic][0] += dedphi*dzibxic + d2edphi2*dphidzib*dphidxic;
            hessx[ic][1] += dedphi*dxibyic + d2edphi2*dphidxib*dphidyic;
            hessy[ic][1] += dedphi*dyibyic + d2edphi2*dphidyib*dphidyic;
            hessz[ic][1] += dedphi*dzibyic + d2edphi2*dphidzib*dphidyic;
            hessx[ic][2] += dedphi*dxibzic + d2edphi2*dphidxib*dphidzic;
            hessy[ic][2] += dedphi*dyibzic + d2edphi2*dphidyib*dphidzic;
            hessz[ic][2] += dedphi*dzibzic + d2edphi2*dphidzib*dphidzic;
            hessx[id][0] += dedphi*dxibxid + d2edphi2*dphidxib*dphidxid;
            hessy[id][0] += dedphi*dyibxid + d2edphi2*dphidyib*dphidxid;
            hessz[id][0] += dedphi*dzibxid + d2edphi2*dphidzib*dphidxid;
            hessx[id][1] += dedphi*dxibyid + d2edphi2*dphidxib*dphidyid;
            hessy[id][1] += dedphi*dyibyid + d2edphi2*dphidyib*dphidyid;
            hessz[id][1] += dedphi*dzibyid + d2edphi2*dphidzib*dphidyid;
            hessx[id][2] += dedphi*dxibzid + d2edphi2*dphidxib*dphidzid;
            hessy[id][2] += dedphi*dyibzid + d2edphi2*dphidyib*dphidzid;
            hessz[id][2] += dedphi*dzibzid + d2edphi2*dphidzib*dphidzid;
         }else if (iatom == ic)
         {
            hessx[ic][0] += dedphi*dxicxic + d2edphi2*dphidxic*dphidxic;
            hessy[ic][0] += dedphi*dxicyic + d2edphi2*dphidxic*dphidyic;
            hessz[ic][0] += dedphi*dxiczic + d2edphi2*dphidxic*dphidzic;
            hessx[ic][1] += dedphi*dxicyic + d2edphi2*dphidxic*dphidyic;
            hessy[ic][1] += dedphi*dyicyic + d2edphi2*dphidyic*dphidyic;
            hessz[ic][1] += dedphi*dyiczic + d2edphi2*dphidyic*dphidzic;
            hessx[ic][2] += dedphi*dxiczic + d2edphi2*dphidxic*dphidzic;
            hessy[ic][2] += dedphi*dyiczic + d2edphi2*dphidyic*dphidzic;
            hessz[ic][2] += dedphi*dziczic + d2edphi2*dphidzic*dphidzic;
            hessx[ia][0] += dedphi*dxiaxic + d2edphi2*dphidxic*dphidxia;
            hessy[ia][0] += dedphi*dxiayic + d2edphi2*dphidyic*dphidxia;
            hessz[ia][0] += dedphi*dxiazic + d2edphi2*dphidzic*dphidxia;
            hessx[ia][1] += dedphi*dyiaxic + d2edphi2*dphidxic*dphidyia;
            hessy[ia][1] += dedphi*dyiayic + d2edphi2*dphidyic*dphidyia;
            hessz[ia][1] += dedphi*dyiazic + d2edphi2*dphidzic*dphidyia;
            hessx[ia][2] += dedphi*dziaxic + d2edphi2*dphidxic*dphidzia;
            hessy[ia][2] += dedphi*dziayic + d2edphi2*dphidyic*dphidzia;
            hessz[ia][2] += dedphi*dziazic + d2edphi2*dphidzic*dphidzia;
            hessx[ib][0] += dedphi*dxibxic + d2edphi2*dphidxic*dphidxib;
            hessy[ib][0] += dedphi*dxibyic + d2edphi2*dphidyic*dphidxib;
            hessz[ib][0] += dedphi*dxibzic + d2edphi2*dphidzic*dphidxib;
            hessx[ib][1] += dedphi*dyibxic + d2edphi2*dphidxic*dphidyib;
            hessy[ib][1] += dedphi*dyibyic + d2edphi2*dphidyic*dphidyib;
            hessz[ib][1] += dedphi*dyibzic + d2edphi2*dphidzic*dphidyib;
            hessx[ib][2] += dedphi*dzibxic + d2edphi2*dphidxic*dphidzib;
            hessy[ib][2] += dedphi*dzibyic + d2edphi2*dphidyic*dphidzib;
            hessz[ib][2] += dedphi*dzibzic + d2edphi2*dphidzic*dphidzib;
            hessx[id][0] += dedphi*dxicxid + d2edphi2*dphidxic*dphidxid;
            hessy[id][0] += dedphi*dyicxid + d2edphi2*dphidyic*dphidxid;
            hessz[id][0] += dedphi*dzicxid + d2edphi2*dphidzic*dphidxid;
            hessx[id][1] += dedphi*dxicyid + d2edphi2*dphidxic*dphidyid;
            hessy[id][1] += dedphi*dyicyid + d2edphi2*dphidyic*dphidyid;
            hessz[id][1] += dedphi*dzicyid + d2edphi2*dphidzic*dphidyid;
            hessx[id][2] += dedphi*dxiczid + d2edphi2*dphidxic*dphidzid;
            hessy[id][2] += dedphi*dyiczid + d2edphi2*dphidyic*dphidzid;
            hessz[id][2] += dedphi*dziczid + d2edphi2*dphidzic*dphidzid;
         }else if (iatom == id)
         {
            hessx[id][0] += dedphi*dxidxid + d2edphi2*dphidxid*dphidxid;
            hessy[id][0] += dedphi*dxidyid + d2edphi2*dphidxid*dphidyid;
            hessz[id][0] += dedphi*dxidzid + d2edphi2*dphidxid*dphidzid;
            hessx[id][1] += dedphi*dxidyid + d2edphi2*dphidxid*dphidyid;
            hessy[id][1] += dedphi*dyidyid + d2edphi2*dphidyid*dphidyid;
            hessz[id][1] += dedphi*dyidzid + d2edphi2*dphidyid*dphidzid;
            hessx[id][2] += dedphi*dxidzid + d2edphi2*dphidxid*dphidzid;
            hessy[id][2] += dedphi*dyidzid + d2edphi2*dphidyid*dphidzid;
            hessz[id][2] += dedphi*dzidzid + d2edphi2*dphidzid*dphidzid;
            hessx[ia][0] += dedphi*dxiaxid + d2edphi2*dphidxid*dphidxia;
            hessy[ia][0] += dedphi*dxiayid + d2edphi2*dphidyid*dphidxia;
            hessz[ia][0] += dedphi*dxiazid + d2edphi2*dphidzid*dphidxia;
            hessx[ia][1] += dedphi*dyiaxid + d2edphi2*dphidxid*dphidyia;
            hessy[ia][1] += dedphi*dyiayid + d2edphi2*dphidyid*dphidyia;
            hessz[ia][1] += dedphi*dyiazid + d2edphi2*dphidzid*dphidyia;
            hessx[ia][2] += dedphi*dziaxid + d2edphi2*dphidxid*dphidzia;
            hessy[ia][2] += dedphi*dziayid + d2edphi2*dphidyid*dphidzia;
            hessz[ia][2] += dedphi*dziazid + d2edphi2*dphidzid*dphidzia;
            hessx[ib][0] += dedphi*dxibxid + d2edphi2*dphidxid*dphidxib;
            hessy[ib][0] += dedphi*dxibyid + d2edphi2*dphidyid*dphidxib;
            hessz[ib][0] += dedphi*dxibzid + d2edphi2*dphidzid*dphidxib;
            hessx[ib][1] += dedphi*dyibxid + d2edphi2*dphidxid*dphidyib;
            hessy[ib][1] += dedphi*dyibyid + d2edphi2*dphidyid*dphidyib;
            hessz[ib][1] += dedphi*dyibzid + d2edphi2*dphidzid*dphidyib;
            hessx[ib][2] += dedphi*dzibxid + d2edphi2*dphidxid*dphidzib;
            hessy[ib][2] += dedphi*dzibyid + d2edphi2*dphidyid*dphidzib;
            hessz[ib][2] += dedphi*dzibzid + d2edphi2*dphidzid*dphidzib;
            hessx[ic][0] += dedphi*dxicxid + d2edphi2*dphidxid*dphidxic;
            hessy[ic][0] += dedphi*dxicyid + d2edphi2*dphidyid*dphidxic;
            hessz[ic][0] += dedphi*dxiczid + d2edphi2*dphidzid*dphidxic;
            hessx[ic][1] += dedphi*dyicxid + d2edphi2*dphidxid*dphidyic;
            hessy[ic][1] += dedphi*dyicyid + d2edphi2*dphidyid*dphidyic;
            hessz[ic][1] += dedphi*dyiczid + d2edphi2*dphidzid*dphidyic;
            hessx[ic][2] += dedphi*dzicxid + d2edphi2*dphidxid*dphidzic;
            hessy[ic][2] += dedphi*dzicyid + d2edphi2*dphidyid*dphidzic;
            hessz[ic][2] += dedphi*dziczid + d2edphi2*dphidzid*dphidzic;
         }
            }
        }
    }
}
