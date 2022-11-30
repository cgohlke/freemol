#define EXTERN extern

#include "pcwin.h"
#include "fix.h"

EXTERN struct t_minim_values {
        int iprint, ndc, nconst;
        float dielc;
        } minim_values;
int find_bond(int,int);
int find_fixed_bond(int,int);

void egeom(int natom,int *use,double *x,double *y,double *z,double *egeom);
void egeom1(int natom,int *use,double *x,double *y,double *z,double *egeom,double **degeom);
void egeom2(int ,int natom,int *use,double *x,double *y,double *z,float **hessx,float **hessy,float **hessz);

// ============================================
void egeom(int natom,int *use,double *x,double *y,double *z,double *egeom)
{
   int i;
   int ia, ib, ic, id;
   double xr, yr, zr, rik, rik2;
   double xia,yia,zia, xib,yib,zib, xic,yic,zic, xid,yid,zid;
   double xab,yab,zab, xcb,ycb,zcb, rab2,rcb2, rcb;
   double dot, cosine,angle, af1, af2;
   double xba,yba,zba, xdc,ydc,zdc, xt,yt,zt,xu,yu,zu;
   double xtu,ytu,ztu, rt2, ru2, rtru, sine, tf1,tf2;
   double t1, t2, df1, df2,target;
   double dt, dt2, e;

   if (minim_values.iprint)
   {
       fprintf(pcmlogfile,"\nFixed Distance Terms \n");
       fprintf(pcmlogfile,"        At1       At2     R       BLen    Bconst      Eb\n");
   }
// restrained atoms
    for (i=0; i < restrain_atom.natom_restrain; i++)
    {
        ia = restrain_atom.katom_restrain[i];
        if (use[ia])
        {
            xr = yr = zr = 0.0;
            xr = x[ia] - restrain_atom.restrain_position[i][0];
            yr = y[ia] - restrain_atom.restrain_position[i][1];
            zr = z[ia] - restrain_atom.restrain_position[i][2];
            dt2 = (xr*xr + yr*yr + zr*zr);
            e = units.bndunit *restrain_atom.restrain_const[i]*dt2;
            *egeom += e;   
        }
    }
// fixed distances
   for (i=0; i < fx_dist.ndfix; i++)
   {
      ia = fx_dist.kdfix[i][0];
      ib = fx_dist.kdfix[i][1];
      if ( use[ia] || use[ib] )
      {
         xr = x[ia] - x[ib];
         yr = y[ia] - y[ib];
         zr = z[ia] - z[ib];
         rik2 = xr*xr + yr*yr + zr*zr;
         rik = sqrt(rik2);
         df1 = fx_dist.min_dist[i];
         df2 = fx_dist.max_dist[i];
         target = rik;
         if (rik < df1) target = df1;
         if (rik > df2) target = df2;
         dt = rik - target;
         
         dt2 = dt*dt;
         e = units.bndunit *fx_dist.fdconst[i]*dt2;

         *egeom += e;
         if (minim_values.iprint)
           fprintf(pcmlogfile,"Fixed Dist: (%-3d) - (%-3d) %-8.3f %-8.3f %-8.3f %-8.3f = %-8.4f\n",ia
              ,ib, rik, fx_dist.min_dist[i],fx_dist.max_dist[i], fx_dist.fdconst[i],e);
      }
   }
// fixed angle
    for (i=0; i < fx_angle.nafix; i++)
    {
        ia = fx_angle.kafix[i][0];
        ib = fx_angle.kafix[i][1];
        ic = fx_angle.kafix[i][2];
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
            xcb = xic - xib;
            ycb = yic - yib;
            zcb = zic - zib;
            rab2 = xab*xab + yab*yab + zab*zab;
            rcb2 = xcb*xcb + ycb*ycb + zcb*zcb;
            if (rab2 > 0.00001 && rcb2 > 0.00001)
            {
                dot = xab*xcb + yab*ycb + zab*zcb;
                cosine = dot /sqrt(rab2*rcb2);
                if (cosine > 1.0)
                   cosine = 1.0;
                if (cosine < -1.0)
                   cosine = -1.0;
                angle = radian*acos(cosine);
                af1 = fx_angle.min_ang[i];
                af2 = fx_angle.max_ang[i];
                target = angle;
                if (angle < af1) target = af1;
                if (angle > af2) target = af2;
                dt = angle - target;
                dt2 = dt*dt;
                e = units.angunit*fx_angle.faconst[i]*dt2;
                *egeom += e;
                if (minim_values.iprint)
                    fprintf(pcmlogfile,"Fixed Angle: %-3d %-3d %-3d %-8.3f %-8.3f %-8.3f %-8.3f = %-8.4f\n",ia,
                    ib,ic, angle, fx_angle.min_ang[i],fx_angle.max_ang[i], fx_angle.faconst[i],e);
            }
        }
    }
// fixed torsion  
    for (i=0; i < fx_torsion.ntfix; i++)
    {
        ia = fx_torsion.ktfix[i][0];
        ib = fx_torsion.ktfix[i][1];
        ic = fx_torsion.ktfix[i][2];
        id = fx_torsion.ktfix[i][3];
        if (use[ia] || use[ib] || use[ic] || use[id])
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
            if (rtru > 0.0001)
            {
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb);
               cosine = (xt*xu + yt*yu + zt*zu) / rtru;
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru);
               if (cosine > 1.0) cosine = 1.0;
               if (cosine < -1.0) cosine = -1.0;
               angle = radian*acos(cosine);
               if (sine < 0.0) angle = -angle;
               tf1 = fx_torsion.min_tor[i];
               tf2 = fx_torsion.max_tor[i];
               if (angle > tf1 && angle < tf2)
                target = angle;
               else if ( (angle > tf1) && (tf1 > tf2))
                target = angle;
               else if ( (angle < tf2) && (tf1 > tf2))
                target = angle;
               else
               {
                   t1 = angle - tf1;
                   t2 = angle - tf2;
                   if (t1 > 180.0)
                    t1 -= 360.0;
                   else if (t1 < -180.0)
                    t1 += 360.0;
                   if (t2 > 180.0)
                    t2 -= 360.0;
                   else if (t2 < -180.0)
                    t2 += 360.0;
                   if (fabs(t1) < fabs(t2))
                    target = tf1;
                   else
                    target = tf2;
               }
               dt = angle - target;
               if (dt > 180.0) dt -= 360.0;
               if (dt < -180.0) dt += 360.0;
               dt2 = dt*dt;
               e = units.torsunit*dt2*fx_torsion.ftconst[i];
               *egeom += e;
            }
        }
    }   
}
// =====================================================
void egeom1(int natom,int *use,double *x,double *y,double *z,double *egeom,double **degeom)
{
   int i;
   int ia, ib, ic, id;
   double xr, yr, zr, rik, rik2;
   double xia,yia,zia, xib,yib,zib, xic,yic,zic, xid,yid,zid;
   double xab,yab,zab, xcb,ycb,zcb, rab2,rcb2, rcb;
   double dot, cosine,angle, af1, af2;
   double xba,yba,zba, xdc,ydc,zdc, xt,yt,zt,xu,yu,zu;
   double xtu,ytu,ztu, rt2, ru2, rtru, sine, tf1,tf2;
   double xp,yp,zp, rp, xca,yca,zca, xdb,ydb,zdb;
   double t1, t2, df1, df2,target;
   double dt, dt2, e;
   double de, deddt, dedx, dedy,dedz, terma,termc;
   double dedphi,dedxt,dedyt,dedzt,dedxu,dedyu,dedzu;
   double dedxia,dedyia,dedzia;
   double dedxib,dedyib,dedzib;
   double dedxic,dedyic,dedzic;
   double dedxid,dedyid,dedzid;

// restrained atoms
    for (i=0; i < restrain_atom.natom_restrain; i++)
    {
        ia = restrain_atom.katom_restrain[i];
        if (use[ia])
        {
            xr = yr = zr = 0.0;
            xr = x[ia] - restrain_atom.restrain_position[i][0];
            yr = y[ia] - restrain_atom.restrain_position[i][1];
            zr = z[ia] - restrain_atom.restrain_position[i][2];
            rik2 = (xr*xr + yr*yr + zr*zr);
            dt2 = rik2;
            e = units.bndunit *restrain_atom.restrain_const[i]*dt2;
            *egeom += e;

            rik = sqrt(rik2);
            dt = rik;
            if (rik < 0.00001) rik = 1.0;
            de = 2.0*units.bndunit *restrain_atom.restrain_const[i]*dt/rik;
            dedx = de*xr;
            dedy = de*yr;
            dedz = de*zr;
            degeom[ia][0] += dedx;
            degeom[ia][1] += dedy;
            degeom[ia][2] += dedz;
        }
    }
// fixed distances
   for (i=0; i < fx_dist.ndfix; i++)
   {
      ia = fx_dist.kdfix[i][0];
      ib = fx_dist.kdfix[i][1];
      if ( use[ia] || use[ib] )
      {
         xr = x[ia] - x[ib];
         yr = y[ia] - y[ib];
         zr = z[ia] - z[ib];
         rik2 = xr*xr + yr*yr + zr*zr;
         rik = sqrt(rik2);
         df1 = fx_dist.min_dist[i];
         df2 = fx_dist.max_dist[i];
         target = rik;
         if (rik < df1) target = df1;
         if (rik > df2) target = df2;
         dt = rik - target;
         
         dt2 = dt*dt;
         e = units.bndunit *fx_dist.fdconst[i]*dt2;
         if (rik == 0.0) rik = 1.0;
         *egeom += e;

         de = 2.0*units.bndunit *fx_dist.fdconst[i]*dt/rik;
         dedx = de*xr;
         dedy = de*yr;
         dedz = de*zr;
         degeom[ia][0] += dedx;
         degeom[ia][1] += dedy;
         degeom[ia][2] += dedz;
         degeom[ib][0] -= dedx;
         degeom[ib][1] -= dedy;
         degeom[ib][2] -= dedz;
      }
   }
// fixed angle
    for (i=0; i < fx_angle.nafix; i++)
    {
        ia = fx_angle.kafix[i][0];
        ib = fx_angle.kafix[i][1];
        ic = fx_angle.kafix[i][2];
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
            xcb = xic - xib;
            ycb = yic - yib;
            zcb = zic - zib;
            rab2 = xab*xab + yab*yab + zab*zab;
            rcb2 = xcb*xcb + ycb*ycb + zcb*zcb;
            if (rab2 > 0.00001 && rcb2 > 0.00001)
            {
               xp = ycb*zab - zcb*yab;
               yp = zcb*xab - xcb*zab;
               zp = xcb*yab - ycb*xab;
               rp = sqrt(xp*xp + yp*yp + zp*zp);
               if (rp < 0.00001) rp = 0.00001;               
                dot = xab*xcb + yab*ycb + zab*zcb;
                cosine = dot /sqrt(rab2*rcb2);
                if (cosine > 1.0)
                   cosine = 1.0;
                if (cosine < -1.0)
                   cosine = -1.0;
                angle = radian*acos(cosine);
                af1 = fx_angle.min_ang[i];
                af2 = fx_angle.max_ang[i];
                target = angle;
                if (angle < af1) target = af1;
                if (angle > af2) target = af2;
                dt = angle - target;
                dt2 = dt*dt;
                e = units.angunit*fx_angle.faconst[i]*dt2;
                *egeom += e;

                deddt = 2.0*units.angunit*fx_angle.faconst[i]*dt*radian;
                terma = -deddt/(rab2*rp);
                termc = deddt/(rcb2*rp);
                dedxia = terma * (yab*zp-zab*yp);
                dedyia = terma * (zab*xp-xab*zp);
                dedzia = terma * (xab*yp-yab*xp);
                dedxic = termc * (ycb*zp-zcb*yp);
                dedyic = termc * (zcb*xp-xcb*zp);
                dedzic = termc * (xcb*yp-ycb*xp);
                dedxib = -dedxia - dedxic;
                dedyib = -dedyia - dedyic;
                dedzib = -dedzia - dedzic;
                degeom[ia][0] += dedxia;
                degeom[ia][1] += dedyia;
                degeom[ia][2] += dedzia;
                degeom[ib][0] += dedxib;
                degeom[ib][1] += dedyib;
                degeom[ib][2] += dedzib;
                degeom[ic][0] += dedxic;
                degeom[ic][1] += dedyic;
                degeom[ic][2] += dedzic;
            }
        }
    }
// fixed torsion  
    for (i=0; i < fx_torsion.ntfix; i++)
    {
        ia = fx_torsion.ktfix[i][0];
        ib = fx_torsion.ktfix[i][1];
        ic = fx_torsion.ktfix[i][2];
        id = fx_torsion.ktfix[i][3];
        if (use[ia] || use[ib] || use[ic] || use[id])
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
            if (rtru > 0.0001)
            {
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb);
               cosine = (xt*xu + yt*yu + zt*zu) / rtru;
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru);
               if (cosine > 1.0) cosine = 1.0;
               if (cosine < -1.0) cosine = -1.0;
               angle = radian*acos(cosine);
               if (sine < 0.0) angle = -angle;
               tf1 = fx_torsion.min_tor[i];
               tf2 = fx_torsion.max_tor[i];
               if (angle > tf1 && angle < tf2)
                target = angle;
               else if ( (angle > tf1) && (tf1 > tf2))
                target = angle;
               else if ( (angle < tf2) && (tf1 > tf2))
                target = angle;
               else
               {
                   t1 = angle - tf1;
                   t2 = angle - tf2;
                   if (t1 > 180.0)
                    t1 -= 360.0;
                   else if (t1 < -180.0)
                    t1 += 360.0;
                   if (t2 > 180.0)
                    t2 -= 360.0;
                   else if (t2 < -180.0)
                    t2 += 360.0;
                   if (fabs(t1) < fabs(t2))
                    target = tf1;
                   else
                    target = tf2;
               }
               dt = angle - target;
               if (dt > 180.0) dt -= 360.0;
               if (dt < -180.0) dt += 360.0;
               dt2 = dt*dt;
               e = units.torsunit*dt2*fx_torsion.ftconst[i];
               *egeom += e;

               dedphi = 2.0*units.torsunit*fx_torsion.ftconst[i]*dt*radian;
               xca = xic - xia;
               yca = yic - yia;
               zca = zic - zia;
               xdb = xid - xib;
               ydb = yid - yib;
               zdb = zid - zib;
               dedxt = dedphi * (yt*zcb - ycb*zt) / (rt2*rcb);
               dedyt = dedphi * (zt*xcb - zcb*xt) / (rt2*rcb);
               dedzt = dedphi * (xt*ycb - xcb*yt) / (rt2*rcb);
               dedxu = -dedphi * (yu*zcb - ycb*zu) / (ru2*rcb);
               dedyu = -dedphi * (zu*xcb - zcb*xu) / (ru2*rcb);
               dedzu = -dedphi * (xu*ycb - xcb*yu) / (ru2*rcb);

               dedxia = zcb*dedyt - ycb*dedzt;
               dedyia = xcb*dedzt - zcb*dedxt;
               dedzia = ycb*dedxt - xcb*dedyt;
               dedxib = yca*dedzt - zca*dedyt + zdc*dedyu - ydc*dedzu;
               dedyib = zca*dedxt - xca*dedzt + xdc*dedzu - zdc*dedxu;
               dedzib = xca*dedyt - yca*dedxt + ydc*dedxu - xdc*dedyu;
               dedxic = zba*dedyt - yba*dedzt + ydb*dedzu - zdb*dedyu;
               dedyic = xba*dedzt - zba*dedxt + zdb*dedxu - xdb*dedzu;
               dedzic = yba*dedxt - xba*dedyt + xdb*dedyu - ydb*dedxu;
               dedxid = zcb*dedyu - ycb*dedzu;
               dedyid = xcb*dedzu - zcb*dedxu;
               dedzid = ycb*dedxu - xcb*dedyu;
               
               degeom[ia][0] += dedxia;
               degeom[ia][1] += dedyia;
               degeom[ia][2] += dedzia;
               degeom[ib][0] += dedxib;
               degeom[ib][1] += dedyib;
               degeom[ib][2] += dedzib;
               degeom[ic][0] += dedxic;
               degeom[ic][1] += dedyic;
               degeom[ic][2] += dedzic;
               degeom[id][0] += dedxid;
               degeom[id][1] += dedyid;
               degeom[id][2] += dedzid;
            }
        }
    }   
}
// ==========================================
void egeom2(int jatm,int natom,int *use,double *x,double *y,double *z,float **hessx,float **hessy,float **hessz)
{
   int i,j, ia,ib,ic,id;
   double xr, yr, zr, rik, rik2;
   double xia,yia,zia, xib,yib,zib, xic,yic,zic, xid,yid,zid;
   double xab,yab,zab, xcb,ycb,zcb, rab2,rcb2, rcb;
   double dot, cosine,angle, af1, af2;
   double xba,yba,zba, xdc,ydc,zdc, xt,yt,zt,xu,yu,zu;
   double xtu,ytu,ztu, rt2, ru2, rtru, sine, tf1,tf2;
   double xp,yp,zp, rp,rp2, xca,yca,zca, xdb,ydb,zdb;
   double t1, t2, df1, df2,target;
   double dt, dt2;
   double de, deddt, terma,termc;
   double dedphi;
   double d2e[3][3], d2eddt2, d2edphi2;
   double term,termx,termy,termz;
   double xrab,yrab,zrab;
   double xrcb,yrcb,zrcb;
   double xabp,yabp,zabp;
   double xcbp,ycbp,zcbp;
      double ddtdxia,ddtdyia,ddtdzia;
      double ddtdxib,ddtdyib,ddtdzib;
      double ddtdxic,ddtdyic,ddtdzic;
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
      double dxiaxia,dyiayia,dziazia;
      double dxibxib,dyibyib,dzibzib;
      double dxicxic,dyicyic,dziczic;
      double dxidxid,dyidyid,dzidzid;
      double dxiayia,dxiazia,dyiazia;
      double dxibyib,dxibzib,dyibzib;
      double dxicyic,dxiczic,dyiczic;
      double dxidyid,dxidzid,dyidzid;
      double dxiaxib,dxiayib,dxiazib;
      double dyiaxib,dyiayib,dyiazib;
      double dziaxib,dziayib,dziazib;
      double dxiaxic,dxiayic,dxiazic;
      double dyiaxic,dyiayic,dyiazic;
      double dziaxic,dziayic,dziazic;
      double dxiaxid,dxiayid,dxiazid;
      double dyiaxid,dyiayid,dyiazid;
      double dziaxid,dziayid,dziazid;
      double dxibxia,dxibyia,dxibzia;
      double dyibxia,dyibyia,dyibzia;
      double dzibxia,dzibyia,dzibzia;
      double dxibxic,dxibyic,dxibzic;
      double dyibxic,dyibyic,dyibzic;
      double dzibxic,dzibyic,dzibzic;
      double dxibxid,dxibyid,dxibzid;
      double dyibxid,dyibyid,dyibzid;
      double dzibxid,dzibyid,dzibzid;
      double dxicxid,dxicyid,dxiczid;
      double dyicxid,dyicyid,dyiczid;
      double dzicxid,dzicyid,dziczid;
 
// restrained atoms
    for (i=0; i < restrain_atom.natom_restrain; i++)
    {
        ia = restrain_atom.katom_restrain[i];
        if (ia == jatm && use[ia])
        {
            xr = yr = zr = 0.0;
            xr = x[ia] - restrain_atom.restrain_position[i][0];
            yr = y[ia] - restrain_atom.restrain_position[i][1];
            zr = z[ia] - restrain_atom.restrain_position[i][2];
            rik2 = (xr*xr + yr*yr + zr*zr);
            dt2 = rik2;
            rik = sqrt(rik2);
            dt = rik;
            deddt = 2.0*units.bndunit *restrain_atom.restrain_const[i];
            if (dt == 0.0) deddt = 0.0;
            if (rik < 0.00001)
            {
                de = deddt;
                term = 0.0;
            } else
            {
                de = deddt *dt/rik;
                term = (deddt-de)/rik2;
            }
            de = deddt*dt/rik;
            term = (deddt-de) / rik2;
            termx = term * xr;
            termy = term * yr;
            termz = term * zr;
         
            d2e[0][0] = termx*xr + de;
            d2e[0][1] = termx*yr;
            d2e[0][2] = termx*zr;
            d2e[1][0] = d2e[0][0];
            d2e[1][1] = termy*yr + de;
            d2e[1][2] = termy*zr;
            d2e[2][0] = d2e[0][2];
            d2e[2][1] = d2e[1][2];
            d2e[2][2] = termz*zr + de;

            for (j=0; j < 3; j++)
            {
               hessx[ia][j] += d2e[j][0];
               hessy[ia][j] += d2e[j][1];
               hessz[ia][j] += d2e[j][2];
            }
        }
    }
// fixed distance
   for (i=0; i < fx_dist.ndfix; i++)
   {
      ia = fx_dist.kdfix[i][0];
      ib = fx_dist.kdfix[i][1];
      if ( jatm == ia || jatm == ib )
      {
         if (jatm == ib)
         {
             ib = ia;
             ia = jatm;
         }
         xr = x[ia] - x[ib];
         yr = y[ia] - y[ib];
         zr = z[ia] - z[ib];
         rik2 = xr*xr + yr*yr + zr*zr;
         rik = sqrt(rik2);
         df1 = fx_dist.min_dist[i];
         df2 = fx_dist.max_dist[i];
         target = rik;
         if (rik < df1) target = df1;
         if (rik > df2) target = df2;
         dt = rik - target;
         dt2 = dt*dt;
         deddt = 2.0*units.bndunit *fx_dist.fdconst[i];
         if (dt == 0.0) deddt = 0.0;
         if (rik < 0.00001)
         {
             rik = 0.0001;
             rik2 = rik*rik;
         }
         de = deddt*dt/rik;
         term = (deddt-de) / rik2;
         termx = term * xr;
         termy = term * yr;
         termz = term * zr;
         
         d2e[0][0] = termx*xr + de;
         d2e[0][1] = termx*yr;
         d2e[0][2] = termx*zr;
         d2e[1][0] = d2e[0][0];
         d2e[1][1] = termy*yr + de;
         d2e[1][2] = termy*zr;
         d2e[2][0] = d2e[0][2];
         d2e[2][1] = d2e[1][2];
         d2e[2][2] = termz*zr + de;

         for (j=0; j < 3; j++)
         {
               hessx[ia][j] += d2e[j][0];
               hessy[ia][j] += d2e[j][1];
               hessz[ia][j] += d2e[j][2];
               hessx[ib][j] -= d2e[j][0];
               hessy[ib][j] -= d2e[j][1];
               hessz[ib][j] -= d2e[j][2];
         }
      }
   } 
// fixed angle
    for (i=0; i < fx_angle.nafix; i++)
    {
        ia = fx_angle.kafix[i][0];
        ib = fx_angle.kafix[i][1];
        ic = fx_angle.kafix[i][2];
        if (ia == jatm || ib == jatm || ic == jatm)
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
            xcb = xic - xib;
            ycb = yic - yib;
            zcb = zic - zib;
            rab2 = xab*xab + yab*yab + zab*zab;
            rcb2 = xcb*xcb + ycb*ycb + zcb*zcb;
            if (rab2 > 0.00001 && rcb2 > 0.00001)
            {
               xp = ycb*zab - zcb*yab;
               yp = zcb*xab - xcb*zab;
               zp = xcb*yab - ycb*xab;
               rp = sqrt(xp*xp + yp*yp + zp*zp);
               if (rp < 0.00001) rp = 0.00001;               
                dot = xab*xcb + yab*ycb + zab*zcb;
                cosine = dot /sqrt(rab2*rcb2);
                if (cosine > 1.0)
                   cosine = 1.0;
                if (cosine < -1.0)
                   cosine = -1.0;
                angle = radian*acos(cosine);
                af1 = fx_angle.min_ang[i];
                af2 = fx_angle.max_ang[i];
                target = angle;
                if (angle < af1) target = af1;
                if (angle > af2) target = af2;
                dt = angle - target;
                dt2 = dt*dt;
                deddt = 2.0*units.angunit*fx_angle.faconst[i]*dt*radian;
                d2eddt2 = 2.0*units.angunit*fx_angle.faconst[i]*radian*radian;
                if (dt == 0.0) d2edphi2 = 0.0;
                
                terma = -1.0/(rab2*rp);
                termc = 1.0/(rcb2*rp);
                ddtdxia = terma * (yab*zp-zab*yp);
                ddtdyia = terma * (zab*xp-xab*zp);
                ddtdzia = terma * (xab*yp-yab*xp);
                ddtdxic = termc * (ycb*zp-zcb*yp);
                ddtdyic = termc * (zcb*xp-xcb*zp);
                ddtdzic = termc * (xcb*yp-ycb*xp);
                ddtdxib = -ddtdxia - ddtdxic;
                ddtdyib = -ddtdyia - ddtdyic;
                ddtdzib = -ddtdzia - ddtdzic;
                xrab = 2.00 * xab / rab2;
                yrab = 2.00 * yab / rab2;
                zrab = 2.00 * zab / rab2;
                xrcb = 2.00 * xcb / rcb2;
                yrcb = 2.00 * ycb / rcb2;
                zrcb = 2.00 * zcb / rcb2;
                rp2 = 1.00 / (rp*rp);
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

                  if (ia == jatm)
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
                  } else if (ib == jatm)
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
                  }else if (ic == jatm)
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
// fixed torsion
    for (i=0; i < fx_torsion.ntfix; i++)
    {
        ia = fx_torsion.ktfix[i][0];
        ib = fx_torsion.ktfix[i][1];
        ic = fx_torsion.ktfix[i][2];
        id = fx_torsion.ktfix[i][3];
        if (ia == jatm || ib == jatm  || ic == jatm || id == jatm )
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
            if (rtru > 0.0001)
            {
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb);
               cosine = (xt*xu + yt*yu + zt*zu) / rtru;
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru);
               if (cosine > 1.0) cosine = 1.0;
               if (cosine < -1.0) cosine = -1.0;
               angle = radian*acos(cosine);
               if (sine < 0.0) angle = -angle;
               tf1 = fx_torsion.min_tor[i];
               tf2 = fx_torsion.max_tor[i];
               if (angle > tf1 && angle < tf2)
                target = angle;
               else if ( (angle > tf1) && (tf1 > tf2))
                target = angle;
               else if ( (angle < tf2) && (tf1 > tf2))
                target = angle;
               else
               {
                   t1 = angle - tf1;
                   t2 = angle - tf2;
                   if (t1 > 180.0)
                    t1 -= 360.0;
                   else if (t1 < -180.0)
                    t1 += 360.0;
                   if (t2 > 180.0)
                    t2 -= 360.0;
                   else if (t2 < -180.0)
                    t2 += 360.0;
                   if (fabs(t1) < fabs(t2))
                    target = tf1;
                   else
                    target = tf2;
               }
               dt = angle - target;
               if (dt > 180.0) dt -= 360.0;
               if (dt < -180.0) dt += 360.0;
               dt2 = dt*dt;
               dedphi = 2.0*units.torsunit*fx_torsion.ftconst[i]*dt;
               d2edphi2 = 2.0*units.torsunit*fx_torsion.ftconst[i]*radian*radian;
               if (dt < 0.000001) d2edphi2 = 0.0;

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
               rcbxt = -2.00 * rcb * dphidxt;
               rcbyt = -2.00 * rcb * dphidyt;
               rcbzt = -2.00 * rcb * dphidzt;
               rcbt2 = rcb * rt2;
               rcbxu = 2.00 * rcb * dphidxu;
               rcbyu = 2.00 * rcb * dphidyu;
               rcbzu = 2.00 * rcb * dphidzu;
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
               dxiaxic = rcbxt*dphidxict + xcb*xt/rcbt2;
               dxiayic = rcbxt*dphidyict - dphidzt - (xba*zcb*xcb+zba*yzcb2)/rcbt2;
               dxiazic = rcbxt*dphidzict + dphidyt + (xba*ycb*xcb+yba*yzcb2)/rcbt2;
               dxiaxid = 0.00;
               dxiayid = 0.00;
               dxiazid = 0.00;
               dyiayia = rcbyt*dphidyia;
               dyiazia = rcbyt*dphidzia - xcb*rcb/rt2;
               dyiaxib = rcbyt*dphidxibt - dphidzt - (yca*zcb*ycb+zca*xzcb2)/rcbt2;
               dyiaxic = rcbyt*dphidxict + dphidzt + (yba*zcb*ycb+zba*xzcb2)/rcbt2;
               dyiayic = rcbyt*dphidyict + ycb*yt/rcbt2;
               dyiazic = rcbyt*dphidzict - dphidxt - (yba*xcb*ycb+xba*xzcb2)/rcbt2;
               dyiaxid = 0.00;
               dyiayid = 0.00;
               dyiazid = 0.00;
               dziazia = rcbzt*dphidzia;
               dziaxib = rcbzt*dphidxibt + dphidyt + (zca*ycb*zcb+yca*xycb2)/rcbt2;
               dziayib = rcbzt*dphidyibt - dphidxt - (zca*xcb*zcb+xca*xycb2)/rcbt2;
               dziaxic = rcbzt*dphidxict - dphidyt - (zba*ycb*zcb+yba*xycb2)/rcbt2;
               dziayic = rcbzt*dphidyict + dphidxt + (zba*xcb*zcb+xba*xycb2)/rcbt2;
               dziazic = rcbzt*dphidzict + zcb*zt/rcbt2;
               dziaxid = 0.00;
               dziayid = 0.00;
               dziazid = 0.00;
               dxibxic = -xcb*dphidxib/(rcb*rcb) - (yca*(zba*xcb+yt)-zca*(yba*xcb-zt))/rcbt2
                  - 2.00*(yt*zba-yba*zt)*dphidxibt/rt2 - (zdc*(ydb*xcb+zu)-ydc*(zdb*xcb-yu))/rcbu2
                  + 2.00*(yu*zdb-ydb*zu)*dphidxibu/ru2;
               dxibyic = -ycb*dphidxib/(rcb*rcb) + dphidzt + dphidzu
                  - (yca*(zba*ycb-xt)+zca*(xba*xcb+zcb*zba))/rcbt2
                  - 2.00*(zt*xba-zba*xt)*dphidxibt/rt2
                  + (zdc*(xdb*xcb+zcb*zdb)+ydc*(zdb*ycb+xu))/rcbu2
                  + 2.00*(zu*xdb-zdb*xu)*dphidxibu/ru2;
               dxibxid = rcbxu*dphidxibu + xcb*xu/rcbu2;
               dxibyid = rcbyu*dphidxibu - dphidzu - (ydc*zcb*ycb+zdc*xzcb2)/rcbu2;
               dxibzid = rcbzu*dphidxibu + dphidyu + (zdc*ycb*zcb+ydc*xycb2)/rcbu2;
               dyibzib = ycb*dphidzib/(rcb*rcb) - (xca*(xca*xcb+zcb*zca)+yca*(ycb*xca+zt))/rcbt2
                  - 2.00*(xt*zca-xca*zt)*dphidzibt/rt2
                  + (ydc*(xdc*ycb-zu)+xdc*(xdc*xcb+zcb*zdc))/rcbu2
                  + 2.00*(xu*zdc-xdc*zu)*dphidzibu/ru2;
               dyibxic = -xcb*dphidyib/(rcb*rcb) - dphidzt - dphidzu + (xca*(zba*xcb+yt)+zca*(zba*zcb+ycb*yba))/rcbt2
                  - 2.00*(yt*zba-yba*zt)*dphidyibt/rt2
                  - (zdc*(zdb*zcb+ycb*ydb)+xdc*(zdb*xcb-yu))/rcbu2
                  + 2.00*(yu*zdb-ydb*zu)*dphidyibu/ru2;
               dyibyic = -ycb*dphidyib/(rcb*rcb)  - (zca*(xba*ycb+zt)-xca*(zba*ycb-xt))/rcbt2
                  - 2.00*(zt*xba-zba*xt)*dphidyibt/rt2
                  - (xdc*(zdb*ycb+xu)-zdc*(xdb*ycb-zu))/rcbu2
                  + 2.00*(zu*xdb-zdb*xu)*dphidyibu/ru2;
               dyibxid = rcbxu*dphidyibu + dphidzu  + (xdc*zcb*xcb+zdc*yzcb2)/rcbu2;
               dyibyid = rcbyu*dphidyibu + ycb*yu/rcbu2;
               dyibzid = rcbzu*dphidyibu - dphidxu - (zdc*xcb*zcb+xdc*xycb2)/rcbu2;
               dzibxic = -xcb*dphidzib/(rcb*rcb) + dphidyt + dphidyu  - (xca*(yba*xcb-zt)+yca*(zba*zcb+ycb*yba))/rcbt2
                  - 2.00*(yt*zba-yba*zt)*dphidzibt/rt2
                  + (ydc*(zdb*zcb+ycb*ydb)+xdc*(ydb*xcb+zu))/rcbu2
                  + 2.00*(yu*zdb-ydb*zu)*dphidzibu/ru2;
               dzibzic = -zcb*dphidzib/(rcb*rcb)  - (xca*(yba*zcb+xt)-yca*(xba*zcb-yt))/rcbt2
                  - 2.00*(xt*yba-xba*yt)*dphidzibt/rt2
                  - (ydc*(xdb*zcb+yu)-xdc*(ydb*zcb-xu))/rcbu2
                  + 2.00*(xu*ydb-xdb*yu)*dphidzibu/ru2;
               dzibxid = rcbxu*dphidzibu - dphidyu - (xdc*ycb*xcb+ydc*yzcb2)/rcbu2;
               dzibyid = rcbyu*dphidzibu + dphidxu + (ydc*xcb*ycb+xdc*xzcb2)/rcbu2;
               dzibzid = rcbzu*dphidzibu + zcb*zu/rcbu2;
               dxicxid = rcbxu*dphidxicu - xcb*(zdb*ycb-ydb*zcb)/rcbu2;
               dxicyid = rcbyu*dphidxicu + dphidzu + (ydb*zcb*ycb+zdb*xzcb2)/rcbu2;
               dxiczid = rcbzu*dphidxicu - dphidyu - (zdb*ycb*zcb+ydb*xycb2)/rcbu2;
               dyicxid = rcbxu*dphidyicu - dphidzu  - (xdb*zcb*xcb+zdb*yzcb2)/rcbu2;
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

               dxiaxib = -dxiaxia - dxiaxic - dxiaxid;
               dxiayib = -dxiayia - dxiayic - dxiayid;
               dxiazib = -dxiazia - dxiazic - dxiazid;
               dyiayib = -dyiayia - dyiayic - dyiayid;
               dyiazib = -dyiazia - dyiazic - dyiazid;
               dziazib = -dziazia - dziazic - dziazid;
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
         if (jatm == ia)
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
         }else if (jatm == ib)
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
         }else if (jatm == ic)
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
         }else if (jatm == id)
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
