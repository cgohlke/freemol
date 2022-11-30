#define EXTERN extern

#include "pcwin.h"
#include "pcmod.h"
#include "utility.h"
#include "asnsym.h"

double TOLER = 0.002;
double TOL2 = 0.00002;

double SIGN(double,double);
static void vibsym(int);
static int tstc3(double **,double **,int,int,int,double *);
static int tstc4(double **,double **,int,int,int,double *);
static int tstc5(double **,double **,int,int,int,double *);
static void secmom(double **, double *, double *, double *, int *, double (*)[3]);
static void orients(double **, double **,double *, double *,int,int);
static void tcopy(double **, double **,int);
static int findcn(double **,double **,double *, double *, int *);
static int findc2(double **,double **, double **, int *);
static int findv(double **,double **, double **, int *);
static int findsn(double **,double **, double **, int *);
static void tform(double (*)[3], double **, double **,int);
static void reflect(double **, double **,int);
static int tcomp(double **,double **);
static void circset(double **, double **,int,int *,int *);
static void rotates(double **a,double **b,int,double);
static void sphere(double **,double **, double **,int *, int *);
static void sym_invert(double **, double **);
static void triang(double **,int,int,int,double *, double *,double *,double *,double *,double *);
static void zero_array(double **,int);
static void hqrii(double *,int,int,double *,double (*)[3]);
static void sphset(double **, double **,int *, int *,int *);

// ========================
void zero_array(double **a,int num)
{
    int i;
    for (i=1; i <= num; i++)
    {
        a[i][0] = 0.0;
        a[i][1] = 0.0;
        a[i][2] = 0.0;
    }
}
void triang(double **a,int iat,int jat,int kat,double *alpha, double *beta,
            double *gamma,double *dij, double *dik, double *djk)
{
    double xi,yi,zi,xj,yj,zj,xk,yk,zk,dotjk,dotik,dotij;

    xi = a[iat][0];
    yi = a[iat][1];
    zi = a[iat][2];
    xj = a[jat][0];
    yj = a[jat][1];
    zj = a[jat][2];
    xk = a[kat][0];
    yk = a[kat][1];
    zk = a[kat][2];
    *dij = sqrt( (xi-xj)*(xi-xj) + (yi-yj)*(yi-yj) + (zi-zj)*(zi-zj) );
    *dik = sqrt( (xi-xk)*(xi-xk) + (yi-yk)*(yi-yk) + (zi-zk)*(zi-zk) );
    *djk = sqrt( (xj-xk)*(xj-xk) + (yj-yk)*(yj-yk) + (zj-zk)*(zj-zk) );
    dotjk = (xj-xi)*(xk-xi) + (yj-yi)*(yk-yi) + (zj-zi)*(zk-zi);
    dotik = (xi-xj)*(xk-xj) + (yi-yj)*(yk-yj) + (zi-zj)*(zk-zj);
    dotij = (xi-xk)*(xj-xk) + (yi-yk)*(yj-yk) + (zi-zk)*(zj-zk);
    *alpha = acos(dotjk/(*dij**dik));
    *beta  = acos(dotik/(*dij**djk));
    *gamma = acos(dotij/(*dik**djk));
    
}
/* ================================================= */
int tstc3(double **a,double **b,int iat,int jat,int kat,double *center)
{
    int itst, iz;
    double p1[3];
    double phi3,theta3,fact3;
    double alpha,beta,gamma,dij,dik,djk;
    double px,py,pz;

    p1[0] = p1[1] = p1[2] = 0.0;
    itst = 0;
    phi3 = (8.0/3.0)*atan(1.0);
    theta3 = 0.5*phi3;
    fact3 = 2.0/3.0;

    triang(a,iat,jat,kat,&alpha,&beta,&gamma,&dij,&dik,&djk);

    if (fabs(alpha-theta3) < TOLER && fabs(dij-dik) < TOLER)
    {
        px = 0.5*(a[jat][0] + a[kat][0]);
        py = 0.5*(a[jat][1] + a[kat][1]);
        pz = 0.5*(a[jat][2] + a[kat][2]);
        center[0] = a[iat][0] + fact3*(px - a[iat][0]);
        center[1] = a[iat][1] + fact3*(py - a[iat][1]);
        center[2] = a[iat][2] + fact3*(pz - a[iat][2]);
    
        orients(a,b,p1,center,3,natom);
        tcopy(b,a,natom);
        rotates(a,b,3,phi3);
        iz = tcomp(a,b);
        return (iz);
    }
    return FALSE;
}
int tstc4(double **a,double **b,int iat,int jat,int kat,double *center)
{
    int itst, iz;
    double p1[3];
    double alpha,beta,gamma,dij,dik,djk;
    double halfpi;

    p1[0] = p1[1] = p1[2] = 0.0;
    itst = 0;
    halfpi = 2.0*atan(1.0);

    triang(a,iat,jat,kat,&alpha,&beta,&gamma,&dij,&dik,&djk);

    if (fabs(alpha-halfpi) < TOLER && fabs(dij-dik) < TOLER)
    {
        center[0] = 0.5*(a[jat][0] + a[kat][0]);
        center[1] = 0.5*(a[jat][1] + a[kat][1]);
        center[2] = 0.5*(a[jat][2] + a[kat][2]);
    } else if (fabs(beta-halfpi) < TOLER && fabs(dij-djk) < TOLER)
    {
        center[0] = 0.5*(a[iat][0] + a[kat][0]);
        center[1] = 0.5*(a[iat][1] + a[kat][1]);
        center[2] = 0.5*(a[iat][2] + a[kat][2]);
    } else if (fabs(gamma-halfpi) < TOLER && fabs(dik-djk) < TOLER)
    {
        center[0] = 0.5*(a[jat][0] + a[iat][0]);
        center[1] = 0.5*(a[jat][1] + a[iat][1]);
        center[2] = 0.5*(a[jat][2] + a[iat][2]);
    } else if (fabs(dij-dik) < TOLER && fabs(dij-djk) < TOLER && TOLER && fabs(dik-djk))
    {
        center[0] = a[iat][0];
        center[1] = a[iat][1];
        center[2] = a[iat][2];
    } else
       return FALSE;

    orients(a,b,p1,center,3,natom);
    tcopy(b,a,natom);
    rotates(a,b,3,halfpi);
    iz = tcomp(a,b);
    return (iz);

}
int tstc5(double **a,double **b,int iat,int jat,int kat,double *center)
{
    int itst, iz;
    double p1[3];
    double piovr4, phi5,theta5,fact5;
    double alpha,beta,gamma,dij,dik,djk;
    double px,py,pz;

    p1[0] = p1[1] = p1[2] = 0.0;
    itst = 0;
    piovr4 = atan(1.0);
    phi5 = 1.6*piovr4;
    theta5 = 2.4*piovr4;
    fact5 = 1.0/(2.0*sin(0.8*piovr4)*sin(0.8*piovr4));

    triang(a,iat,jat,kat,&alpha,&beta,&gamma,&dij,&dik,&djk);

    if (fabs(alpha-theta5) < TOLER && fabs(dij-dik) < TOLER)
    {
        px = 0.5*(a[jat][0] + a[kat][0]);
        py = 0.5*(a[jat][1] + a[kat][1]);
        pz = 0.5*(a[jat][2] + a[kat][2]);
        center[0] = a[iat][0] + fact5*(px - a[iat][0]);
        center[1] = a[iat][1] + fact5*(py - a[iat][1]);
        center[2] = a[iat][2] + fact5*(pz - a[iat][2]);
    } else if (fabs(beta-theta5) < TOLER && fabs(dij-djk) < TOLER)
    {
        px = 0.5*(a[iat][0] + a[kat][0]);
        py = 0.5*(a[iat][1] + a[kat][1]);
        pz = 0.5*(a[iat][2] + a[kat][2]);
        center[0] = a[jat][0] + fact5*(px - a[jat][0]);
        center[1] = a[jat][1] + fact5*(py - a[jat][1]);
        center[2] = a[jat][2] + fact5*(pz - a[jat][2]);
    } else if (fabs(gamma-theta5) < TOLER && fabs(dik-djk) < TOLER)
    {
        px = 0.5*(a[jat][0] + a[iat][0]);
        py = 0.5*(a[jat][1] + a[iat][1]);
        pz = 0.5*(a[jat][2] + a[iat][2]);
        center[0] = a[kat][0] + fact5*(px - a[kat][0]);
        center[1] = a[kat][1] + fact5*(py - a[kat][1]);
        center[2] = a[kat][2] + fact5*(pz - a[kat][2]);
    } else
       return FALSE;
    orients(a,b,p1,center,3,natom);
    tcopy(b,a,natom);
    rotates(a,b,3,phi5);
    iz = tcomp(a,b);
    return (iz);
}
/* ======================================  */
void sphset(double **a, double **aset,int *npop, int *nset,int *numset)
{
    int i,j1,j,iset, num;
    double xtmp,ytmp,ztmp;
    
    for (i=1; i <= natom; i++)
    {
        aset[i][0] = i;
        aset[i][1] = atom.atomwt[i];
        xtmp = a[i][0];
        ytmp = a[i][1];
        ztmp = a[i][2];
        aset[i][2] = sqrt(xtmp*xtmp+ytmp*ytmp+ztmp*ztmp);
    }
    iset = 0;
    for (i=1; i < natom; i++)
    {
        num = 0;
        if (aset[i][0] != 0.0)
        {
            iset++;
            num++;
            nset[i] = iset;
            aset[i][0] = 0.0;
            j1 = i + 1;
            for (j=j1; j <= natom; j++)
            {
                if (aset[j][0] != 0.0)
                {
                    if (fabs(aset[j][1]-aset[i][1]) <= TOLER &&
                        fabs(aset[j][2]-aset[i][2]) <= TOLER)
                    {
                        num++;
                        nset[j] = iset;
                        aset[j][0] = 0.0;
                    }
                }
            }
            npop[iset] = num;
        }
    }
    *numset = iset; 
}
void sphere(double **a, double **b, double **c, int *nset, int *norder)
{
    int i,j,k,keyset, numset, minnum;
    int mpop, num3, i2, j1, j2, k1, iat=0,jat,kat;
    int iz,iz1,iz2,key,ihop,ixyz;
    double savez,curz,x,y,z,cmax,theta,tst;
    int *npop;
    double halfpi,pi;
    double p1[3],center[3],save1[3],save2[3];

    halfpi = 2.0*atan(1.0);
    pi = 2.0*halfpi;
    p1[0] = p1[1] = p1[2] = 0.0;
    save2[0] = save2[1] = save2[2] = 0.0;
    savez = 0.0;
    jat = 0;
    keyset = 0;
    npop = ivector(0,natom+1);
    sphset(a,c,npop,nset,&numset);

    minnum = natom;
    *norder = 0;
    for (i=1; i <= numset; i++)
    {
        j = npop[i];
        if (npop[i] <= minnum && npop[i] > 3)
        {
            minnum = npop[i];
            keyset = i;
        }
    }

    key = natom;
    mpop = 0;
    for (i=1; i <= natom; i++)
    {
        if (nset[i] == keyset)
        {
            mpop++;
            if (key < i)
               key = i;
            nset[mpop] = i;
        }
    }

    if (mpop < 3)
    {
        *norder = 0;
        free_ivector(npop,0,natom+1);
        return;
    }
    num3 = 0;
    i2 = mpop-2;
    j2 = mpop-1;
    for (i=1; i <= i2; i++)
    {
        iat = nset[i];
        j1 = i + 1;
        for (j=j1; j <= j2; j++)
        {
            jat = nset[j];
            k1 = j+1;
            for (k=k1; k <= mpop; k++)
            {
                kat = nset[k];
                if (*norder <= 3)
                {
                    tcopy(a,c,natom);
                    iz = tstc5(c,b,iat,jat,kat,center);
                    if (iz == TRUE)
                    {
                        *norder = 5;
                        savez = fabs(c[key][2]);
                        save1[0] = center[0];
                        save1[1] = center[1];
                        save1[2] = center[2];
                    }else
                    {
                        tcopy(a,c,natom);
                        zero_array(b,natom);
                        iz1 = tstc4(c,b,iat,jat,kat,center);
                        if (iz1 == TRUE)
                        {
                            *norder = 4;
                            tcopy(c,a,natom);
                        } else
                        {
                            if (num3 != 2)
                            {
                                tcopy(a,c,natom);
                                iz2 = tstc3(c,b,iat,jat,kat,center);
                                if (iz2 == TRUE)
                                {
                                    ihop = num3+1;
                                    if (ihop == 1)
                                    {
                                        *norder = 3;
                                        num3 = 1;
                                        save1[0] = center[0];
                                        save1[1] = center[1];
                                        save1[2] = center[2];
                                    } else if (ihop == 2)
                                    {
                                        num3 = 2;
                                        save2[0] = center[0];
                                        save2[1] = center[1];
                                        save2[2] = center[2];
                                    }
                                }
                            }
                        }
                    }
                } else if (*norder == 5)
                {
                    tcopy(a,c,natom);
                    iz = tstc5(c,b,iat,jat,kat,center);
                    if (iz == TRUE)
                    {
                        curz = fabs(c[key][2]);
                        if (fabs(curz-savez) > TOLER && (savez < curz ))
                        {
                            savez = curz;
                            save1[0] = center[0];
                            save1[1] = center[1];
                            save1[2] = center[2];
                        }
                    }
                } else if (*norder == 4)
                {
                    tcopy(a,c,natom);
                    iz = tstc4(c,b,iat,jat,kat,center);
                    if (iz == TRUE)
                    {
                        if (fabs(center[2]) < TOLER)
                        {
                            x = center[0];
                            y = center[1];
                            theta = halfpi;
                            if (fabs(y) > TOLER) theta = atan(x/y);
                            rotates(c,a,3,theta);
                            goto L_5;
                        }
                    }
                }
            }
        }
    }
L_5:
    if (*norder == 0)
    {
        center[0] = 0.5*(a[iat][0] + a[jat][0]);
        center[1] = 0.5*(a[iat][1] + a[jat][1]);
        center[2] = 0.5*(a[iat][2] + a[jat][2]);
        orients(a,b,p1,center,3,natom);
        x = b[iat][0];
        y = b[iat][1];
        theta = halfpi;
        if (fabs(y) > TOLER) theta = atan(x/y);
        rotates(b,a,3,theta);
    } else if (*norder == 3)
    {
        center[0] = 0.5*(save1[0] + save2[0]);
        center[1] = 0.5*(save1[1] + save2[1]);
        center[2] = 0.5*(save1[2] + save2[2]);
        orients(a,b,p1,center,3,natom);
        tcopy(b,a,natom);
        b[0][0] = save1[0];
        b[0][1] = save1[1];
        b[0][2] = save1[2];
        b[1][0] = save2[0];
        b[1][1] = save2[1];
        b[1][2] = save2[2];
        orients(b,c,p1,center,3,2);
        tcopy(c,b,2);
        save1[0] = b[0][0];
        save1[1] = b[0][1];
        save1[2] = b[0][2];
        save2[0] = b[1][0];
        save2[1] = b[1][1];
        save2[2] = b[1][2];
        x = 0.5*(save1[0]-save2[1]);
        y = 0.5*(save1[1]+save2[0]);
        theta = halfpi;
        if (fabs(y) > TOLER) theta = atan(x/y);
        rotates(a,b,3,theta);
        tcopy(b,a,natom);
        
    } else if (*norder == 4)
    {
        x = a[key][0];
        y = a[key][1];
        z = a[key][2];
        cmax = fabs(z);
        ixyz = 3;
        for (i=1; i <= 2; i++)
        {
            tst = fabs(a[key][i]);
            if (fabs(tst-cmax) > TOLER && (cmax <= tst))
            {
                ixyz = i;
                cmax = tst;
            }
        }
        if (ixyz < 3)
        {
            if (ixyz == 2)
               ixyz = 1;
            else
               ixyz = 2;
            rotates(a,b,ixyz,halfpi);
            tcopy(b,a,natom);
        }
        if (a[key][3] < 0.0)
        {
            rotates(a,b,1,pi);
            tcopy(b,a,natom);
        }
    } else if (*norder == 5)
    {
        orients(a,b,p1,save1,3,natom);
        tcopy(b,a,natom);
        if (a[key][3] < 0.0)
        {
            rotates(a,b,1,pi);
            tcopy(b,a,natom);
        }
    }
    
        
    free_ivector(npop,0,natom+1);
}
void sym_invert(double **a, double **b)
{
    int i,j;
    double t[3][3];
    for (i=0; i < 3; i++)
    {
        for (j=0; j <3; j++)
          t[i][j] = 0.0;
    }
    t[0][0] = -1.0;
    t[1][1] = -1.0;
    t[2][2] = -1.0;
    tform(t,a,b,natom);
}
void reflect(double **a, double **b, int iaxis)
{
    int i,j;
    double t[3][3];
    for (i=0; i < 3; i++)
    {
        for (j=0; j <3; j++)
          t[i][j] = 0.0;
        t[i][i] = 1.0;
    }
    t[iaxis-1][iaxis-1] = -1.0;
    tform(t,a,b,natom);
}
/* ====================================  */
int tcomp(double **a, double **b)
{
    int i,j, k;
    int *nrs;
    double dsum, x1, x2;

    nrs = ivector(0,natom+1);
    for (i=0; i <= natom; i++)
       nrs[i] = 0;
       
    for (i=1; i <= natom; i++)
    {
        for (j=1; j <= natom; j++)
        {
            if (fabs(atom.atomwt[i]-atom.atomwt[j]) <= TOLER)
            {
                if (nrs[i] == 0)
                {
                    dsum = 0;
                    for (k=0; k < 3; k++)
                    {
                        x1 = a[i][k];
                        x2 = b[j][k];
                        dsum += (a[i][k]-b[j][k])*(a[i][k]-b[j][k]);
                     }
                     if (dsum <= TOLER)
                     {
                         nrs[i] = j;
                         goto L_5;
                     }
                }
            }
        }
        free_ivector(nrs,0,natom+1);
        return FALSE;
L_5:
        continue;
    }
    free_ivector(nrs,0,natom+1);
    return TRUE;
}
void tcopy(double **a, double **b, int num)
{
    int i,j;
    for (i=1; i <= num; i++)
    {
        for (j=0; j < 3; j++)
          b[i][j] = a[i][j];
    }
}
void tform(double t[][3], double **a, double **b, int num)
{
    int i;
    for (i=1; i <= num; i++)
    {
        b[i][0] = t[0][0]*a[i][0] + t[0][1]*a[i][1] + t[0][2]*a[i][2];
        b[i][1] = t[1][0]*a[i][0] + t[1][1]*a[i][1] + t[1][2]*a[i][2];
        b[i][2] = t[2][0]*a[i][0] + t[2][1]*a[i][1] + t[2][2]*a[i][2];
    }
}
/* ====================================  */
void circset(double **a, double **aset,int ixyz,int *nset,int *numset)
{
    int i,j, i1, i2, i3, iset,j1;
    double q2, q3, an,ap,ad;

    i1 = ixyz;
    i2 = ixyz%3+1;
    i3 = i2%3 +1;
    i1--;
    i2--;
    i3--;
    for (i=1; i <= natom; i++)
    {
        aset[i][0] = atom.atomwt[i];
        aset[i][1] = a[i][i1];
        q2 = a[i][i2];
        q3 = a[i][i3];
        aset[i][2] = sqrt(q2*q2+q3*q3);
    }
// set 0  on axis atoms
    for (i=1; i <= natom; i++)
    {
        if (fabs(aset[i][2]) <= TOLER)
        {
            nset[i] = 0;
            aset[i][0] = 0.0;
        }
    }
// remaining sets
    iset = 0;
    for (i=1; i < natom; i++)  // loop is natom-1
    {
        if (aset[i][0] != 0.0)
        {
            iset++;
            nset[i] = iset;
            an = aset[i][0];
            ap = aset[i][1];
            ad = aset[i][2];
            aset[i][0] = 0.0;
            j1 = i+1;
            for (j=j1; j <= natom; j++)
            {
                if (fabs(aset[j][0]-an) <= TOL2 &&
                    fabs(aset[j][1]-ap) <= TOLER &&
                    fabs(aset[j][2]-ad) <= TOLER)
                    {
                       nset[j] = iset;
                       aset[j][0] = 0.0;
                    }
            }
        }
    }
    *numset = iset;
    for (i=1; i <= natom; i++)
       aset[i][0] = atom.atomwt[i];
}
/* =============================================== */            
void rotates(double **a,double **b,int iaxis,double theta)
{
    int i,j,i1,i2,i3;
    double s,c;
    double t[3][3];

    i1 = iaxis;
    i2 = iaxis%3 +1;
    i3 = i2%3 + 1;
    i1--;
    i2--;
    i3--;
    s = sin(theta);
    c = cos(theta);
    for (i=0; i < 3; i++)
    {
        for (j=0; j < 3; j++)
            t[i][j] = 0.0;
    }
    t[i1][i1] = 1.0;
    t[i2][i2] = c;
    t[i2][i3] = s;
    t[i3][i2] = -s;
    t[i3][i3] = c;
    tform(t,a,b,natom); 
}
/* =============================================== */            
void orients(double **a, double **b, double *p1, double *p2, int nset, int num)
{
    int i,j, i1, i2, i3;
    double trvec[3], t[3][3];
    double v1, v2, v3, vnorm;
    double alph, beta, gamm, v2v2, v3v3, v2233;
    double tol = 1.0e-17;

    trvec[0] = p1[0];
    trvec[1] = p1[1];
    trvec[2] = p1[2];

    for (i=1; i <= num; i++)
    {
        for (j=0; j < 3; j++)
           a[i][j] -= trvec[j];
    }

    i1 = nset;
    i2 = i1%3 + 1;
    i3 = i2%3 + 1;
    i1--;
    i2--;
    i3--;
    v1 = p2[i1] - p1[i1];
    v2 = p2[i2] - p1[i2];
    v3 = p2[i3] - p1[i3];
    vnorm = sqrt(v1*v1 + v2*v2 + v3*v3);
    if (vnorm == 0.0)
        return;
    if ( fabs(v1) < TOLER && fabs(v2) < TOLER && fabs(v3) < TOLER)
        return;
    if (fabs(v2) < tol && fabs(v3) < tol)
    {
        if (v1 > 0.0)
        {
            t[i1][i1] = 1.0;
            t[i2][i2] = 1.0;
            t[i3][i3] = 1.0;
            t[i1][i2] = t[i2][i1] = 0.0;
            t[i1][i3] = t[i3][i1] = 0.0;
            t[i3][i2] = t[i2][i3] = 0.0;
        } else
        {
            t[i1][i1] = -1.0;
            t[i2][i2] = 1.0;
            t[i3][i3] = 1.0;
            t[i1][i2] = t[i2][i1] = 0.0;
            t[i1][i3] = t[i3][i1] = 0.0;
            t[i3][i2] = t[i2][i3] = 0.0;
        }
    } else
    {
        alph = v1/vnorm;
        beta = v2/vnorm;
        gamm = v3/vnorm;
        v2v2 = v2*v2;
        v3v3 = v3*v3;
        v2233 = 1.0/(v2v2+v3v3);
        t[i1][i1] = alph;
        t[i1][i2] = beta;
        t[i1][i3] = gamm;
        t[i2][i1] = -t[i1][i2];
        t[i3][i1] = -t[i1][i3];
        t[i2][i3] = v2*v3*(alph-1.0)*v2233;
        t[i3][i2] = t[i2][i3];        
        t[i2][i2] = (v2v2*alph +v3v3)*v2233;
        t[i3][i3] = (v3v3*alph +v2v2)*v2233;
    }
    tform(t,a,b,num);
    for (i=1; i <= num; i++)
    {
        for (j=0; j < 3; j++)
           b[i][j] += trvec[j];
    }        
}
/* ==================================  */
int findi(double **a, double **b)
{
    int iz;
    sym_invert(a,b);
    iz = tcomp(a,b);
    return (iz);
}
int findh(double **a, double **b, int nset)
{
    int iz;
    reflect(a,b,nset);
    iz = tcomp(a,b);
    return (iz);
}
int findcn(double **a, double **b,double *p1, double *p2, int *maxcn)
{
    int i,j,k,l, iz, iaxis;
    int i1, i2;
    double halfpi, pi,theta;
    double t[3][3];
    
    halfpi = 2.0*atan(1.0);
    pi = 2.0*halfpi;
    iaxis = 0;
    for (i=6; i >= 2; i--)
    {
        theta = 2.0*pi/(float)i;
        for (j=0; j < 3; j++)
        {
            for (k=0; k < 3; k++)
            {
                for (l=0; l < 3; l++)
                {
                    t[k][l] = 0.0;
                }
            }
            i1 = j+1;
            i2 = j+2;
            if (i1 > 2) i1 -= 3;
            if (i2 > 2) i2 -= 3;
            t[i1][i1] = cos(theta);
            t[i2][i2] = cos(theta);
            t[i1][i2] = sin(theta);
            t[i2][i1] = -sin(theta);
            t[j][j] = 1.0;
            tform(t,a,b,natom);
            iz = tcomp(a,b);
            if  (iz == TRUE)
            {
                *maxcn = i;
                for (k=0; k < 3; k++)
                {
                    p1[k] = 0.0;
                    p2[k] = 0.0;
                    if (k == j) p2[k] = 1.0;
                }
                orients(b,a,p1,p2,3,natom);
                return TRUE;
            }
        }
    }
    *maxcn = 1;
    return FALSE;
}
int findsn(double **a, double **b, double **c, int *maxcn)
{
    int iz;
    double pi,theta;
    pi = 4.0*atan(1.0);
    theta = pi/(*maxcn);
    rotates(a,b,3,theta);
    reflect(b,c,3);
    iz = tcomp(a,c);
    return iz;
}
int findc2(double **a, double **b, double **aset, int *nset)
{
    int i,j,k, j1, numset, iset, jset, iz;
    double halfpi, pi, x,y,theta, xi,yi;
    double proi, ani,disi;
    
    halfpi = 2.0*atan(1.0);
    pi = 2.0*halfpi;

    circset(a,aset,3,nset,&numset);
    for (i=1; i <= numset; i++)
    {
        for (j=1; j <= natom; j++)
        {
            if(nset[j] == i)  // atom belongs to set
            {
                if ( fabs(aset[j][1]) < TOLER )
                {
                    j1 = j+1;
                    for (k=j1; k <= natom; k++)
                    {
                        if (nset[k] == i)
                        {
                            x = (a[j][0]+a[k][0])*0.5;
                            y = (a[j][1]+a[k][1])*0.5;
                            theta = halfpi;
                            if (fabs(y) > TOLER) theta = -atan(x/y);
                            rotates(a,b,3,theta);
                            rotates(b,aset,2,pi);
                            iz = tcomp(b,aset);
                            if (iz == TRUE)
                            {
                                tcopy(b,a,natom);
                                return TRUE;
                            }
                        }
                    }
                }
            }
        }
    }
//
    circset(a,aset,3,nset,&numset);
    iset = 1;
    for (i=1; i <= natom; i++)
      if (nset[i] == iset) goto L_10;
    return FALSE;
//
L_10:
    proi = aset[i][1];
    ani = aset[i][0];
    disi = aset[i][2];
    xi = a[i][0];
    yi = a[i][1];
    j1 = iset+1;
    for (j=j1; j <= numset; j++)
    {
        for (k=1; k <= natom; k++)
        {
            if (nset[k] == j)
            {
                jset = j;
                goto L_20;
            }
        }
         return FALSE;
L_20:
       if ( fabs(proi+aset[k][1]) > TOLER || fabs(ani-aset[k][0]) > TOL2 ||
            fabs(disi-aset[k][2]) > TOLER )
          goto L_30;

       for (k=1; k <= natom; k++)
       {
         if (nset[k] == jset)
         {
            x = (xi+a[k][0])*0.5;
            y = (yi+a[k][1])*0.5;
            theta = halfpi;
            if (fabs(y) > TOL2) theta = -atan(x/y);
            rotates(a,b,3,theta);
            rotates(b,aset,2,pi);
            iz = tcomp(b,aset);
            if (iz == TRUE)
            {
                tcopy(b,a,natom);
                return TRUE;
            }
        }
       }
L_30:
       continue;
    }
    return FALSE;
}
int findv(double **a, double **b, double **aset, int *nset)
{
    int i,j, j1, numset, iset, jset, iz;
    double halfpi, pi, x,y,theta;
    
    halfpi = 2.0*atan(1.0);
    pi = 2.0*halfpi;

    circset(a,aset,3,nset,&numset);
    iset = 1;
    for (i=1; i < natom; i++)
    {
        if (nset[i] == iset)
        {
           jset = i;
           goto L_5;
        }
    }
    return FALSE;
L_5:
    j1 = jset + 1;
    for (j=j1; j <= natom; j++)
    {
        if (nset[j] == iset)
        {
            x = (a[i][0]+a[j][0])*0.5;
            y = (a[i][1]+a[j][1])*0.5;
            if (fabs(x) <= TOLER && fabs(y) <= TOLER)
            {
                x = 0.5*a[j][0];
                y = 0.5*a[j][1];
            }
            theta = halfpi;
            if (fabs(y) > TOLER) theta = -atan(x/y);
            rotates(a,b,3,theta);
            reflect(b,aset,1);
            iz = tcomp(b,aset);
            if (iz == TRUE)
            {
                tcopy(b,a,natom);
                return TRUE;
            }
        }
    }
    return FALSE;
}
/* ======================================================= */
// compute moments of inertia
void secmom(double **xyz, double *xi,double *yi, double *zi, int *ishape,double t[3][3])
{
    int i,j,k,ii;
    float b[3];
    double a[6], eigval[3], eigvec[3][3], atemp;
    double summ, wti, xx,yy,zz;
    
    for (i=0; i < 3; i++)
    {
        b[i] = 0.0;
        eigval[i] = 0.0;
        for (j=0; j < 3; j++)
          eigvec[i][j] = 0.0;
    }
    for (i=0; i < 6; i++)
       a[i] = 0.0;

    summ = 0.0;
    for (i=1; i <= natom; i++)
    {
        wti = atom.atomwt[i];
        for (j=0; j < 3; j++)
            b[j] += xyz[i][j]*wti;
        summ += wti;
    }
    for (i=0; i < 3; i++)
       b[i] /= summ;
    for (i=1; i <= natom; i++)
    {
        for (j=0; j < 3; j++)
           xyz[i][j] -= b[j];
    }
   // compute moments of inertia , sort and place smallest moment on x axis
   for (i=1; i <= natom; i++)
   {
       wti = atom.atomwt[i];
       *xi = xyz[i][0];
       *yi = xyz[i][1];
       *zi = xyz[i][2];
       xx = *xi* *xi;      
       yy = *yi* *yi;      
       zz = *zi* *zi;
       a[0] += wti*(yy+zz);
       a[2] += wti*(zz+xx);
       a[5] += wti*(xx+yy);
       a[1] -= wti**xi**yi; 
       a[3] -= wti**xi**zi; 
       a[4] -= wti**yi**zi;
   }
   summ = fabs(a[1]) + fabs(a[3]) + fabs(a[4]);
   if (summ <= 1.0e-15)
   {
       eigval[0] = a[0];
       eigval[1] = a[2];
       eigval[2] = a[5];
       eigvec[0][0] = 1.0;
       eigvec[1][1] = 1.0;
       eigvec[2][2] = 1.0;
   } else
   {
       hqrii(a,3,3,eigval,eigvec);
   }

   for (i=0; i < 2; i++)
   {
       ii = i + 1;
       for (j=ii; j < 3; j++)
       {
           if (eigval[i] > eigval[j])
           {
               atemp = eigval[i];
               eigval[i] = eigval[j];
               eigval[j] = atemp;
               for (k=0; k < 3; k++)
               {
                   atemp = eigvec[k][i];
                   eigvec[k][i] = eigvec[k][j];
                   eigvec[k][j] = atemp;
               }
           }
       }
   }        
   *xi = eigval[0] / 6.0228;
   *yi = eigval[1] / 6.0228;
   *zi = eigval[2] / 6.0228;
// now rotate the molecule
    for (i=1; i <= natom; i++)
    {
        for (j=0; j < 3; j++)
        {
            b[j] = 0.0;
            for (k=0; k < 3; k++)
                b[j] += eigvec[j][k]*xyz[i][k];
        }
        xyz[i][0] = b[0];
        xyz[i][1] = b[1];
        xyz[i][2] = b[2];
    }
// assign shape
    *ishape = 0;
    if (*xi <= 0.01)
      *ishape = 1;
    else
    {
        if (fabs(*xi+*yi-*zi) <= 0.01)
           *ishape = 6;
        else
        {
            if (fabs( (*zi/(*yi))-(*yi/(*xi))) <= 0.01)
               *ishape = 5;
            else
            {
                atemp = (2.0/(*yi) - 1.0/(*zi) - 1.0/(*xi))/ (1.0/(*xi) - 1.0/(*zi));
                if (atemp <= -0.999)
                   *ishape = 3;
                else
                {
                    if (atemp >= 0.999)
                        *ishape = 4;
                    else
                        *ishape = 2;
                }
            }
        }
    }
}
// assign point group, symmetry number and chirality
void vibsym(int mode)
{
    int i;
    double **a;
    char ngrp[4];
    
    a = dmatrix (0, natom+1, 0,3);
    for (i=1; i <= natom; i++)
    {
        a[i][0] = atom.x[i];
        a[i][1] = atom.y[i];
        a[i][2] = atom.z[i];
    }

    if (mode == 0)
    {
        strcpy(ngrp,"");
        ptgrp(mode, ngrp,a);
    }
    free_dmatrix(a ,0, natom+1, 0,3);
}
    
void ptgrp(int mode,char *ngrp,double **a)
{
    int i, j, k, ishape, numsym, iaxis, maxcn, xxx;
    int iz, iz1, iz2, iopt;
    int chiral, *nset;
    double **b, **c;
    double t[3][3], p1[3], p2[3];
    double xi, yi, zi;
    float halfpi, pi;
    char string[256];

    numsym = 0;
    halfpi = 2.0*atan(1.0);
    pi = 2.0*halfpi;
    chiral = FALSE;
    strcpy(ngrp,"Cs");
    
       
    b = dmatrix (0, natom+1, 0,3);
    c = dmatrix (0, natom+1, 0,3);
    nset = ivector(0,natom+1);

    secmom(a, &xi, &yi, &zi, &ishape, t);
    
    if (ishape == 1)   // linear
    {
      p1[0] = 0.0;
      p1[1] = 0.0;
      p1[2] = 0.0;
      p2[0] = 1.0;
      p2[1] = 0.0;
      p2[1] = 0.0;
      orients(a,b,p1,p2,3,natom);
      tcopy(b,a,natom);
      iaxis = findh(a,b,3);
      if (iaxis == TRUE)
      {
          strcpy(ngrp,"D*h");
          numsym = 2;
      } else
      {
          strcpy(ngrp,"C*v");
          numsym = 1;
      }
    } else if (ishape == 2) // asymmetric
     /*------------------------*
      ASyMMETRiC TOP MOLECULES
     *------------------------*

     These molecules can have no axes of order greater than 2.  Thus
     the possible point groups are:  D2H, D2, C2V, C2H, C2, Ci, CS, 
     and C1.  FiNDCN routine finds the principal axes and align the
     principal axes with the z-axis. */
    {
        iaxis = findcn(a,b,p1,p2,&maxcn);
        if (maxcn == 1)
        {
            j = findh(a,b,3);
            if (j == TRUE)
            {
                strcpy(ngrp,"Cs");
                numsym = 1;
            }else
            {
                i = findh(a,b,2);
                if (i == TRUE)
                {
                    strcpy(ngrp,"Cs");
                    numsym = 1;
                } else
                {
                   k = findh(a,b,1);
                   if (k == TRUE)
                   {
                      strcpy(ngrp,"Cs");
                      numsym = 1;
                   } else
                   {
                       xxx = findi(a,b);
                       if (xxx == TRUE)
                       {
                           strcpy(ngrp,"Ci");
                           numsym = 1;
                       } else
                       {
                           strcpy(ngrp,"C1");
                           numsym = 1;
                           chiral = TRUE;
                       }
                   }
                }
            }
        } else  // maxcn == 2
        {
            iaxis = findc2(a,b,c,nset);
            if (iaxis == FALSE)
            {
                i = findv(a,b,c,nset);
                if (i == TRUE)
                {
                    strcpy(ngrp,"C2v");
                    numsym = 2;
                } else
                {
                    j = findh(a,b,3);
                    if (j == TRUE)
                    {
                        strcpy(ngrp,"C2h");
                        numsym = 2;
                    } else
                    {
                        strcpy(ngrp,"C2");
                        numsym = 2;
                        chiral = TRUE;
                    }
                }
            } else
            {
                i = findh(a,b,3);
                if (i == TRUE)
                {
                    strcpy(ngrp,"D2h");
                    numsym = 4;
                } else
                {
                    j = findv(a,b,c,nset);
                    if (j == TRUE)
                    {
                        strcpy(ngrp,"D2d");
                        numsym = 4;
                    } else
                    {
                        strcpy(ngrp,"D2");
                        numsym = 4;
                        chiral = TRUE;
                    }
                }
            }
        }
    } else if (ishape == 3 || ishape == 4) //prolate
    {
     /*-----------------------*
      SyMMETRiC TOP MOLECULES
     *-----------------------*

     These molecules can belong to any axial point group, thus only
     the cubic point groups (T, Td, Th, O, Oh, i, ih) are impossible.
     However, except in rare cases the unique axis is a rotation axis
     of order 3 or greater. */

        iz = findcn(a,b,p1,p2,&maxcn);
        iz = findsn(a,b,c,&maxcn);
        if (iz == TRUE)
        {
            iz1 = findv(a,b,c,nset);
            if (iz1 == FALSE)
            {
               sprintf(ngrp,"S%1d",2*maxcn);
               numsym = maxcn;
               goto L_2000;
            }
        }
        iz = findc2(a,b,c,nset);
        if (iz == FALSE)
        {
           iz1 = findv(a,b,c,nset);
           if (iz1 == TRUE)
           {
                 sprintf(ngrp,"C%1dv",maxcn);
                 numsym = maxcn;
           } else
           {
                    iz2 = findh(a,b,3);
                    if (iz2 == TRUE)
                    {
                        sprintf(ngrp,"C%1dh",maxcn);
                        numsym = maxcn;
                    } else
                    {
                        sprintf(ngrp,"C%1d",maxcn);
                        numsym = maxcn;
                        chiral = TRUE;
                    }
           }
        }else
        {
                iz1 = findh(a,b,3);
                if (iz1 == TRUE)
                {
                    sprintf(ngrp,"D%1dh",maxcn);
                    numsym = 2*maxcn;
                } else
                {
                    iz2 = findv(a,b,c,nset);
                    if (iz2 == TRUE)
                    {
                        sprintf(ngrp,"D%1dd",maxcn);
                        numsym = 2*maxcn;
                    } else
                    {
                        sprintf(ngrp,"D%1d",maxcn);
                        numsym = 2*maxcn;
                        chiral = TRUE;
                    }
                }
        }
        // SPHERiCAL TOP MOLECULES       
    } else if (ishape == 5) // spherical
    {
        sphere(a,b,c,nset,&maxcn);
        iopt = maxcn - 2;
        if (iopt == 1)
        {
            iz = findi(a,b);
            if (iz == TRUE)
            {
                sprintf(ngrp,"Th");
                numsym = 12;
            } else
            {
                iz1 = findv(a,b,c,nset);
                if (iz1 == TRUE)
                {
                    sprintf(ngrp,"Td");
                    numsym = 12;
                } else
                {
                    sprintf(ngrp,"T");
                    numsym = 12;
                    chiral = TRUE;
                }
            }
        } else if (iopt == 2)
        {
            iz = findi(a,b);
            if (iz == TRUE)
            {
                sprintf(ngrp,"Oh");
                numsym = 24;
            } else
            {
                sprintf(ngrp,"O");
                numsym = 24;
                chiral = TRUE;
            }
        } else if (iopt == 3)
        {
            iz = findi(a,b);
            if (iz == TRUE)
            {
                sprintf(ngrp,"Ih");
                numsym = 60;
            } else
            {
                sprintf(ngrp,"I");
                numsym = 60;
                chiral = TRUE;
            }
        }
    } else if (ishape == 6) // planar
    {
        iz = findh(a,b,1);
        iz1 = findh(a,b,2);
        iz2 = findcn(a,b,p1,p2,&maxcn);
        if (maxcn == 2)
        {
            if (iz == TRUE && iz1 == TRUE)
            {
                sprintf(ngrp,"D%1dh",maxcn);
                numsym = maxcn*2;
            } else
            {
                if (iz == TRUE || iz1 == TRUE)
                {
                    sprintf(ngrp,"C%1dv",maxcn);
                    numsym = maxcn;
                }else
                {
                    i = findh(a,b,3);
                    if (i == TRUE)
                    {
                        sprintf(ngrp,"C%1dh",maxcn);
                        numsym = maxcn;
                    } else
                    {
                        sprintf(ngrp,"C%1d",maxcn);
                        numsym = maxcn;
                        chiral = TRUE;
                    }
                }
            }
        } else
        {
            iz = findsn(a,b,c,&maxcn);
            if (iz == TRUE)
            {
                iz1 = findv(a,b,c,nset);
                if (iz1 == FALSE)
                {
                   sprintf(ngrp,"S%1d",2*maxcn);
                   numsym = maxcn;
                   goto L_2000;
                }
            }
            iz = findc2(a,b,c,nset);
            if (iz == FALSE)
            {
                iz1 = findv(a,b,c,nset);
                if (iz1 == TRUE)
                {
                    sprintf(ngrp,"C%1dv",maxcn);
                    numsym = maxcn;
                } else
                {
                    iz2 = findh(a,b,3);
                    if (iz2 == TRUE)
                    {
                        sprintf(ngrp,"C%1dh",maxcn);
                        numsym = maxcn;
                    } else
                    {
                        sprintf(ngrp,"C%1d",maxcn);
                        numsym = maxcn;
                        chiral = TRUE;
                    }
                }
            } else
            {
                iz1 = findh(a,b,3);
                if (iz1 == TRUE)
                {
                    sprintf(ngrp,"D%1dh",maxcn);
                    numsym = 2*maxcn;
                } else
                {
                    iz2 = findv(a,b,c,nset);
                    if (iz2 == TRUE)
                    {
                        sprintf(ngrp,"D%1dd",maxcn);
                        numsym = 2*maxcn;
                    } else
                    {
                        sprintf(ngrp,"D%1d",maxcn);
                        numsym = 2*maxcn;
                        chiral = TRUE;
                    }
                }
            }       
        }         
    }
L_2000:   
    free_dmatrix(b ,0, natom+1, 0,3);
    free_dmatrix(c ,0, natom+1, 0,3);
    free_ivector(nset, 0, natom+1);

    if (mode == 2)
    {
       symmetry.chiral = chiral;
       symmetry.numsym = numsym;
       symmetry.ishape = ishape;
       symmetry.xi = xi;
       symmetry.yi = yi;
       symmetry.zi = zi;
    }
    
    if (mode == 0)
    {
       sprintf(string,"Symmetry group: %s\nMoments of Inertia\nIx: %7.2f\nIy: %7.2f\nIz: %7.2f",ngrp,xi,yi,zi);
       message_alert(string,"Symmetry Information");
    }
}
/* ========================================== */
void hqrii(double *a,int n,int m,double *e,double v[3][3])
{
#define V(I_,J_)        (*(v+(I_)*(n)+(J_)))
        int _do0, _do1, _do10, _do11, _do12, _do13, _do14, _do15, 
         _do16, _do17, _do18, _do19, _do2, _do20, _do21, _do22, _do3, 
         _do4, _do5, _do6, _do7, _do8, _do9, i, i_, ig, ii, im1, iord, 
         ip1, irank, itere, itere_, j, j_, jrank, k, k_, kk, kk_, kp1, 
         kpiv, krank, l, ll, nm1, nm2=0;
        double c, del, ee, eps, eps1, eps2, eps3, ff, fn, gersch, h, r, 
         ra, rn, s, seps, sinv, sorter, sum, summ, t, u, vn, w[5][5], 
         ww, z, zero;

        double *const A = &a[0] - 1;
        double *const E = &e[0] - 1;

        /*************************************************************
         *
         * HQRII IS A DIAGONALISATION ROUTINE, WRITTEN BY YOSHITAKA BEPPU OF
         *       NAGOYA UNIVERSITY, JAPAN.
         *       FOR DETAILS SEE 'COMPUTERS & CHEMISTRY' VOL.6 1982. PAGE 000.
         *
         * ON INPUT    A       = MATRIX TO BE DIAGONALISED (PACKED CANONICAL)
         *             N       = SIZE OF MATRIX TO BE DIAGONALISED.
         *             M       = NUMBER OF EIGENVECTORS NEEDED.
         *             E       = ARRAY OF SIZE AT LEAST N
         *             V       = ARRAY OF SIZE AT LEAST NMAX*M
         *
         * ON OUTPUT   E       = EIGENVALUES
         *             V       = EIGENVECTORS IN ARRAY OF SIZE NMAX*M
         *
         ************************************************************************ */

        /* EPS3 AND EPS ARE MACHINE-PRECISION DEPENDENT
         * */
        eps3 = 1.e-30;
        zero = 0.e0;
        ll = (n*(n + 1))/2 + 1;
        eps = 1.e-8;
        iord = -1;
        nm1 = n - 1;
        if( n == 2 )
                goto L_90;
        nm2 = n - 2;
        krank = 0;
        /*     HOUSEHOLDER TRANSFORMATION */
        for( k = 1, _do0 = nm2; k <= _do0; k++ )
        {
                k_ = k - 1;
                kp1 = k + 1;
                krank += k;
                w[k_][1] = A[krank];
                sum = 0.;
                jrank = krank;
                for( j = kp1, _do1 = n; j <= _do1; j++ )
                {
                        j_ = j - 1;
                        w[j_][1] = A[jrank + k];
                        jrank += j;
                        sum += (w[j_][1])*(w[j_][1]);
                }
                s = SIGN( sqrt( sum ), w[kp1 - 1][1] );
                w[k_][0] = -s;
                w[kp1 - 1][1] += s;
                A[k + krank] = w[kp1 - 1][1];
                h = w[kp1 - 1][1]*s;
                if( fabs( h ) < eps3 )
                        goto L_80;
                summ = 0.e0;
                irank = krank;
                for( i = kp1, _do2 = n; i <= _do2; i++ )
                {
                        i_ = i - 1;
                        sum = 0.e0;
                        for( j = kp1, _do3 = i; j <= _do3; j++ )
                        {
                                j_ = j - 1;
                                sum += A[j + irank]*w[j_][1];
                        }
                        if( i >= n )
                                goto L_40;
                        ip1 = i + 1;
                        jrank = i*(i + 3)/2;
                        for( j = ip1, _do4 = n; j <= _do4; j++ )
                        {
                                j_ = j - 1;
                                sum += A[jrank]*w[j_][1];
                                jrank += j;
                        }
L_40:
                        w[i_][0] = sum/h;
                        irank += i;
                        summ += w[i_][0]*w[i_][1];
                }
                u = summ*0.5e0/h;
                jrank = krank;
                for( j = kp1, _do5 = n; j <= _do5; j++ )
                {
                        j_ = j - 1;
                        w[j_][0] = w[j_][1]*u - w[j_][0];
                        for( i = kp1, _do6 = j; i <= _do6; i++ )
                        {
                                i_ = i - 1;
                                A[i + jrank] += w[i_][0]*w[j_][1] + w[j_][0]*w[i_][1];
                        }
                        jrank += j;
                }
L_80:
                A[krank] = h;
        }
L_90:
        w[nm1 - 1][1] = A[(nm1*(nm1 + 1))/2];
        w[n - 1][1] = A[(n*(n + 1))/2];
        w[nm1 - 1][0] = A[nm1 + (n*(n - 1))/2];
        w[n - 1][0] = 0.e0;
        gersch = fabs( w[0][1] ) + fabs( w[0][0] );
        for( i = 1, _do7 = nm1; i <= _do7; i++ )
        {
                i_ = i - 1;
                if ( fabs( w[i_ + 1][1] ) + fabs( w[i_][0] ) + fabs( w[i_ + 1][0] ) > gersch)
                   gersch = fabs( w[i_ + 1][1] ) + fabs( w[i_][0] ) + fabs( w[i_ + 1][0] );
        }
        del = eps*gersch;
        for( i = 1, _do8 = n; i <= _do8; i++ )
        {
                i_ = i - 1;
                w[i_][2] = w[i_][0];
                E[i] = w[i_][1];
                v[m-1][i_] = E[i];
        }
        if( fabs( del ) < eps3 )
                goto L_220;
        /*     QR-METHOD WITH ORIGIN SHIFT */
        k = n;
L_120:
        l = k;
L_130:
        if( fabs( w[l - 2][2] ) < del )
                goto L_140;
        l -= 1;
        if( l > 1 )
                goto L_130;
L_140:
        if( l == k )
                goto L_170;
        ww = (E[k - 1] + E[k])*0.5e0;
        r = E[k] - ww;
        z = SIGN( sqrt( (w[k - 2][2])*(w[k - 2][2]) + r*r ), r ) + ww;
        ee = E[l] - z;
        E[l] = ee;
        ff = w[l - 1][2];
        r = sqrt( ee*ee + ff*ff );
        j = l;
        goto L_160;
L_150:
        r = sqrt( (E[j])*(E[j]) + (w[j - 1][2])*(w[j - 1][2]) );
        w[j - 2][2] = s*r;
        ee = E[j]*c;
        ff = w[j - 1][2]*c;
L_160:
        c = E[j]/r;
        s = w[j - 1][2]/r;
        ww = E[j + 1] - z;
        E[j] = (ff*c + ww*s)*s + ee + z;
        E[j + 1] = c*ww - s*ff;
        j += 1;
        if( j < k )
                goto L_150;
        w[k - 2][2] = E[k]*s;
        E[k] = E[k]*c + z;
        goto L_120;
L_170:
        k -= 1;
        if( k > 1 )
                goto L_120;
        /*    *    *    *    *    *    *    *    *    *    *    *    *
         *
         *   AT THIS POINT THE ARRAY 'E' CONTAINS THE UN-ORDERED EIGENVALUES
         *
         *    *    *    *    *    *    *    *    *    *    *    *    *
         *     STRAIGHT SELECTION SORT OF EIGENVALUES */
        sorter = 1.e0;
        if( iord < 0 )
                sorter = -1.e0;
        j = n;
L_180:
        l = 1;
        ii = 1;
        ll = 1;
        for( i = 2, _do9 = j; i <= _do9; i++ )
        {
                i_ = i - 1;
                if( (E[i] - E[l])*sorter > 0.e0 )
                        goto L_190;
                l = i;
                goto L_200;
L_190:
                ii = i;
                ll = l;
L_200:
                ;
        }
        if( ii == ll )
                goto L_210;
        ww = E[ll];
        E[ll] = E[ii];
        E[ii] = ww;
L_210:
        j = ii - 1;
        if( j >= 2 )
                goto L_180;
L_220:
        if( !m )
                return;
        /***************
         *  ORDERING OF EIGENVALUES COMPLETE.
         ***************
         *      INVERSE-ITERATION FOR EIGENVECTORS */
        fn = (float)( n );
        eps1 = 1.e-5;
        seps = sqrt( eps );
        eps2 = 0.05e0;
        rn = 0.e0;
        ra = eps*0.6180339887485e0;
        /*    0.618... IS THE FIBONACCI NUMBER (-1+SQRT(5))/2. */
        ig = 1;
        for( i = 1, _do10 = m; i <= _do10; i++ )
        {
                i_ = i - 1;
                im1 = i - 1;
                for( j = 1, _do11 = n; j <= _do11; j++ )
                {
                        j_ = j - 1;
                        w[j_][2] = 0.e0;
                        w[j_][3] = w[j_][0];
                        w[j_][4] = v[m-1][j_] - E[i];
                        rn += ra;
                        if( rn >= eps )
                                rn -= eps;
                        v[i_][j_] = rn;
                }
                for( j = 1, _do12 = nm1; j <= _do12; j++ )
                {
                        j_ = j - 1;
                        if( fabs( w[j_][4] ) >= fabs( w[j_][0] ) )
                                goto L_240;
                        w[j_][1] = -w[j_][4]/w[j_][0];
                        w[j_][4] = w[j_][0];
                        t = w[j_ + 1][4];
                        w[j_ + 1][4] = w[j_][3];
                        w[j_][3] = t;
                        w[j_][2] = w[j_ + 1][3];
                        if( fabs( w[j_][2] ) < eps3 )
                                w[j_][2] = del;
                        w[j_ + 1][3] = 0.e0;
                        goto L_250;
L_240:
                        if( fabs( w[j_][4] ) < eps3 )
                                w[j_][4] = del;
                        w[j_][1] = -w[j_][0]/w[j_][4];
L_250:
                        w[j_ + 1][3] += w[j_][2]*w[j_][1];
                        w[j_ + 1][4] += w[j_][3]*w[j_][1];
                }
                if( fabs( w[n - 1][4] ) < eps3 )
                        w[n - 1][4] = del;
                for( itere = 1; itere <= 5; itere++ )
                {
                        itere_ = itere - 1;
                        if( itere == 1 )
                                goto L_280;
                        for( j = 1, _do13 = nm1; j <= _do13; j++ )
                        {
                                j_ = j - 1;
                                if( fabs( w[j_][2] ) < eps3 )
                                        goto L_270;
                                t = v[i_][j_];
                                v[i_][j_] = v[i_][j_+1];
                                v[i_][j_+1] = t;
L_270:
                                v[i_][j_+1] += v[i_][j_]*w[j_][1];
                        }
L_280:
                        v[i_][n-1] /= w[n - 1][4];
                        v[i_][nm1 - 1] = (v[i_][nm1 - 1] - v[i_][n-1]*w[nm1 - 1][3])/w[nm1 - 1][4];
                        vn = 1.0e-20;
                        if (fabs( v[i_][n-1] ) > vn)
                           vn = fabs( v[i_][n-1] );
                        if ( fabs( v[i_][nm1 - 1] ) > vn)
                           vn =  fabs( v[i_][nm1 - 1] );
                        if( n == 2 )
                                goto L_300;
                        k = nm2;
L_290:
                        v[i_][k - 1] = (v[i_][k - 1] - v[i_][k]*w[k - 1][3] - v[i_][k + 1]*
                         w[k - 1][2])/w[k - 1][4];
                        if (vn < 1.0e-20)
                            vn = 1.0e-20;
                        if (fabs( v[i_][k - 1] ) > vn)
                           vn = fabs( v[i_][n-1] );
                        k -= 1;
                        if( k >= 1 )
                                goto L_290;
L_300:
                        s = eps1/vn;
                        for( j = 1, _do14 = n; j <= _do14; j++ )
                        {
                                j_ = j - 1;
                                v[i_][j_] *= s;
                        }
                        if( itere > 1 && vn > 1 )
                                goto L_330;
                }
                /*     TRANSFORMATION OF EIGENVECTORS */
L_330:
                if( n == 2 )
                        goto L_380;
                krank = nm2*(n + 1)/2;
                kpiv = nm2*nm1/2;
                for( k = nm2; k >= 1; k-- )
                {
                        k_ = k - 1;
                        kp1 = k + 1;
                        if( fabs( A[kpiv] ) <= eps3 )
                                goto L_360;
                        sum = 0.e0;
                        for( kk = kp1, _do15 = n; kk <= _do15; kk++ )
                        {
                                kk_ = kk - 1;
                                sum += A[krank]*v[i_][kk_];
                                krank += kk;
                        }
                        s = -sum/A[kpiv];
                        for( kk = n, _do16 = kp1; kk >= _do16; kk-- )
                        {
                                kk_ = kk - 1;
                                krank -= kk;
                                v[i_][kk_] += A[krank]*s;
                        }
L_360:
                        kpiv -= k;
                        krank -= kp1;
                }
L_380:
                for( j = ig, _do17 = i; j <= _do17; j++ )
                {
                        j_ = j - 1;
                        if( fabs( E[j] - E[i] ) < eps2 )
                                goto L_400;
                }
                j = i;
L_400:
                ig = j;
                if( ig == i )
                        goto L_430;
                /*     RE-ORTHOGONALISATION */
                for( k = ig, _do18 = im1; k <= _do18; k++ )
                {
                        k_ = k - 1;
                        sum = 0.e0;
                        for( j = 1, _do19 = n; j <= _do19; j++ )
                        {
                                j_ = j - 1;
                                sum += v[k_][j_]*v[i_][j_];
                        }
                        s = -sum;
                        for( j = 1, _do20 = n; j <= _do20; j++ )
                        {
                                j_ = j - 1;
                                v[i_][j_] += v[k_][j_]*s;
                        }
                }
                /*     NORMALISATION */
L_430:
                sum = 1.e-24;
                for( j = 1, _do21 = n; j <= _do21; j++ )
                {
                        j_ = j - 1;
                        sum += (v[i_][j_])*(v[i_][j_]);
                }
                sinv = 1.e0/sqrt( sum );
                for( j = 1, _do22 = n; j <= _do22; j++ )
                {
                        j_ = j - 1;
                        v[i_][j_] *= sinv;
                }
        }
        return;
#undef  V
}

