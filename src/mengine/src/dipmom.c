#define EXTERN extern

#include "pcwin.h"
#include "pcmod.h"
#include "dipmom.h"

double get_dipole_moment(void);
void charge_dipole(int natom,double *x,double *y,double *z,double *atomwt,double *charge);

// =====================
double get_dipole_moment()
{
  return dipolemom.total;
}
// =====================
void charge_dipole(int natom,double *x,double *y,double *z,double *atomwt,double *charge)
{
    int i;
    double weight, xcenter, ycenter, zcenter;
    double xsum, ysum, zsum, debye;

    dipolemom.xdipole = 0.0;
    dipolemom.ydipole = 0.0;
    dipolemom.zdipole = 0.0;
    debye = 4.8033324;

    weight = 0.0;
    xcenter = 0.0;
    ycenter = 0.0;
    zcenter = 0.0;
    for (i=1; i <= natom; i++)
    {
        weight += atom.atomwt[i];
        xcenter += atom.x[i]*atom.atomwt[i];
        ycenter += atom.y[i]*atom.atomwt[i];
        zcenter += atom.z[i]*atom.atomwt[i];
    }
    xcenter /= weight;
    ycenter /= weight;
    zcenter /= weight;

    xsum = ysum = zsum = 0.0;
    for(i=1; i <= natom; i++)
    {
        xsum += (atom.x[i]-xcenter)*atom.charge[i];
        ysum += (atom.y[i]-ycenter)*atom.charge[i];
        zsum += (atom.z[i]-zcenter)*atom.charge[i];
    }
    dipolemom.xdipole = xsum*debye;
    dipolemom.ydipole = ysum*debye;
    dipolemom.zdipole = zsum*debye;
    dipolemom.total = sqrt(dipolemom.xdipole*dipolemom.xdipole +
          dipolemom.ydipole*dipolemom.ydipole + dipolemom.zdipole*dipolemom.zdipole);
}
