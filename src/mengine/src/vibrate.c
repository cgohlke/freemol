#define EXTERN extern

#include "pcwin.h"

#include "asnsym.h"
#include "utility.h"
#include "bonds_ff.h"
#include "vibrate.h"

char *symname[] = {"A","A'","A''","Ag","Au",
                "A1","A2","A1'","A2'","A1''","A2''",
                "A1g","A1u","A2g","A2u","B",
                "B1","B2","B3","Bg","Bu","B1g",
                "B1u","B2g","B2u","B3g","B3u","E",
                "E1","E2","E3","E'","E''",
                "Eg","Eu","E1'","E1''","E2'",
                "E2''","E1g","E1u","E2g","E2u","T","Tg","Tu" };

              
EXTERN double hesscut;

// FILE * fopen_path ( char * , char * , char * ) ;

static void tred2(double **,int , double *, double *);
static void tqli(double *,double *,int ,double **);
static void assign_vibsym(int, char *,char *,int *,double **);
double get_total_energy(void);
void hessian(int, double *,int *, int *, int *,double *);
void vibrate(int natom,double *atomwt,double *x,double *y,double *z,double *charge);
double SIGN(double,double);
//void vibcharge_dipole(double **);
//void vibdipole_dipole(double **);

void vibrate(int natom,double *atomwt,double *x,double *y,double *z,double *charge)
{
    int i, k, j, ii, jj, nsym,imix,mm;
    int maxvar, maxhess, ihess, nfreq, maxvib;
    int *hinit, *hstop, *hindex;
    int *vibsym;
    double mass, efact,cnosym,cimix,temp,qptf,wdt,efwt,cpfact,efwt2;
    double xyzi,ssym,entrov1,entrov2,ezpe;
    double etrans,htrans,strans,gtrans,cptrans;
    double erot,hrot,srot,grot,cprot;
    double evib,hvib,svib,gvib,cpvib;
    double emix,hmix,smix,gmix,cpmix;
    double etot,htot,stot,gtot,cptot;
    double *hdiag, *h;
    double **matrix, *a;
    double **acoord;
    double *eigen, **vects, *mass2;
    double factor, convert, lightspd, vnorm;
    double *dipole;
    double dqx,dqy,dqz,dqxx,dqyy,dqzz,dqxy,dqxz,dqyx,dqyz,dqzx,dqzy;
    char ngrp[4], vsym[4];

    convert = 4.184e2;
    lightspd = 2.99792458e-2;
    factor = sqrt(convert)/(2.0*PI*lightspd);

    for (i=1; i <= natom; i++)
    {
        if (atomwt[i] <= 0.0)
           atomwt[i] = 0.001;
    }
          
    maxvar = 3*natom;
    maxhess = (3*natom*(3*natom-1))/2;
    maxvib = 3*natom;

    vibsym = ivector(0,maxvar);
    hinit = ivector(0,maxvar);
    hstop = ivector(0,maxvar);
    hindex = ivector(0,maxhess);
    hdiag = dvector(0,maxvar);
    h = dvector(0,maxhess);
    matrix = dmatrix(0, maxvib, 0, maxvib);
    a = dvector(0, maxvib);
    mass2 = dvector(0, natom+1);
    eigen = dvector(0, 3*natom);
    vects = dmatrix(0, 3*natom, 0, 3*natom);
    acoord = dmatrix (0, natom+1, 0,3);
    dipole = dvector(0,maxvar);

    for (i=1; i <= natom; i++)
        mass2[i] = sqrt(atomwt[i]);

     hesscut = 0.0;
     hessian(maxhess, h, hinit, hstop, hindex, hdiag);
     ihess = 0;
     for (i=0; i < maxvar; i++)
     {
        matrix[i][i] = hdiag[i];
        ii = i;
        for (k = hinit[i]; k < hstop[i]; k++)
        {
           ii++;
           matrix[i][ii] = h[k];
           matrix[ii][i] = h[k];
        }
     }
     nfreq = maxvar;

      tred2(matrix, nfreq, eigen, a);
      tqli(eigen,a,nfreq,matrix);
      
// now do mass weighted
    ihess = 0;
    for (i=0; i < maxvar; i++)
    {
        jj = i/3 + 1;
        matrix[i][i] = hdiag[i]/atomwt[jj];
        ii = i;
        for (k = hinit[i]; k < hstop[i]; k++)
        {
            ii++;
            j = ii/3 + 1;
            matrix[i][ii] = h[k]/(mass2[jj]*mass2[j]);
            matrix[ii][i] = matrix[i][ii];
        }
    }

      nfreq = maxvar;
      tred2(matrix, nfreq, eigen, a);
      tqli(eigen,a,nfreq,matrix);
      
//      vibfile = fopen_path(Savebox.path,Savebox.fname,"w");

/*      if (vibfile == NULL)
      {
          message_alert("Error opening vibration output file.","Vibration Setup");
      } */

//      fprintf(vibfile,"\nCurrent Structure\n");

//      for (i=1; i <= natom; i++)
//         fprintf(vibfile,"%d  %d %f %f %f\n",i,atom[i].atomnum,atom[i].x, atom[i].y, atom[i].z);


 //     fprintf(vibfile,"\n\nMass Weighted Eigenvectors\n");

      for (i=0; i < maxvar; i++)
         eigen[i] = factor*sqrt(fabs(eigen[i]));
         
      for (i=0; i < maxvar; i++)
      {
          vnorm = 0.0;
          for (j=0; j < maxvar; j++)
          {
              k = j/3 + 1;
              matrix[j][i] /= mass2[k];
              vnorm += matrix[j][i]*matrix[j][i];
          }
          vnorm = sqrt(vnorm);
          for (j=0; j < maxvar; j++)
             matrix[j][i] /= vnorm;
      }
///
// get overall molecular symmetry
      for (i=1; i <= natom; i++)
      {
          acoord[i][0] = x[i];
          acoord[i][1] = y[i];
          acoord[i][2] = z[i];
      }
      strcpy(ngrp,"   ");
      ptgrp(2,ngrp,acoord);
      // compute dipole moment
      
      for (i=0; i < maxvar; i++)
      {
          dqx = 0.0;
          dqy = 0.0;
          dqz = 0.0;
          dqxx = dqyy = dqzz = 0.0;
          dqxy = dqxz = 0.0;
          dqyx = dqyz = 0.0;
          dqzx = dqzy = 0.0;
          
          for (j=0; j < maxvar; j+= 3)
          {
              strcpy(vsym,"   ");
              k = j/3 + 1;
              acoord[k][0] = x[k] + matrix[j][i];
              acoord[k][1] = y[k] + matrix[j+1][i];
              acoord[k][2] = z[k] + matrix[j+2][i];
              dqx += charge[k]*matrix[j][i];
              dqy += charge[k]*matrix[j+1][i];
              dqz += charge[k]*matrix[j+2][i];
                            
          }

          ptgrp(1,vsym,acoord);
          assign_vibsym(natom,ngrp,vsym,&nsym,acoord);
          vibsym[i] = nsym;
          
      }
      
/*      icycle = 3*natom/5;
      for (i=0; i < (5*icycle); i += 5)
      {    
         fprintf(vibfile,"\nFreq:     %8.3f %8.3f %8.3f %8.3f %8.3f\n",eigen[i],eigen[i+1],eigen[i+2],eigen[i+3],eigen[i+4]);
         fprintf(vibfile,"Sym:         %s       %s       %s       %s       %s\n",symname[vibsym[i]],symname[vibsym[i+1]],symname[vibsym[i+2]],
                   symname[vibsym[i+3]],symname[vibsym[i+4]]);
         fprintf(vibfile,"Int:      %8.4f %8.4f %8.4f %8.4f %8.4f\n",dipole[i],dipole[i+1],dipole[i+2],dipole[i+3],dipole[i+4]);
         fprintf(vibfile,"\n");
      
         for (j=0; j < maxvar; j++)
          fprintf(vibfile,"%-d :       %8.3f %8.3f %8.3f %8.3f %8.3f\n",j,matrix[j][i],matrix[j][i+1],matrix[j][i+2],
             matrix[j][i+3],matrix[j][i+4]);
      }
      if (icycle*5 < 3*natom)
      {
          iz = (3*natom-icycle*5);
          i = icycle*5;
          if (iz == 1)
          {
            fprintf(vibfile,"\nFreq:     %8.3f\n",eigen[i]);
            fprintf(vibfile,"Sym:         %s\n",symname[vibsym[i]]);
            fprintf(vibfile,"Int:      %8.4f\n",dipole[i]);
            fprintf(vibfile,"\n");
      
           for (j=0; j < maxvar; j++)
             fprintf(vibfile,"%-d :       %8.3f %8.3f\n",j,matrix[j][i],matrix[j][i+1]);
          } else if (iz == 2)
          {
            fprintf(vibfile,"\nFreq:     %8.3f %8.3f\n",eigen[i],eigen[i+1]);
            fprintf(vibfile,"Sym:         %s       %s\n",symname[vibsym[i]],symname[vibsym[i+1]]);
            fprintf(vibfile,"Int:      %8.4f %8.4f\n",dipole[i],dipole[i+1]);
            fprintf(vibfile,"\n");
      
           for (j=0; j < maxvar; j++)
             fprintf(vibfile,"%-d :       %8.3f %8.3f\n",j,matrix[j][i],matrix[j][i+1]);
          } else if (iz == 3)
          {
            fprintf(vibfile,"\nFreq:     %8.3f %8.3f %8.3f\n",eigen[i],eigen[i+1],eigen[i+2]);
            fprintf(vibfile,"Sym:         %s       %s       %s\n",symname[vibsym[i]],symname[vibsym[i+1]],symname[vibsym[i+2]]);
            fprintf(vibfile,"Int:      %8.4f %8.4f %8.4f\n",dipole[i],dipole[i+1],dipole[i+2]);
            fprintf(vibfile,"\n");
      
           for (j=0; j < maxvar; j++)
             fprintf(vibfile,"%-d :       %8.3f %8.3f %8.3f\n",j,matrix[j][i],matrix[j][i+1],matrix[j][i+2]);
          } else if (iz == 4)
          {
            fprintf(vibfile,"\nFreq:     %8.3f %8.3f %8.3f %8.3f\n",eigen[i],eigen[i+1],eigen[i+2],eigen[i+3]);
            fprintf(vibfile,"Sym:         %s       %s       %s       %s\n",symname[vibsym[i]],symname[vibsym[i+1]],symname[vibsym[i+2]],
                   symname[vibsym[i+3]]);
            fprintf(vibfile,"Int:      %8.4f %8.4f %8.4f %8.4f\n",dipole[i],dipole[i+1],dipole[i+2],dipole[i+3]);
            fprintf(vibfile,"\n");
      
           for (j=0; j < maxvar; j++)
             fprintf(vibfile,"%-d :       %8.3f %8.3f %8.3f %8.3f\n",j,matrix[j][i],matrix[j][i+1],matrix[j][i+2],
               matrix[j][i+3]);
          }
      }
      fclose(vibfile); */

/*      fprintf(pcmlogfile,"\n==========================================================\n");
      fprintf(pcmlogfile,"Normal Mode Vibrational Analysis\n\n");
      if (symmetry.ishape == 1)
         fprintf(pcmlogfile,"Molecular Shape: LINEAR\n");
      else if (symmetry.ishape == 2)
         fprintf(pcmlogfile,"Molecular Shape: ASYMMETRIC ROTOR\n");
      else if (symmetry.ishape == 3)
         fprintf(pcmlogfile,"Molecular Shape: PROLATE SYMMETRIC ROTOR\n");
      else if (symmetry.ishape == 4)
         fprintf(pcmlogfile,"Molecular Shape: OBLATE SYMMETRIC ROTOR\n");
      else if (symmetry.ishape == 5)
         fprintf(pcmlogfile,"Molecular Shape: SPHERICAL ROTOR\n");
      else if (symmetry.ishape == 6)
         fprintf(pcmlogfile,"Molecular Shape: PLANAR\n");
      
      fprintf(pcmlogfile,"Symmetry Point Group: %s\n",ngrp);
      fprintf(pcmlogfile,"Symmetry Number: %d\n",symmetry.numsym);
      if (symmetry.chiral)
          fprintf(pcmlogfile,"Chirality: Yes\n");
      else
          fprintf(pcmlogfile,"Chirality: No\n"); */

      if (symmetry.ishape == 1)
         mm = maxvib - 5;
      else
         mm = maxvib - 6;

/*      fprintf(pcmlogfile,"==========================================================\n");
      fprintf(pcmlogfile,"\nFundmental Normal Vibrational Frequencies\n");
      for (i=0; i < mm; i++)
         fprintf(pcmlogfile,"   %-3d            %8.3f      %s    %8.3f\n",i,eigen[i],symname[vibsym[i]],dipole[i]);
         
      fprintf(pcmlogfile,"\nMoments of Inertia\n");
      fprintf(pcmlogfile," IX = %8.3f   IY = %8.3f   IZ = %8.3f\n",symmetry.xi,symmetry.yi,symmetry.zi); */
            
      strcpy(vibdata.ptgrp,ngrp);
      vibdata.mom_ix = symmetry.xi;
      vibdata.mom_iy = symmetry.yi;
      vibdata.mom_iz = symmetry.zi;

      if (symmetry.chiral)
         imix = 2;
      else
         imix = 1;          
//  set values
      temp = 298.15;
      mass = 0.0;
      for (i=1; i <= natom; i++)
         mass += atomwt[i];
      efact = 1.987;
      cnosym = symmetry.numsym;
      cimix = (double)imix;
      etrans=htrans=strans=gtrans=cptrans=0.0;
      erot=hrot=srot=grot=cprot=0.0;
      evib=hvib=svib=gvib=cpvib=0.0;
      emix=hmix=smix=gmix=cpmix=0.0;
      etot=htot=stot=gtot=cptot=0.0;

//  compute entropy, internal energy
      strans = 2.2868*(5.0*log10(temp)+3.0*log10(mass))-2.3135;
      etrans = 0.001987*1.5*temp;
      if (strans < 0.0) strans = 0.0;
      if (symmetry.ishape == 1)
      {
         xyzi = sqrt(symmetry.xi*symmetry.xi+symmetry.yi*symmetry.yi+symmetry.zi*symmetry.zi);
         srot = 4.5736*(log10(temp)+log10(xyzi)-0.6049)+ efact;
         erot = 0.001987*temp;
      }else
      {
         xyzi = symmetry.xi*symmetry.yi*symmetry.zi;
         srot = 2.2868*(3.0*log10(temp)+log10(xyzi)-1.3176)+ 2.9805;
         erot = 0.001987*temp*1.5;
      }
      if (srot < 0.0) srot = 0.0;
      ssym = -2.2868*2.0*log10(cnosym);
      srot += ssym;

      qptf = 0.0;
      ezpe = 0.0;
      entrov1 = 0.0;
      entrov2 = 0.0;
      for (i=0; i < mm; i++)
      {
          wdt = fabs(eigen[i])/temp;
          efwt = exp(-1.43893*wdt);
          qptf += efwt;
          entrov1 += eigen[i]*efwt/(1.0-efwt);
          entrov2 += log(1.0-efwt);

          evib += eigen[i]*0.0028591*efwt/(1.0-efwt);
          ezpe += 0.5*eigen[i]*0.0028591;
      }
      svib = efact*(1.438903*entrov1/temp - entrov2);
      smix = 2.2868*2.0*log10(cimix);
      stot = strans + srot + svib + smix;
      evib += ezpe;
      etot = etrans + erot + evib + emix + get_total_energy();
      
// compute enthalpy
      htrans = etrans + 0.001987*temp;
      hrot = erot;
      hvib = evib;
      hmix = emix;
      htot = htrans + hrot + hvib + hmix +  get_total_energy();
// compute free energy
      gtrans = htrans - temp*strans*0.001;
      grot = hrot - temp*srot*0.001;
      gvib = hvib - temp*svib*0.001;
      gmix = hmix - temp*smix*0.001;
      gtot = gtrans + grot + gvib + gmix +  get_total_energy();
//  compute heat capacity
      cptrans = 1.987 + 1.987*1.5;
      if (symmetry.ishape == 1)
         cprot = 1.987;
      else
         cprot = 1.987*1.5;
      cpfact = 4.111/(temp*temp);
      for (i=0; i < mm; i++)
      {
          wdt = fabs(eigen[i])/temp;
          efwt = exp(-1.43893*wdt);
          efwt2 = 1.0-efwt;
          cpvib += eigen[i]*eigen[i]*cpfact*efwt/(efwt2*efwt2);
      }
      cptot = cptrans + cprot + cpvib + cpmix;
//      
      vibdata.etot = etot;
      vibdata.htot = htot;
      vibdata.stot = stot;
      vibdata.gtot = gtot;
      vibdata.cptot = cptot;
// print out
/*      fprintf(pcmlogfile,"\nTemperatur: %8.3f K \n",temp);
      fprintf(pcmlogfile,"\nZero Point Energy: %8.3f kcal/mol \n",ezpe);
      fprintf(pcmlogfile,"\n                    Energy       Enthalpy       Entropy      Free Energy      Heat Capacity\n");
      fprintf(pcmlogfile,"\n                    kcal/mol     kcal/mol       eu           kcal/mol          cal/mol/deg\n");
      fprintf(pcmlogfile,"Translational: %10.3f     %10.3f     %10.3f    %10.3f       %10.3f\n",etrans,htrans,strans,gtrans,cptrans);
      fprintf(pcmlogfile,"Rotational:    %10.3f     %10.3f     %10.3f    %10.3f       %10.3f\n",erot,hrot,srot,grot,cprot);
      fprintf(pcmlogfile,"Vibrational:   %10.3f     %10.3f     %10.3f    %10.3f       %10.3f\n",evib,hvib,svib,gvib,cpvib);
      fprintf(pcmlogfile,"Mixing:        %10.3f     %10.3f     %10.3f    %10.3f       %10.3f\n",emix,hmix,smix,gmix,cpmix);
      fprintf(pcmlogfile,"Total:         %10.3f     %10.3f     %10.3f    %10.3f       %10.3f\n",etot,htot,stot,gtot,cptot); */
      
//      
    free_dvector(eigen, 0, 3*natom);
    free_dmatrix(vects, 0, 3*natom, 0, 3*natom);
    free_ivector(vibsym ,0,maxvar);
    free_ivector(hinit ,0,maxvar);
    free_ivector(hstop ,0,maxvar);
    free_ivector(hindex ,0,maxhess);
    free_dvector(hdiag ,0,maxvar);
    free_dvector(h ,0,maxhess);

    free_dmatrix(matrix ,0, maxvib, 0, maxvib);
    free_dvector(a ,0, maxvib);
    free_dvector(mass2, 0 , natom+1);
    free_dvector(dipole,0,maxvar);
}
/* =============================================== */
void assign_vibsym(int natom, char *ngrp,char *vsym,int *nsym,double **a)
{
    int iz;
    double  **b;
    
    b = dmatrix (0, natom+1, 0,3);
        
    *nsym = 0;
    if (strcmp(ngrp,"Cs") == 0)
    {
        if (strcmp(vsym,"Cs") == 0)
           *nsym = 1; // A'
        else
           *nsym = 2; // a''
    } else if (strcmp(ngrp,"Ci") == 0)
    {
        if (strcmp(vsym,"Ci") == 0)
           *nsym = 3; // ag
        else
           *nsym = 4; // au
    } else if (strcmp(ngrp,"C1") == 0)
    {
        *nsym = 0; // a
    } else if (strcmp(ngrp,"C2") == 0)
    {
        if (strcmp(vsym,"C2") == 0)
           *nsym = 0; // a
        else
           *nsym = 15; // b
    } else if (strcmp(ngrp,"C3") == 0)
    {
        if (strcmp(vsym,"C3") == 0)
           *nsym = 0; //  a
        else
           *nsym = 27; // e
    } else if (strcmp(ngrp,"C4") == 0)
    {
        if (strcmp(vsym,"C4") == 0)
           *nsym = 0; // a
        else 
        {
            if (strcmp(vsym,"C2") == 0)
               *nsym = 15; // b
            else
               *nsym = 27; // e
        }
    } else if (strcmp(ngrp,"C5") == 0)
    {
        if (strcmp(vsym,"C5") == 0)
           *nsym = 0; // a
        else
           *nsym = 28; // e?
    } else if (strcmp(ngrp,"C6") == 0)
    {
        if (strcmp(vsym,"C6") == 0)
           *nsym = 0; // a
        else
        {
            if (strcmp(vsym,"C3") == 0)
               *nsym = 15; // b
            else
               *nsym = 27; // e
        }
    } else if (strcmp(ngrp,"C2v") == 0)
    {
        if (strcmp(vsym,"C2v") == 0)
           *nsym = 5; // a1
        else
        {
            if (strcmp(vsym,"C2") == 0)
               *nsym = 6; //a2
            else
            {
               iz = findh(a,b,3);
               if (iz)
                   *nsym = 16; // b1
               else
                   *nsym = 17; // b2
            }
        }
    } else if (strcmp(ngrp,"C3v") == 0)
    {
        if (strcmp(vsym,"C3v") == 0)
           *nsym = 5; // a1
        else
        {
            if (strcmp(vsym,"C3") == 0)
               *nsym = 6; // a2
            else
               *nsym = 27; // e
        }
    } else if (strcmp(ngrp,"C4v") == 0)
    {
        if (strcmp(vsym,"C4v") == 0)
           *nsym = 5; // a1
        else
        {
            if (strcmp(vsym,"C4") == 0)
               *nsym = 6; // a2
            else
            {
                if (strcmp(vsym,"C2v") == 0)
                   *nsym = 16; //b?
                else
                   *nsym = 27; // e
            }
        }
    } else if (strcmp(ngrp,"C5v") == 0)
    {
        if (strcmp(vsym,"C5v") == 0)
           *nsym = 5; // a1
        else
        {
            if (strcmp(vsym,"C5") == 0)
               *nsym = 6; // a2
            else
               *nsym = 28; // e?
        }
    } else if (strcmp(ngrp,"C6v") == 0)
    {
        if (strcmp(vsym,"C6v") == 0)
           *nsym = 5; // a1
        else
        {
            if (strcmp(vsym,"C6") == 0)
               *nsym = 6; // a2
            else
            {
                if (strcmp(vsym,"C3v") == 0)
                   *nsym = 16; // b?
                else
                   *nsym = 28; // e?
            }
        }
    } else if (strcmp(ngrp,"C2h") == 0)
    {
        if (strcmp(vsym,"C2h") == 0)
           *nsym = 3; // ag
        else
        {
            if (strcmp(vsym,"C2") == 0)
               *nsym = 4; // au
            else
            {
                if (strcmp(vsym,"Ci") == 0)
                   *nsym = 19; // bg
                else
                   *nsym = 20; // bu
            }
        }
    } else if (strcmp(ngrp,"C3h") == 0)
    {
        if (strcmp(vsym,"C3h") == 0)
           *nsym = 1; // a'
        else
        {
            if (strcmp(vsym,"C3") == 0)
               *nsym = 2; // a"
            else
            {
                if (strcmp(vsym,"Cs") == 0)
                   *nsym = 31; // e' 
                else
                   *nsym = 32; // e"
            }
        }
    } else if (strcmp(ngrp,"C4h") == 0)
    {
        if (strcmp(vsym,"C4h") == 0)
           *nsym = 3; // ag
        else
        {
            if (strcmp(vsym,"C2h") == 0)
               *nsym = 19; // Bg
            else
            {
                if (strcmp(vsym,"C4") == 0)
                   *nsym = 4; // au 
                else
                {
                    if (strcmp(vsym,"S4") == 0)
                        *nsym = 20; // bu
                    else
                    {
                        iz = findi(a,b);
                        if (iz)
                           *nsym = 33 ; // eg
                        else
                           *nsym = 34 ; // eu
                    }
                }
            }
        }
    } else if (strcmp(ngrp,"C5h") == 0)
    {
        if (strcmp(vsym,"C5h") == 0)
           *nsym = 1; // a'
        else
        {
            if (strcmp(vsym,"C5") == 0)
               *nsym = 2; // a"
            else
            {
                if (strcmp(vsym,"Cs") == 0)
                   *nsym = 35; // e?' 
                else
                   *nsym = 36; // e?"
            }
        }
    } else if (strcmp(ngrp,"C6h") == 0)
    {
        if (strcmp(vsym,"C6h") == 0)
           *nsym = 3; // ag
        else
        {
            if (strcmp(vsym,"C6") == 0)
               *nsym = 4; // au
            else
            {
                if (strcmp(vsym,"S6") == 0)
                   *nsym = 19; // bg 
                else
                {
                    if (strcmp(vsym,"C3h") == 0)
                        *nsym = 20; // bu
                    else
                    {
                        if (strcmp(vsym,"C2") == 0)
                        {
                           iz = findi(a,b);
                           if (iz)
                             *nsym = 41 ; // e2g
                           else
                              *nsym = 42 ; // e2u
                        } else
                        {
                            iz = findi(a,b);
                            if (iz)
                              *nsym = 39; // e1g
                            else
                              *nsym = 40; // e1u
                        }
                    }
                }
            }
        }
    } else if (strcmp(ngrp,"D2") == 0)
    {
        if (strcmp(vsym,"D2") == 0)
          *nsym = 0; // a
        else
        {
            if (strcmp(vsym,"C2") == 0)
               *nsym = 16; // b1 or b2
        }
    } else if (strcmp(ngrp,"D3") == 0)
    {
        if (strcmp(vsym,"D3") == 0)
          *nsym = 5; // a1
        else
        {
            if (strcmp(vsym,"C3") == 0)
               *nsym = 6; // a2
            else
               *nsym = 27; // e
        }
    } else if (strcmp(ngrp,"D4") == 0)
    {
        if (strcmp(vsym,"D4") == 0)
          *nsym = 5; // a1
        else
        {
            if (strcmp(vsym,"C4") == 0)
               *nsym = 6; // a2
            else
            {
                if (strcmp(vsym,"C1") == 0)
                   *nsym = 27; // e
                else
                {
                    if (strcmp(vsym,"D2") == 0)
                       *nsym = 16; // b1
                    else
                       *nsym = 17; // b2
                }
            }
        }
    } else if (strcmp(ngrp,"D5") == 0)
    {
        if (strcmp(vsym,"D5") == 0)
          *nsym = 5; // a1
        else
        {
            if (strcmp(vsym,"C5") == 0)
               *nsym = 6; // a2
            else
               *nsym = 28; // e?
        }
    } else if (strcmp(ngrp,"D6") == 0)
    {
        if (strcmp(vsym,"D6") == 0)
          *nsym = 5; // a1
        else
        {
            if (strcmp(vsym,"C6") == 0)
               *nsym = 6; // a2
            else
            {
                if (strcmp(vsym,"D3") == 0)
                   *nsym = 16; // b1
                else
                {
                    if (strcmp(vsym,"C3") == 0)
                       *nsym = 17; // b2
                    else
                    {
                        if (strcmp(vsym,"C2") == 0)
                           *nsym = 29; // e2
                        else
                           *nsym = 28; // e1
                    }
                }
            }
        }
    } else if (strcmp(ngrp,"D2h") == 0)
    {
        if (strcmp(vsym,"D2h") == 0)
          *nsym = 3; // ag
        else
        {
            if (strcmp(vsym,"D2") == 0)
               *nsym = 4; // au
            else
            {
                iz = findi(a,b);
                if (iz)
                   *nsym = 21; // b1g
                else
                   *nsym = 22; // b1u?
            }
        }
    } else if (strcmp(ngrp,"D3h") == 0)
    {
        if (strcmp(vsym,"D3h") == 0)
          *nsym = 5; // a1
        else
        {
            if (strcmp(vsym,"C3h") == 0)
               *nsym = 6; // a2
            else
            {
               if (strcmp(vsym,"Cs") == 0)
                  *nsym = 31; // e'
               else
               {
                  if (strcmp(vsym,"D3") == 0)
                     *nsym = 5; // a1
                  else
                  {
                     if (strcmp(vsym,"C3v") == 0)
                        *nsym = 10; // a2"
                     else
                        *nsym = 32; // e"
                  }
               }
            }
        }
    } else if (strcmp(ngrp,"D4h") == 0)
    {
        if (strcmp(vsym,"D4h") == 0)
          *nsym = 11; // a1g
        else
        {
            if (strcmp(vsym,"C4h") == 0)
               *nsym = 13; // a2g
            else
            {
               if (strcmp(vsym,"D4") == 0)
                  *nsym = 12; // a1u
               else
               {
                  if (strcmp(vsym,"C4v") == 0)
                     *nsym = 14; // a2u
                  else
                  {
                      iz = findi(a,b);
                      if (iz)
                      {
                          if (strcmp(vsym,"D2h") == 0)
                            *nsym = 21; // b1g
                          else
                            *nsym = 33; // eg
                      } else
                      {
                          if (strcmp(vsym,"D2d") == 0)
                             *nsym = 22; // b?u
                          else
                             *nsym = 34; //eu
                      }                              
                  }
               }
            }
        }
    } else if (strcmp(ngrp,"D5h") == 0)
    {
        if (strcmp(vsym,"D5h") == 0)
          *nsym = 7; // a1'
        else
        {
            if (strcmp(vsym,"C5h") == 0)
               *nsym = 8; // a2'
            else
            {
               if (strcmp(vsym,"D5") == 0)
                  *nsym = 9; // a1''
               else
               {
                  if (strcmp(vsym,"C5v") == 0)
                     *nsym = 10; // a2''
                  else
                     *nsym = 36; // e"
               }
            }
        }
    } else if (strcmp(ngrp,"D6h") == 0)
    {
        if (strcmp(vsym,"D6h") == 0)
          *nsym = 11; // a1g
        else
        {
            if (strcmp(vsym,"C6h") == 0)
               *nsym = 13; // a2g
            else
            {
               if (strcmp(vsym,"D3d") == 0)
                  *nsym = 21; // b?g
               else
               {
                  if (strcmp(vsym,"D6") == 0)
                     *nsym = 12; // a1u
                  else
                  {
                     if (strcmp(vsym,"C6v") == 0)
                        *nsym = 14; // a2u
                     else
                     {
                        if (strcmp(vsym,"D3h") == 0)
                           *nsym = 22; // b?u
                        else
                        {
                            if ( strcmp(vsym,"C2") == 0 || strcmp(vsym,"D2") == 0 ||
                                 strcmp(vsym,"C2h") == 0 )
                            {
                                iz = findi(a,b);
                                if (iz)
                                   *nsym = 41; // e2g
                                else
                                   *nsym = 42; // e2u
                            } else
                            {
                                iz = findi(a,b);
                                if (iz)
                                   *nsym = 39; // e1g
                                else
                                   *nsym = 40; // e1u
                            }
                        }
                     }
                  }
               }
            }
        }
    } else if (strcmp(ngrp,"D2d") == 0)
    {
        if (strcmp(vsym,"D2d") == 0 )
          *nsym = 5; // a1
        else
        {
           if (strcmp(vsym,"S4") == 0)
             *nsym = 6; // a2
           else
           {
              if (strcmp(vsym,"D2") == 0) 
                *nsym = 16; // b1
              else
              {
                 if (strcmp(vsym,"C2v") == 0)
                   *nsym = 17; // b2
                 else
                   *nsym = 27; // e
              }
           }
        }
    } else if (strcmp(ngrp,"D3d") == 0)
    {
        iz = findi(a,b);
        if (iz)
        {
           if (strcmp(vsym,"D3d") == 0 )
             *nsym = 11; // a1g
           else
           {
              if (strcmp(vsym,"S6") == 0 )
                *nsym = 13; // a2g
              else
                *nsym = 33; // eg
           }
        } else
        {
           if (strcmp(vsym,"D3") == 0 )
             *nsym = 12; // a1u
           else
           {
              if (strcmp(vsym,"C3v") == 0 )
                *nsym = 14; // a2u
              else
                *nsym = 34; // eu
           }
        }
    } else if (strcmp(ngrp,"D4d") == 0)
    {
        if (strcmp(vsym,"D4d") == 0 )
          *nsym = 5; // a1
        else
        {
           if (strcmp(vsym,"S8") == 0)
             *nsym = 6; // a2
           else
           {
              if (strcmp(vsym,"D4") == 0) 
                *nsym = 16; // b1
              else
              {
                 if (strcmp(vsym,"C4v") == 0)
                   *nsym = 17; // b2
                 else
                   *nsym = 27; // e
              }
           }
        }
    } else if (strcmp(ngrp,"D5d") == 0)
    {
        iz = findi(a,b);
        if (iz)
        {
           if (strcmp(vsym,"D5d") == 0 )
             *nsym = 11; // a1g
           else
           {
              if (strcmp(vsym,"S10") == 0 )
                *nsym = 13; // a2g
              else
                *nsym = 33; // eg
           }
        } else
        {
           if (strcmp(vsym,"D5") == 0 )
             *nsym = 12; // a1u
           else
           {
              if (strcmp(vsym,"C5v") == 0 )
                *nsym = 14; // a2u
              else
                *nsym = 34; // eu
           }
        }
    } else if (strcmp(ngrp,"D6d") == 0)
    {
        if (strcmp(vsym,"D6d") == 0 )
          *nsym = 5; // a1
        else
        {
           if (strcmp(vsym,"S12") == 0)
             *nsym = 6; // a2
           else
           {
              if (strcmp(vsym,"D6") == 0) 
                *nsym = 16; // b1
              else
              {
                 if (strcmp(vsym,"C6v") == 0)
                   *nsym = 17; // b2
                 else
                   *nsym = 27; // e
              }
           }
        }
    } else if (strcmp(ngrp,"S4") == 0)
    {
        if (strcmp(vsym,"S4") == 0 )
          *nsym = 0; // a
        else
        {
           if (strcmp(vsym,"C2") == 0) 
             *nsym = 15; // b
           else
             *nsym = 27; // e
        }
    } else if (strcmp(ngrp,"S6") == 0)
    {
        iz = findi(a,b);
        if (iz)
        {
           if (strcmp(vsym,"S6") == 0 )
             *nsym = 3; // ag
           else
             *nsym = 33; // eg
        } else
        {
           if (strcmp(vsym,"C3") == 0 )
             *nsym = 4; // au
           else
             *nsym = 34; // eu
        }
    } else if (strcmp(ngrp,"S8") == 0)
    {
        if (strcmp(vsym,"S8") == 0 )
          *nsym = 0; // a
        else
        {
           if (strcmp(vsym,"C4") == 0) 
             *nsym = 15; // b
           else
             *nsym = 27; // e
        }
    } else if (strcmp(ngrp,"T") == 0)
    {
        if (strcmp(vsym,"T") == 0 || strcmp(vsym,"C3") == 0)
          *nsym = 0; // a
        else
        {
           if (strcmp(vsym,"D2") == 0 || strcmp(vsym,"C2") == 0)
             *nsym = 27; // e
           else
             *nsym = 43; // t
        }
    } else if (strcmp(ngrp,"Th") == 0)
    {
        if (strcmp(vsym,"Th") == 0)
          *nsym = 3; // ag
        else
        {
           if (strcmp(vsym,"T") == 0)
             *nsym = 4; // au
           else
           {
              iz = findi(a,b);
              if (iz)
              {
                  if (strcmp(vsym,"D2h") == 0 || strcmp(vsym,"D2") == 0)
                     *nsym = 33; // Eg
                  else
                     *nsym = 44; // Tg
              } else
              {
                 if (strcmp(vsym,"D2") == 0 || strcmp(vsym,"C2") == 0)
                   *nsym = 34; // eu
                 else
                   *nsym = 45; // Tu
              }
           }
        }   
    } else if (strcmp(ngrp,"Td") == 0)
    {
        if (strcmp(vsym,"Td") == 0)
          *nsym = 5; // a1
        else
        {
           if (strcmp(vsym,"T") == 0)
             *nsym = 6; // a2
           else
           {
              if (strcmp(vsym,"D2") == 0 || strcmp(vsym,"D2d") == 0)
                *nsym = 27; // e
              else
                *nsym = 43; // T
           }
        }   
    } else if (strcmp(ngrp,"O") == 0)
    {
        if (strcmp(vsym,"O") == 0)
          *nsym = 5; // a1
        else
        {
           if (strcmp(vsym,"T") == 0)
             *nsym = 6; // a2
           else
           {
              if (strcmp(vsym,"D2") == 0)
                *nsym = 27; // e
              else
                *nsym = 43; // T
           }
        }   
    }
    free_dmatrix(b ,0, natom+1, 0,3);
}    
// ======================
void tred2(double **a,int n, double *d, double *e)
{

    int i,j,k,l, nn;
    double scale, hh, h, g, f;

    nn = n-1;
    for (i= nn; i >= 1; i--)
    {
        l = i - 1;
        h = scale = 0.0;
        if (l > 0)
        {
            for (k=0; k <= l; k++)
                scale += fabs(a[i][k]);
            if (scale == 0.0)
                e[i] = a[i][l];
            else
            {
                for (k=0; k <= l; k++)
                {
                    a[i][k] /= scale;
                    h += a[i][k]*a[i][k];
                }
                f = a[i][l];
                g = f>0 ? -sqrt(h) : sqrt(h);
                e[i] = scale*g;
                h -= f*g;
                a[i][l] = f-g;
                f = 0.0;
                for (j=0; j <= l; j++)
                {
                    a[j][i] = a[i][j]/h;
                    g = 0.0;
                    for (k=0; k <= j; k++)
                        g += a[j][k]*a[i][k];
                    for (k = j+1; k <= l; k++)
                        g += a[k][j]*a[i][k];
                    e[j] = g/h;
                    f += e[j]*a[i][j];
                }
                hh = f/(h+h);
                for (j=0; j <= l; j++)
                {
                    f = a[i][j];
                    e[j] = g = e[j] - hh*f;
                    for (k=0; k <= j; k++)
                      a[j][k] -= (f*e[k]+g*a[i][k]);
                }
            }
        }else
           e[i] = a[i][l];
        d[i] = h;
    }
    d[0] = 0.0;
    e[0] = 0.0;
    for (i=0; i < n; i++)
    {
        l = i - 1;
        if (i > 0)
        {
            for (j=0; j <= l; j++)
            {
                g = 0.0;
                for (k=0; k <= l; k++)
                    g += a[i][k]*a[k][j];
                for (k=0; k <= l; k++)
                    a[k][j] -= g*a[k][i];
            }
        }
        d[i] = a[i][i];
        a[i][i] = 1.0;
        for (j=0; j <= l; j++)
        {
            a[j][i] = a[i][j] = 0.0;
        }
    }
}
/* =========================================== */
void tqli(double *d, double *e, int n, double **z)
{
    int m,l,iter,i,k,j;
    double s,r,p,g,f,dd,c,b;

    for (i=1; i < n; i++)
       e[i-1] = e[i];
    e[n-1] = 0.0;

    for (l=0; l < n; l++)
    {
        iter = 0;
        do
        {
            for (m=l; m < n-1; m++)
            {
                dd = fabs(d[m]) + fabs(d[m+1]);
                if (fabs(e[m])+dd == dd) break;
            }
            if (m!= l)
            {
                if (iter++ == 30)
                {
                    fprintf(pcmlogfile,"Error in tqli\n");
                    message_alert("Error in Tqli","Error");
                    return;
                }
                g = (d[l+1]-d[l])/(2.0*e[l]);
                r = sqrt((g*g)+1.0);
                g = d[m]-d[l]+e[l]/(g+SIGN(r,g));
                s = c = 1.0;
                p = 0.0;
                for (i=m-1; i >= l; i--)
                {
                    f = s*e[i];
                    b = c*e[i];
                    if (fabs(f) >= fabs(g))
                    {
                        c = g/f;
                        r = sqrt((c*c)+1.0);
                        e[i+1] = f*r;
                        c *= (s=1.0/r);
                    } else
                    {
                        s = f/g;
                        r = sqrt((s*s)+1.0);
                        e[i+1] = g*r;
                        s *= (c=1.0/r);
                    }
                    g = d[i+1]-p;
                    r = (d[i]-g)*s+2.0*c*b;
                    p = s*r;
                    d[i+1] = g+p;
                    g = c*r-b;
                    for (k=0; k < n; k++)
                    {
                        f = z[k][i+1];
                        z[k][i+1] = s*z[k][i]+c*f;
                        z[k][i] = c*z[k][i]-s*f;
                    }
                }
                d[l] = d[l]-p;
                e[l] = g;
                e[m] = 0.0;
            }
        } while (m!= l);
    }
// sort d and z
    for (i=0; i < n; i++)
    {
        p  = d[k=i];
        for (j=i+1; j < n; j++)
            if (d[j] >= p) p = d[k=j];
        if (k != i)
        {
            d[k] = d[i];
            d[i] = p;
            for (j=0; j < n; j++)
            {
                p = z[j][i];
                z[j][i] = z[j][k];
                z[j][k] = p;
            }
        }
    }  
}
