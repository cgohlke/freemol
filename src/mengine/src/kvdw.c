#define EXTERN extern
#include "pcwin.h"
#include "nonbond.h"
#include "atom_k.h"

void numeral(int,char *,int);
int get_field(void);
char * get_radiustype(void);
char * get_radiussize(void);
char * get_radiusrule(void);
char * get_epsrule(void);
void kvdw(int natom,int *type,int *atomnum);
       
EXTERN struct t_vdw1 {
        int  nvdw;
        float rad[MAXVDWCONST], eps[MAXVDWCONST];
        int lpd[MAXVDWCONST], ihtyp[MAXVDWCONST], ihdon[MAXVDWCONST];
        float alpha[MAXVDWCONST],n[MAXVDWCONST],a[MAXVDWCONST],g[MAXVDWCONST];
        char da[MAXVDWCONST][2];
        } vdw1;

EXTERN struct t_vdwpr_k {
        int nvdwpr;
        int  ia1[MAXBONDCONST],ia2[MAXBONDCONST];
        char kv[MAXBONDCONST][7];
        float radius[MAXBONDCONST], eps[MAXBONDCONST];
        } vdwpr_k;

struct t_mm3hbond {
        int npair;
        int htype[MAXBONDCONST];
        int otype[MAXBONDCONST];
        } mm3hbond;
        
EXTERN int Missing_constants;
        
void kvdw(int natom,int *type,int *atomnum)
{
  int i, j, k,field;
   int it, ita, itb, found, ianum,ibnum;
   int dapair;
   char pa[4],pb[4],pt[7];
   float rad_ia, rad_ib;
   float eps1, eps2, seps1, seps2;
   float rii,rjj,gammaij, rij, epsij;
   
   field = get_field();
//                
   nonbond.npair = 0;
   for (i=0; i <= nonbond.maxnbtype; i++)
   {
       for (j=0; j <= nonbond.maxnbtype; j++)
           nonbond.iNBtype[j][i] = 0;
   }
   
   for (i=1; i < natom; i++)
   {
       ita = type[i];
       ianum = atomnum[i];
       for (j= i+1; j <= natom; j++)
       {
           itb = type[j];
           ibnum = atomnum[j];
           if (ita < itb)
              it = 1000*ita + itb;
           else
              it = 1000*itb + ita;
              
           found = FALSE;
           if (nonbond.iNBtype[ita][itb] != 0)
              found = TRUE;

           if (found == FALSE)
           {
               if (field == MMFF94)
               {
                        dapair = FALSE;
                        rii = vdw1.a[ita]*pow(vdw1.alpha[ita],0.25);
                        rjj = vdw1.a[itb]*pow(vdw1.alpha[itb],0.25);
                        if ((atomnum[i] == 1 && ita != 5) || (atomnum[j] ==1  && itb != 5) )  // polar hydrogens
                        {
                                rij = 0.5*(rii+rjj);
                                if ( (strcmp(vdw1.da[ita],"A") == 0) && (strcmp(vdw1.da[itb],"D") == 0) )
                                {
                                        dapair = TRUE;
                                }
                                if ( (strcmp(vdw1.da[itb],"A") == 0) && (strcmp(vdw1.da[ita],"D") == 0) )
                                {
                                        dapair = TRUE;
                                }
                        } else
                        {
                                gammaij = (rii-rjj)/(rii+rjj);
                                rij = 0.5*(rii+rjj)*(1.00+0.2*(1.00- exp(-12.0*gammaij*gammaij)));
                        }
                        epsij = ((181.16*vdw1.g[ita]*vdw1.g[itb]*vdw1.alpha[ita]*vdw1.alpha[itb])/
                           ( sqrt(vdw1.alpha[ita]/vdw1.n[ita]) + sqrt(vdw1.alpha[itb]/vdw1.n[itb])) ) / pow(rij,6.0);

                        if (dapair == TRUE)
                        {
                                rij *= 0.8;
                                epsij *= 0.5;
                        }

                   nonbond.vrad[ita][itb] = rij;
                   nonbond.veps[ita][itb] = epsij;
                   nonbond.ipif[ita][itb] = 1.0;
                   nonbond.vrad[itb][ita] = rij;
                   nonbond.veps[itb][ita] = epsij;
                   nonbond.ipif[itb][ita] = 1.0;

                   nonbond.vrad14[ita][itb] = nonbond.vrad[ita][itb];
                   nonbond.vrad14[itb][ita] = nonbond.vrad[itb][ita];
                   nonbond.veps14[ita][itb] = nonbond.veps[ita][itb];
                   nonbond.veps14[itb][ita] = nonbond.veps[itb][ita];

               } else
               {
                  rad_ia = vdw1.rad[ita];
                  rad_ib = vdw1.rad[itb];
                  if (strcmp(get_radiustype(),"SIGMA") == 0)
                  {
                      rad_ia *= 1.122462048;
                      rad_ib *= 1.122462048;
                  }  
                  if (strcmp(get_radiussize(),"DIAMETER") == 0)
                  {
                      rad_ia /= 2.0;
                      rad_ib /= 2.0;
                  }
                  eps1 = vdw1.eps[ita];
                  eps2 = vdw1.eps[itb];
       
                  if (strcmp(get_radiusrule(),"ARITHMETIC") == 0)
                  {
                     nonbond.vrad[ita][itb] = rad_ia + rad_ib;
                     nonbond.vrad[itb][ita] = rad_ia + rad_ib;
                  } else if (strcmp(get_radiusrule(),"GEOMETRIC") == 0)
                  {
                     nonbond.vrad[ita][itb] = 2.0*( sqrt(rad_ia)*sqrt(rad_ib));
                     nonbond.vrad[itb][ita] = 2.0*( sqrt(rad_ia)*sqrt(rad_ib));
                  } else if (strcmp(get_radiusrule(),"CUBIC-MEAN") == 0)
                  {
                     if ((atomnum[i] == 1 && ita != 5) || (atomnum[j] ==1  && itb != 5) )  // polar hydrogens
                     {
                        nonbond.vrad[ita][itb] = rad_ia + rad_ib;
                        nonbond.vrad[itb][ita] = rad_ia + rad_ib;
                     }else
                     {
                        nonbond.vrad[ita][itb] = 2.0*( ((rad_ia*rad_ia*rad_ia) + (rad_ib*rad_ib*rad_ib))/(rad_ia*rad_ia + rad_ib*rad_ib) );
                        nonbond.vrad[itb][ita] = 2.0*( ((rad_ia*rad_ia*rad_ia) + (rad_ib*rad_ib*rad_ib))/(rad_ia*rad_ia + rad_ib*rad_ib) );
                     }
                  } else
                  {
                     nonbond.vrad[ita][itb] = rad_ia + rad_ib;
                     nonbond.vrad[itb][ita] = rad_ia + rad_ib;
                  }
      
                  if (strcmp(get_epsrule(),"ARITHMETIC") == 0)
                  {
                      nonbond.veps[ita][itb] = 0.5*(eps1 + eps2);
                      nonbond.veps[itb][ita] = 0.5*(eps1 + eps2);
                  } else if (strcmp(get_epsrule(),"GEOMETRIC") == 0)
                  {
                      nonbond.veps[ita][itb] = sqrt( eps1*eps2 );
                      nonbond.veps[itb][ita] = sqrt( eps1*eps2 );
                  } else if (strcmp(get_epsrule(),"HARMONIC") == 0)
                 {
                      nonbond.veps[ita][itb] = 2.0* (eps1*eps2)/(eps1+eps2);
                      nonbond.veps[itb][ita] = 2.0* (eps1*eps2)/(eps1+eps2);
                  } else if (strcmp(get_epsrule(),"HHG") == 0)
                  {
                     if ((atomnum[i] == 1 && ita != 5) || (atomnum[j] ==1  && itb != 5) )  // polar hydrogens
                     {
                       nonbond.veps[ita][itb] = 0.5*(eps1 + eps2);
                       nonbond.veps[itb][ita] = 0.5*(eps1 + eps2);
                     } else
                     {
                        seps1 = sqrt(eps1);
                        seps2 = sqrt(eps2);
                        nonbond.veps[ita][itb] = 4.0*(eps1*eps2)/( (seps1+seps2)*(seps1+seps2));
                        nonbond.veps[itb][ita] = 4.0*(eps1*eps2)/( (seps1+seps2)*(seps1+seps2));
                     }
                  }else
                  {
                      nonbond.veps[ita][itb] = sqrt( eps1*eps2 );
                      nonbond.veps[itb][ita] = sqrt( eps1*eps2 );
                  }

                  if (field == MMX)
                  {
                      nonbond.ipif[ita][itb] = 1;
                      nonbond.ipif[itb][ita] = 1;
                      if( (ita == 2 || ita == 4 || ita == 40 || ita == 48) && 
                          (itb == 2 || itb == 4 || itb == 40 || itb == 48) )
                          {
                               nonbond.ipif[ita][itb] = 0;
                               nonbond.ipif[itb][ita] = 0;
                          }
                  }
                   nonbond.vrad14[ita][itb] = nonbond.vrad[ita][itb];
                   nonbond.vrad14[itb][ita] = nonbond.vrad[itb][ita];
                   nonbond.veps14[ita][itb] = nonbond.veps[ita][itb];
                   nonbond.veps14[itb][ita] = nonbond.veps[itb][ita];
               } 
  // look for any special vdw pairs
               numeral(type[i],pa,3);
               numeral(type[j],pb,3);
      
              if (type[i] < type[j])
              {
                strcpy(pt,pa);
                strcat(pt,pb);
              } else
              {
                strcpy(pt,pb);
                strcat(pt,pa);
              }

// look for vdw pairs
                for (k= 0; k < vdwpr_k.nvdwpr ; k++)
                {
                   if (strcmp(vdwpr_k.kv[k],pt) == 0)
                   {
                      nonbond.vrad[ita][itb] = vdwpr_k.radius[k];
                      nonbond.vrad[itb][ita] = vdwpr_k.radius[k];
                      nonbond.veps[ita][itb] = vdwpr_k.eps[k];
                      nonbond.veps[itb][ita] = vdwpr_k.eps[k];
                      if (ianum != 1 && ibnum != 1 )
                      {
                         nonbond.veps[ita][itb] /= units.dielec;
                         nonbond.veps[itb][ita] /= units.dielec;
                      }
                      if (ita == 1 && itb == 5)
                      {
                            nonbond.vrad14[ita][itb] = nonbond.vrad[ita][itb];
                            nonbond.vrad14[itb][ita] = nonbond.vrad[itb][ita];
                            nonbond.veps14[itb][ita] = nonbond.veps[itb][ita];
                            nonbond.veps14[ita][itb] = nonbond.veps[ita][itb];
                      }
                      break;
                   }
                }

             nonbond.iNBtype[ita][itb] = it;
             nonbond.iNBtype[itb][ita] = it;
             
             nonbond.npair++;
           }
       }
   }
   for (i=1; i <= natom ; i++)
   {
       ita = type[i];
       if (vdw1.rad[ita] == 0)
       {
            rii = vdw1.a[ita]*pow(vdw1.alpha[ita],0.25);
            vdw1.rad[ita] = 0.5*rii;
       }
   }
}

