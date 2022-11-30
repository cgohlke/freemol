#define EXTERN extern

#include "pcwin.h"

int get_hybrid(int);
int check_ring1(int);
float xlogp(int natom,int *atomnum,int **iat,int **bo,long int *flags,int *type);
// =========================
float xlogp(int natom,int *atomnum,int **iat,int **bo,long int *flags,int *type)
{
   int i,j,k,ia;
   int iz,nh,nx,npi,nc,itype,dbond,nox,ndouble,nf;
   int npi_c,npi_x;
   int nhydrophobic, namino_acid, nhbond, n13FX, n13XX;
   int atype[100];
   long int pi_mask;
   float total,result;
//   char astring[60];
   float datype[] = {
    0.000, 0.484, 0.168,-0.181, 0.358, 0.009,-0.344,-0.439, 0.051,-0.138,
   -0.417,-0.454,-0.378, 0.223,-0.598,-0.396,-0.699,-0.362, 0.395, 0.236,
   -0.166, 1.726, 0.098,-0.108, 1.637, 1.774, 0.281, 0.142, 0.715, 0.302,
   -0.064, 0.079, 0.200, 0.869, 0.316, 0.054, 0.347, 0.046,-0.399,-0.029,
   -0.330, 0.397, 0.068, 0.327,-2.057, 0.218,-0.582,-0.449,-0.774, 0.040,
   -0.381, 0.443,-0.117,-2.052,-1.716, 0.321,-0.921,-0.704, 0.119, 1.192,
    0.434, 0.587, 0.668,-0.791,-0.212, 0.016, 0.752, 1.071, 0.964,-1.817,
   -1.214,-0.778, 0.493, 1.010, 1.187, 1.489,-0.802,-0.256, 1.626, 0.077,
    0.264, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
    0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 }; 

   pi_mask = (1L << PI_MASK);
   for (i=0;i < 100; i++)
      atype[i] = 0;
//
//  loop through atoms and classify
    for (i=1; i <= natom; i++)
    {
        if (atomnum[i] == 6) // carbon
        {
            iz = get_hybrid(i);
            nh = nx = npi = nc = 0;
            for (j=0; j < MAXIAT; j++)
            {
                if (iat[i][j] != 0)
                {
                    if (atomnum[iat[i][j]] == 1) nh++;
                    if (flags[iat[i][j]] & pi_mask) npi++;
                    itype = atomnum[iat[i][j]];
                    if (itype == 6) nc++;
                    if (itype == 7 || itype == 8 || itype == 15 || itype == 16) nx++;
                    itype = type[iat[i][j]];
                    if (itype == 11 || itype == 12 || itype == 13 || itype == 14) nx++;
                }
            }
            if (iz == 1) // sp3
            {
                if (nh == 3 && nc == 1 && npi == 0)
                  atype[1]++;
                else if (nh == 3 && nc == 1 && npi > 0)
                  atype[2]++;
                else if (nh == 3 && nx == 1)
                  atype[3]++;
                else if (nh == 2 && nc == 2 && npi == 0)
                  atype[4]++;
                else if (nh == 2 && nc == 2 && npi > 0)
                  atype[5]++;
                else if (nh == 2 && nc == 1 && nx == 1)
                  atype[6]++;
                else if (nh == 2 && nx == 2)
                  atype[7]++;
                else if (nh == 1 && nc == 3 && npi == 0)
                  atype[8]++;
                else if (nh == 1 && nc == 3 && npi > 0)
                  atype[9]++;
                else if (nh == 1 && nc == 2 && nx == 1)
                  atype[10]++;
                else if (nh == 1 && nc == 1 && nx == 2)
                  atype[11]++;
                else if (nh == 1 && nx == 3)
                  atype[11]++;
                else if (nh == 0 && nc == 4 && npi == 0)
                  atype[12]++;
                else if (nh == 0 && nc == 4 && npi > 0)
                  atype[13]++;
                else if (nh == 0 && nc == 3 && nx == 1)
                  atype[14]++;
                else if (nh == 0 && nc == 2 && nx == 2)
                  atype[15]++;
                else if (nh == 0 && nc == 1 && nx == 3)
                  atype[16]++;
                else if (nh == 0 && nc == 0 && nx == 4)
                  atype[17]++;   
            } else if (iz == 2) // sp2
            {
                dbond = 0;
                // find double bond atom
                for (j=0; j < 4; j++)
                {
                    if (iat[i][j] != 0 && bo[i][j] == 2)
                    {
                        dbond = iat[i][j];
                        break;
                    }
                }
                if (dbond == 0) goto L_10; // bad atoms
                if (nh == 2 && nc == 1)
                {
                  atype[18]++;
                  goto L_10;
                }
                // if not a pi atom
                if (!(flags[i] & pi_mask))  // non pi atoms
                {
                    if (atomnum[dbond] == 6)  // c=c
                    {
                        if (nh == 1 && nc == 2)
                        {
                           atype[19]++;
                           goto L_10;
                        } else if (nh == 1 && nx == 1)
                        {
                            atype[20]++;
                            goto L_10;
                        } else if (nh == 0 && nc == 3)
                        {
                            atype[22]++;
                            goto L_10;
                        } else if (nh == 0 && nx >= 1)
                        {
                            atype[23]++;
                            goto L_10;
                        }
                    }else   // c=x
                    {
                        if (nh == 1)
                        {
                            atype[21]++;
                            goto L_10;
                        } else if (nc >= 1)
                        {
                            atype[24]++;
                            goto L_10;
                        } else if (nx == 3)
                        {
                            atype[25]++;
                            goto L_10;
                        }
                    }
                } else   // pi conjugated systems
                {
                    nh = nc = nx = npi = 0;
                    npi_c = 0;
                    npi_x = 0;
                    for (j = 0; j < 3; j++)
                    {
                        if (atomnum[iat[i][j]] == 1)
                          nh++;
                        if (atomnum[iat[i][j]] == 6)
                        {
                            nc++;
                            if (flags[iat[i][j]] & pi_mask)
                            {
                                npi++;
                                npi_c++;
                            }
                        }
                        if (atomnum[iat[i][j]] >= 6)
                        {
                          nx++;
                          if (flags[iat[i][j]] & pi_mask)
                          {
                              npi++;
                              npi_x++;
                          }
                        }
                    }
                    if (nh == 1 && nc == 2)
                      atype[26]++;
                    else if (nh == 1 && nc ==1 && nx == 1)
                      atype[27]++;
                    else if (nh == 1 && nx == 2)
                      atype[28]++;
                    else if (nc == 3 && npi == 2)
                      atype[29]++;
                    else if (nc == 2 && npi_x == 2 && nx == 1)
                      atype[30]++;
                    else if (nc == 2 && npi_x == 1 && nx == 1 && npi_x == 1)
                      atype[31]++;
                    else if (nc == 1 && nx == 2 && npi_x == 1 && npi_c == 1)
                      atype[32]++;
                    else if (nc == 1 && nx == 2 && npi_x == 2)
                      atype[33]++;
                    else if (npi == 3)
                      atype[34]++;
                }
            } else if (iz == 3)
            {
                for (j=0; j < 4; j++)
                {
                    if (iat[i][j] != 0 && bo[i][j] == 3)
                    {
                        if (atomnum[iat[i][j]] == 7)
                        {
                            atype[77]++;
                            goto L_10;
                        }
                    }
                }
                if (nh == 1)
                  atype[35]++;
                else if (nc == 1)
                  atype[36]++;
            }
        } else if (atomnum[i] == 1) // hydrogen
        {
            atype[37]++;
        } else if (atomnum[i] == 8) // oxygen
        {
            nh = npi = nx = nc = 0;
            for (j=0; j < 4; j++)
            {
                if (atomnum[iat[i][j]] == 1)
                  nh++;
                if (atomnum[iat[i][j]] == 6)
                {
                    nc++;
                    if (flags[iat[i][j]] & pi_mask) npi++;
                }
                if (atomnum[iat[i][j]] > 6)
                {
                    nx++;
                    if (flags[iat[i][j]] & pi_mask) npi++;
                }
            }
            if (type[i] == 6)
            {
                if (nh == 1 && nc == 1 && npi == 0)
                  atype[38]++;
                else if (nh == 1 && nc == 1 && npi == 1)
                  atype[39]++;
                else if (nh == 1 && nx == 1)
                  atype[40]++;
                else if (nc == 2 && npi < 2)
                  atype[41]++;
                else if ( (nc == 1 && nx == 1) || nx == 2)
                  atype[42]++;
                else if (nc == 2 && npi == 2)
                  atype[43]++;
            } else if (type[i] == 7)
            {
                if (nc == 1)
                  atype[44]++;
                else if (nx == 1)
                  atype[45]++;
            }   
        } else if (atomnum[i] == 7) // nitrogen
        {
            nh = nx = nc = npi = 0;
            for (j=0; j < MAXIAT; j++)
            {
                if (iat[i][j] != 0)
                {
                    if (atomnum[iat[i][j]] == 1) nh++;
                    if (atomnum[iat[i][j]] == 6) nc++;
                    if (atomnum[iat[i][j]] > 6) nx++;
                    if (flags[iat[i][j]] & pi_mask) npi++;
                }
            }
            if (type[i] == 8)
            {
                if (nh == 2 && nc == 1 && npi == 0) atype[46]++;
                else if (nh == 2 && nc == 1 && npi > 1)  atype[47]++;
                else if (nh == 2 && nx == 1) atype[48]++;
                else if (nh == 1 && nc == 2) atype[49]++;
                else if (nh == 1 && nc == 1 && nx == 1) atype[50]++;
                else if (nh == 1 && nx == 2) atype[50]++;
                else if (nh == 0 && nc == 3) atype[51]++;
                else if (nh == 0 && nx >= 1) atype[52]++;
            } else if (type[i] == 9)
            {
                for (j=0; j < 4; j++)  
                {
                    if (iat[i][j] != 0 )
                    {
                        if (bo[i][j] == 2)  // imine type
                        {
                            if (atomnum[iat[i][j]] == 6)
                            {
                                if (nx == 0)
                                   atype[53]++;
                                else
                                   atype[54]++;
                                goto L_10;
                            } else if (atomnum[iat[i][j]] > 6)
                            {
                                if (nc == 1)
                                  atype[55]++;
                                else
                                  atype[56]++;
                                goto L_10;
                            }
                        }
                    }
                }
                for (j=0; j < 4; j++)  // amides
                {
                    if (iat[i][j] != 0 && atomnum[iat[i][j]] == 6)
                    {
                        ia = iat[i][j];
                        for (k=0; k < 4; k++)
                        {
                            if (iat[ia][k] != 0 && bo[ia][k] == 2)
                            {
                                if (atomnum[iat[ia][k]] == 8) // amide
                                {
                                    if (nh == 2 && nc == 1) atype[63]++;
                                    else if (nh == 1) atype[64]++;
                                    else atype[65]++;
                                    goto L_10;
                                }
                            }
                        }
                    }
                }
                // picked up amides and imines - everything else must be aromatic
                if (nh == 1 && nc == 2) atype[58]++;
                else if (nh == 1 && nx >= 1) atype[59]++;
                else if (nc == 3) atype[61]++;
                goto L_10;
            }
        } else if (atomnum[i] == 16) // sulfur
        {
            nh = 0;
            nox = 0;
            ndouble = 0;
            for (j=0; j < MAXIAT; j++)
            {
                if (iat[i][j] != 0)
                {
                    if (atomnum[iat[i][j]] == 1) nh++;
                    if (atomnum[iat[i][j]] == 8) nox++;
                    if (bo[i][j] == 2) ndouble++;
                }
            }
            if (nox == 2)
            {
               atype[71]++;
                goto L_10;
            } else if (nox == 1)
            {
              atype[70]++;
                goto L_10;
            }
            if (nox == 0 && ndouble == 1)
            {
                atype[69]++;
                goto L_10;
            }
            if (nh == 1)
            {
                atype[66]++;
                goto L_10;
            } else
            {
                if (check_ring1(i))
                {
                    if ((flags[iat[i][0]] & pi_mask) && (flags[iat[i][1]] & pi_mask))
                    {
                        atype[68]++;
                        goto L_10;
                    }
                } else
                {
                    atype[67]++;
                    goto L_10;
                }
            }
        } else if (atomnum[i] == 15) // phosphorus
        {
            atype[76]++;
        } else if (atomnum[i] == 11) // fluorine
        {
            atype[72]++;
        } else if (atomnum[i] == 17) // chlorine
        {
            atype[73]++;
        } else if (atomnum[i] == 35) // bromine
        {
            atype[74]++;
        } else if (atomnum[i] == 53) // iodine
        {
            atype[75]++;
        } else
        {
//            sprintf(astring,"Atom %d type %d not supported\0",i,atom[i].type);
//            message_alert(astring,"Error in LogP");
        }
L_10:
        continue;
    }
// done assigning types = compute logP
    total = 0.0;
    for (i=0; i < 100; i++)
       total += atype[i]*datype[i];

// now do correction factors
    nhydrophobic = namino_acid = nhbond = n13FX = n13XX = 0;

    nhydrophobic = atype[1] + atype[2] + atype[4] + atype[5] + atype[8] + atype[9] + atype[12] + atype[13] + atype[18] +
                   atype[19] + atype[22];

// 1,3 halogens
         for (i=1; i <= natom; i++)
         {
             nf = 0;
             nx = 0;
             if (atomnum[i] == 6)
             {
                 for (j=0; j < MAXIAT; j++)
                 {
                     if (atomnum[iat[i][j]] == 9)
                       nf++;
                     if (atomnum[iat[i][j]] == 17)
                       nx++;
                     if (atomnum[iat[i][j]] == 35)
                       nx++;
                     if (atomnum[iat[i][j]] == 53)
                       nx++;
                 }
                 if ((nf >= 1 && nx >= 1) || nf >= 2)
                    n13FX++;
                 if (nx >=2 && nf == 0)
                    n13XX++;
             }
         }
// amino acids
// hydrogen bonds

    
    result = total + nhydrophobic*0.19F + n13FX*0.08 + n13XX*(-0.26F);
    return result;
}
// ==================
