#define EXTERN extern
#include "pcwin.h"
#include "pcmod.h"
#include "atom_k.h"

double dihdrl(int i1,int i2,int i3,int i4);
void type_mmx(void);
void hdel(int);
void set_atomdata(int,int,int,int);
int get_field(void);
void set_atomtypes(int); 
void type(void);       
void quick_type(void);

void set_atomtypes(int newtype)
{
    int i,badtype;
    char hatext[256];
    
    if (newtype == MMX)
    {
        for (i=1; i <= natom; i++)
        {
            atom.type[i] = atom.mmx_type[i];
            if (atom.type[i] < atom_k.natomtype)
              atom.tclass[i] = atom_k.tclass[atom.type[i]];
            else
              atom.tclass[i] = atom.type[i];
        }   
    } else if (newtype == GAFF)
    {
        for (i=1; i <= natom; i++)
        {
            atom.type[i] = atom.gaff_type[i];
            if (atom.type[i] < atom_k.natomtype)
              atom.tclass[i] = atom_k.tclass[atom.type[i]];
            else
              atom.tclass[i] = atom.type[i];
        }
    } else if (newtype == MMFF94)
    {
        for (i=1; i <= natom; i++)
        {
            atom.type[i] = atom.mmff_type[i];
            if (atom.type[i] < atom_k.natomtype)
              atom.tclass[i] = atom_k.tclass[atom.type[i]];
            else
              atom.tclass[i] = atom.type[i];
        }
    } 
    badtype = FALSE;
    
    for (i=1; i <= natom; i++)
         if (atom.type[i] == 0)
             badtype = TRUE;
     if (badtype == TRUE)
     {
           sprintf(hatext,"Error setting atom types in FField: %d \n",newtype);
           message_alert(hatext,"ERROR");
     }
}
void type()
{
  int i, j, jjbo, icount,ncarbon, nhyd,field;

    icount = 0;
    ncarbon = 0;
    nhyd = 0;
    field = get_field();
    // quick check for hydrogens on carbon
    for (i=1; i <= natom; i++)
      {
	if (atom.atomnum[i] == 1 && atom.mmx_type[i] == 5 ) // got C-H
	  {
	    type_mmx();
	    set_atomtypes(field);
	    return;
	  }
      }
    // no CH - check more closely
    for (i=1; i <= natom; i++)
    {
        if (atom.atomnum[i] == 6)
        {
            jjbo = 0;
            for (j=0; j < MAXIAT; j++)
            {
                if (atom.iat[i][j] != 0 && atom.bo[i][j] != 9)
                   jjbo += atom.bo[i][j];
                if (atom.atomnum[atom.iat[i][j]] == 1)
                   nhyd++;
            }
            if ( jjbo < 3 && atom.mmx_type[i] == 40 && icount > 0  )
            {
                quick_type();
                set_atomtypes(field);
                return;
            } else if (jjbo < 4 && atom.mmx_type[i] != 40 && atom.mmx_type[i] != 4 && icount > 0)
            {
                quick_type();
                set_atomtypes(field);
                return;
            }else
            {
                icount++;
                ncarbon++;
            }
            if (icount > 15)
               break; 
        } else if (atom.atomnum[i] == 8)
        {
            jjbo = 0;
            for (j=0; j < MAXIAT; j++)
            {
                if (atom.iat[i][j] != 0 && atom.bo[i][j] != 9)
                   jjbo += atom.bo[i][j];
            }
            if (jjbo < 1 && icount > 0)
            {
                quick_type();
                set_atomtypes(field);
                return;
            } else if (jjbo == 1 && icount > 0 && atom.mmx_type[i] == 6)
            {
                quick_type();
                set_atomtypes(field);
                return;                
            }else
                icount++;
            if (icount > 15)
               break; 
        }
        
    }
    if (ncarbon > 0 && nhyd == 0)  // pathological cases where there are no carbons, or no carbons that normally have hydrogens 
    {
        type_mmx();
        set_atomtypes(field);
        return;
    }
// normal typing routines   
    if (field == MMX)
    {
        type_mmx();
        set_atomtypes(MMX);
    } else if (field == GAFF)
    {
        for (i=1; i <= natom; i++)
        {
            if (atom.mmx_type[i] == 20)
            {
                hdel(1);
                break;
            }
        }
        type_mmx();
        set_atomtypes(GAFF);
    } else if (field == MMFF94)
    {
        for (i=1; i <= natom; i++)
        {
            if (atom.mmx_type[i] == 20)
            {
                hdel(1);
                break;
            }
        }
        type_mmx();
        set_atomtypes(MMFF94);
    }
    else
       message_alert("No typing rules for this force field are implemented!","ERROR");  
} 
/* ------------------------ */
double dihdrl(int i1,int i2,int i3,int i4)
{
        float           a1, a2, ang, b1, b2, c1, c2, siign,
                        x1, x3, x4, y1, y3, y4, z1, z3, z4;
        double          dihdrl_v,r1,r2;

        /****the arguments of this function are the atom indices of atoms a-b-c-d
         *     which determine dihedral angle dihdrl.
         *     dihdrl returns the degree value of the dihedral('torsional')
         *     angle defined by atoms i1,i2,i3, &i4.
         */
         
        dihdrl_v = 0.0;
        if (i1 <= 0 || i2 <= 0 || i3 <= 0 || i4 <= 0)
           return (dihdrl_v);
        if (i1 > natom || i2 > natom || i3 > natom || i4 > natom)
           return (dihdrl_v);        
           
        a1 = atom.x[i2];
        b1 = atom.y[i2];
        c1 = atom.z[i2];
        x1 = atom.x[i1] - a1;
        y1 = atom.y[i1] - b1;
        z1 = atom.z[i1] - c1;
        x3 = atom.x[i3] - a1;
        y3 = atom.y[i3] - b1;
        z3 = atom.z[i3] - c1;
        x4 = atom.x[i4] - a1;
        y4 = atom.y[i4] - b1;
        z4 = atom.z[i4] - c1;
        a1 = y1 * z3 - y3 * z1;
        b1 = x3 * z1 - x1 * z3;
        c1 = x1 * y3 - x3 * y1;
        a2 = y4 * z3 - y3 * z4;
        b2 = x3 * z4 - x4 * z3;
        c2 = x4 * y3 - x3 * y4;
        r1 = sqrt(a1 * a1 + b1 * b1 + c1 * c1);
        r2 = sqrt(a2 * a2 + b2 * b2 + c2 * c2);
        if (r1 != 0.0 && r2 != 0.0)
           ang = (a1 * a2 + b1 * b2 + c1 * c2) /(r1*r2);
        else
           ang = 0.0;

        if (fabs(ang) > 1)
                ang = fabs(ang) / ang;
        ang = acos(ang);
        dihdrl_v = 57.29578F * ang;
        siign = x1 * a2 + y1 * b2 + z1 * c2;
        if (siign < 0.0F)
                dihdrl_v = -dihdrl_v;
        return (dihdrl_v);
}                               /* end of function */
/* --------------------------------- */
/* ------------------------------------------------------------- */       
void quick_type()
{
   int    i, iatik, iatkkk, iatomt, iatyp, ibt[MAXIAT], ibtk,
           ibut, idc, itia, j, ja, jji, k, k1, k2, katm, ki,
           kk, knt, l, lat, latm, lnt, k3, matm, jatm;
   int nox;
   int    mmfftype = 0;   

        /* *** loop over all atoms and determine type *** */
        for (i = 1; i <= natom; i++)
        {
            if (atom.mmx_type[i] != 0)
            {
                iatyp = atom.atomnum[i];
                iatomt = atom.mmx_type[i];
                /* *** load all attached bond types into ibt() *** except for
                 * coordinated bonds */
                for (j = 0; j < MAXIAT; j++) 
                {
                        if (atom.bo[i][j] == 9 || atom.iat[i][j] == 0) 
                            ibt[j] = 0;
                        else 
                            ibt[j] = atom.bo[i][j];
                 }
                /* *** now branch on the atom type *** */
                if (iatyp == 6)
                        goto L_20;
                if (iatyp == 1)
                        goto L_100;
                if (iatyp == 8)
                        goto L_140;
                if (iatyp == 7)
                        goto L_170;
                if (iatyp == 16)
                        goto L_230;
                if (iatyp == 34)
                        goto L_250;
                if (iatyp == 15)
                        goto L_270;
                if (iatyp == 13)  // Aluminum
                        goto L_271;
                if (iatyp == 17)  // chlorine
                        goto L_272;
                /* *** all other atoms have only one type *** */
                goto L_300;

                /* *** carbon *** */
L_20:
                iatomt = 1;
                mmfftype = 1;
                for (k = 0; k < MAXIAT; k++)
                {
                    if (atom.iat[i][k] != 0 && atom.bo[i][k] >= 2 && atom.bo[i][k] != 9)
                    {
		      iatomt = 2;
		      mmfftype = 2;
                        /* *** sp carbon *** */
                        if (atom.bo[i][k] == 3)
                        {
                            iatomt = 4;
                            mmfftype = 4;
                            goto L_290;
                        }
                        /* *** sp2 carbon *** */
                        ja = atom.iat[i][k];
                        if (atom.atomnum[ja] == 8 )
                        {
                                iatomt = 3;
                                mmfftype = 3;
                                /* check for ketenes and co2 */
                                kk = k + 1;
                                for (ki = kk; ki <= 10; ki++)
                                {
                                        if (atom.bo[i][ki] == 2)
                                        {
                                                iatomt = 4;
                                                mmfftype = 4;
                                                goto L_290;
                                        }
                                }
                                goto L_290;
                        }
                        /* check for allene */
                        kk = k + 1;
                        for (ki = kk; ki <= 8; ki++)
                        {
                                if (ibt[ki] == 2)
                                {
                                        iatomt = 4;
                                        mmfftype = 4;
                                        goto L_290;
                                }
                        }
                    }
                }
                /* *** 3 and 4 -mem ring carbons *** */
                if (iatomt != 1)
                        goto L_290;
                for (k = 0; k < MAXIAT; k++) 
                {
                        if (atom.iat[i][k] != 0 && atom.bo[i][k] != 9)
                        {       
                                katm = atom.iat[i][k];
                                for (k1 = 0; k1 < MAXIAT; k1++) 
                                {
                                        if (atom.iat[katm][k1] != 0 && atom.iat[katm][k1] != i && atom.bo[katm][k1] != 9)
                                        {
                                                latm = atom.iat[katm][k1];
                                                for (k2 = 0; k2 < MAXIAT; k2++) 
                                                {
                                                        if (atom.iat[latm][k2] != 0 && atom.iat[latm][k2] != katm && atom.bo[latm][k2] != 9)
                                                        {
                                                                if (atom.iat[latm][k2] == i) 
                                                                {
                                                                        iatomt = 22;
                                                                        mmfftype = 22;
                                                                        ibut = 1;
                                                                        goto L_290;
                                                                }
                                                        }
                                                }
                                        }
                                }
                        }
                }
//  now check for four membered rings
                for (k = 0; k < MAXIAT; k++) 
                {
                        if (atom.iat[i][k] != 0 && atom.bo[i][k] != 9)
                        {       
                                katm = atom.iat[i][k];
                                for (k1 = 0; k1 < MAXIAT; k1++) 
                                {
                                        if (atom.iat[katm][k1] != 0 && atom.iat[katm][k1] != i && atom.bo[katm][k1] != 9)
                                        {
                                                latm = atom.iat[katm][k1];
                                                for (k2 = 0; k2 < MAXIAT; k2++) 
                                                {
                                                        if (atom.iat[latm][k2] != 0 && atom.iat[latm][k2] != katm && atom.iat[latm][k2] != i  && atom.bo[latm][k2] != 9)
                                                        {
                                                                matm = atom.iat[latm][k2];
                                                                for (k3 = 0; k3 < MAXIAT; k3++)
                                                                {
                                                                        if (atom.iat[matm][k3] != 0 && atom.iat[matm][k3] != latm && atom.bo[matm][k3] != 9)
                                                                        {
                                                                                if (atom.iat[matm][k3] == i) 
                                                                                {
                                                                                        iatomt = 56;
                                                                                        mmfftype = 20;
                                                                                        goto L_290;
                                                                                }
                                                                        }
                                                                }
                                                        }
                                                }
                                        }
                                }
                        }
                }
                goto L_290;

                /* *** hydrogen *** */
L_100:
                iatomt = 5;
                mmfftype = 5;
                knt = 0;
                for (k=0;k < MAXIAT; k++)
                  if (atom.iat[i][k] != 0)
                     knt++;
                if (knt == 2)
                {
                   iatomt = 70;
                   goto L_290;
                }
                for (k = 0; k < MAXIAT; k++) 
                {
                   if (atom.iat[i][k] != 0) 
                   {
                       ja = atom.mmx_type[atom.iat[i][k]];
                       goto L_120;
                   }
                }
                goto L_290;
L_120:
                if (((ja == 8 || ja == 9) || ja == 37) || ja == 55)
                {
                        iatomt = 23;
                        mmfftype = 23;
                } else if (ja == 41 || ja == 46)
                {
                        iatomt = 24;
                        mmfftype = 36;
                } else if (ja == 6 || ja == 53)
                {
                        iatomt = 21;
                        mmfftype = 21;
                        for (j = 0; j < MAXIAT; j++) {
                                if (atom.iat[atom.iat[i][k]][j] == 0)
                                        goto L_130;
                                if (atom.iat[atom.iat[i][k]][j] == i)
                                        goto L_130;
                                if (atom.mmx_type[atom.iat[atom.iat[i][k]][j]] == 3)
                                {
                                        iatomt = 24;
                                        mmfftype = 24;
                                } else if (atom.mmx_type[atom.iat[atom.iat[i][k]][j]] == 2 ||
                                           atom.mmx_type[atom.iat[atom.iat[i][k]][j]] == 40)
                                {
                                        iatomt = 28;
                                        mmfftype = 28;
                                }
                L_130:
                                ;
                        }
                }
                goto L_290;

                /* *** oxygen *** */
L_140:
                if (iatomt == 42)
                {
                    mmfftype = 32;
                    goto L_290;
                }
                iatomt = 6;
                mmfftype = 6;
                for (k = 0; k < MAXIAT; k++)
                {
                    if (atom.iat[i][k] != 0)
                    {
                        if (ibt[k] == 2)
                        {
                            iatomt = 7;
                            if (atom.atomnum[atom.iat[i][0]] == 16)
                                mmfftype = 7;
			    if (atom.atomnum[atom.iat[i][0]] == 7 && atom.input_charge[atom.iat[i][0]] == 1)
			      mmfftype = 32;
                            goto L_290;
                        } else if (ibt[k] == 3)
                        {
                            for (lat = 0; lat < MAXIAT; lat++)
                            {
                                 if (atom.iat[i][lat] != 0)
                                 {
                                     if (atom.mmx_type[atom.iat[i][lat]] == 4)
                                     {
                                         iatomt = 46;
                                         mmfftype = 49;
                                         goto L_290;
                                     }
                                 }
                            }
                        }
                    }
                }
                goto L_290;

                /* *** nitrogen *** */
L_170:
                iatomt = 8;
                mmfftype = 8;
                ibtk = 0;
                /* determine maximum number of ligands attached to nitrogen */
                jji = 0;
                for (k = 0; k < MAXIAT; k++) 
                {
                        if (atom.iat[i][k] != 0 && atom.bo[i][k] != 9)
                                jji += 1;
                }
                for (k = 0; k < MAXIAT; k++) 
                {
                        if (ibt[k] != 9) 
                        {
                             ibtk += ibt[k];
                        }
                        if (atom.iat[i][k] != 0) 
                        {
                                if (atom.mmx_type[atom.iat[i][k]] == 20 ||
                                   (atom.mmx_type[atom.iat[i][k]] >= 300) )
                                        goto L_190;
                        }
                }
                /* ammonium nitrogen - no lone pairs */
                if (ibtk == 4) 
                {
                        iatomt = 41;
                        mmfftype = 34;
			nox = 0;
			for (k=0; k < 4; k++)
			  {
			    if (atom.atomnum[ibt[k]] == 8)
			      nox++;
			  }
			if (nox == 2) // nitro nitrogen
			  mmfftype = 45;
                        goto L_290;
                } else if (ibtk < 4 && atom.input_charge[i] == 1)
		  {
		    iatomt = 41;
		    mmfftype = 34;
		    goto L_290;
		  }
                /* got a lone pair can not be ammonium */
L_190:
                for (k = 0; k < MAXIAT; k++) 
                {
                        /* azo, imine */
                        if (ibt[k] == 2) 
                        {
                                iatomt = 37;
                                mmfftype = 9;
                                goto L_290;
                                /* nitrile */
                        } else if (ibt[k] == 3) 
                        {
                                iatomt = 10;
                                mmfftype = 42;
                                goto L_290;
                        }
                }
                /* count the number of double bonds to an adjacent Sulfur or
                 * Selenium */
                idc = 0;
                for (k = 0; k < MAXIAT; k++) 
                {
                     iatik = atom.iat[i][k];
                     if (iatik != 0 && atom.bo[i][k] != 9) 
                     {
                         if (atom.atomnum[iatik] == 16) 
                         {
                            for (kk = 0; kk < MAXIAT; kk++) 
                            {
                                 iatkkk = atom.iat[iatik][kk];
                                 if (iatkkk != 0 && atom.bo[iatik][kk] != 9) 
                                 {
                                     if (atom.bo[iatik][kk] > 1)
                                          idc += 1;
                                 }
                            }
                         }
                     }
                }
                if (idc == 1) 
                {
                    iatomt = 9;
                    mmfftype = 9;
                } else if (idc == 2) 
                {
                    iatomt = 8;
                    mmfftype = 9;
                }
                /* *** loop over all adjacent atoms to find enamines, amides */
                for (k = 0; k < MAXIAT; k++) 
                {
                   if (atom.iat[i][k] != 0 && atom.bo[i][k] != 9 )
                   {
                        itia = atom.mmx_type[atom.iat[i][k]];
                        if ( ((itia >= 300) || itia == 29 || itia == 30 || itia == 40  || itia == 26) && jji <= 3)
                        {
                            iatomt = 9;
                            mmfftype = 9;
                            goto L_290;
                        }
                        if (atom.atomnum[atom.iat[i][k]] == 6) // adjacent carbon
                        {
                            jatm = atom.iat[i][k];
                            for (j=0; j < MAXIAT; j++)
                            {
                                if (atom.iat[jatm][j] != 0 && atom.bo[jatm][j] >= 2)
                                {
                                    iatomt = 9;
                                    mmfftype = 10;
                                    goto L_290;
                                }
                            }
                        }
                        if (atom.atomnum[atom.iat[i][k]] == 16) // adjacent sulfur
                        {
                            jatm = atom.iat[i][k];
                            for (j=0; j < MAXIAT; j++)
                            {
                                if (atom.iat[jatm][j] != 0 && atom.bo[jatm][j] >= 2)
                                {
                                    iatomt = 9;
                                    mmfftype = 9;
                                    goto L_290;
                                }
                            }
                        }
                        if (atom.atomnum[atom.iat[i][k]] == 7) // adjacent nitrogen
                        {
                            jatm = atom.iat[i][k];
                            for (j=0; j < MAXIAT; j++)
                            {
                                if (atom.iat[jatm][j] != 0 && atom.bo[jatm][j] >= 2)
                                {
                                    iatomt = 9;
                                    mmfftype = 9;
                                    goto L_290;
                                } else if (atom.iat[jatm][j] != 0 && atom.iat[jatm][j] != i)
                                {
                                    katm = atom.iat[jatm][j];
                                    if (atom.atomnum[katm] == 6)
                                    {
                                        for (l=0; l < MAXIAT; l++)
                                        {
                                            if (atom.iat[katm][l] != 0 && atom.iat[katm][l] != jatm && atom.bo[katm][l] == 2)
                                            {
                                                iatomt = 9;
                                                mmfftype = 9;
                                                goto L_290;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                   }
                }
                goto L_290;

                /* *** sulfur *** */
L_230:
                if (iatomt == 16)
                {
                        mmfftype = 16;
                        goto L_290;
                }
                iatomt = 15;
                mmfftype = 15;
                knt = 0;
                lnt = 0;
                for (k = 0; k < MAXIAT; k++) {
                        if (atom.iat[i][k] != 0) {
                                if (ibt[k] > 0 && atom.mmx_type[atom.iat[i][k]] != 20)
                                {
                                        lnt += 1;
                                }
                                if (ibt[k] > 1)
                                {
                                        knt += 1;
                                        if (atom.mmx_type[atom.iat[i][k]] <= 4)
                                        {
                                                iatomt = 38;
                                                mmfftype = 16;
                                                goto L_290;
                                        }
                                }
                        }
                }
                if (knt == 1)
                {
                    iatomt = 17;
                    mmfftype = 17;
                    if (lnt == 1)
                        iatomt = 38;
                }
                if (knt >= 2)
                {
                    iatomt = 18;
                    mmfftype = 18;
                }
                goto L_290;
                /* ** selenium */
L_250:
                iatomt = 34;
                for (k = 0; k < MAXIAT; k++) {
                        if (atom.iat[i][k] != 0) {
                                if (ibt[k] == 2) {
                                        iatomt = 39;
                                        goto L_290;
                                }
                        }
                }
                goto L_290;
                /* ** phosphorous */
L_270:
                iatomt = 25;
                mmfftype = 25;
                knt = 0;
                for (j = 0; j < MAXIAT; j++) {
                        if (atom.iat[i][j] != 0)
                                knt += 1;
                }
                if (knt >= 5)
                        iatomt = 47;
                goto L_290;
                /*  aluminum  */
L_271:
               iatomt = 44;
               knt = 0;
                for (j = 0; j < MAXIAT; j++)
                {
                   if (atom.iat[i][j] != 0)
                      knt += 1;
                }
                if (knt == 4)
                   iatomt = 58;
                goto L_290;
           /*  Chlorine  */
L_272:
               iatomt = 12;
               mmfftype = 12;
               knt = 0;
                for (j = 0; j < MAXIAT; j++)
                {
                   if (atom.iat[i][j] != 0)
                      knt += 1;
                }
                if (knt == 2)
                   iatomt = 74;
                goto L_290;
                /* *** typing is now done, load value and return *** */
L_290:
//                if (field == MMX && iatomt < 300 )
                set_atomdata(i,iatomt,0,mmfftype);
L_300:
                continue;
            }
        }
        return;
}
                 
