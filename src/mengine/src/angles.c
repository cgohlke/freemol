#define EXTERN extern

#include "pcwin.h"
#include "pcmod.h"
#include "angles.h"

void get_angles()
{
    int i,j,icount;
    int ia, ib, ic, id, ie;
    int jj, icoord;
    int itemp[5];

    ia = ib = 0;
    angles.nang = 0;
    angles.nopb = 0;

    for (i=0; i < MAXANG; i++)
    {
        angles.i13[i][0] = 0;
        angles.i13[i][1] = 0;
        angles.i13[i][2] = 0;
        angles.i13[i][3] = 0;
    }     
//
    for (i = 1; i <= natom; i++)
    {
        jj = 0;
        icoord = 0;
        for (j=0; j < MAXIAT; j++)
        {
            if (atom.iat[i][j] != 0 && atom.bo[i][j] != 9)
              jj++;
            if (atom.bo[i][j] == 9)
              icoord++;
        }

        if (jj == 2 && atom.type[i] < 300)
        {
            icount = 0;
            for (j=0; j < MAXIAT; j++)
            {
                if (atom.iat[i][j] != 0 && atom.bo[i][j] != 9)
                {
                    icount++;
                    if (icount == 1)
                       ia = atom.iat[i][j];
                    else if (icount == 2)
                    {
                        ib = atom.iat[i][j];
                        break;
                    }
                }
            }
            angles.i13[angles.nang][0] = ia;
            angles.i13[angles.nang][1] = i;
            angles.i13[angles.nang][2] = ib;
            angles.i13[angles.nang][3] = 0;
            angles.nang++;            
        } else if (jj == 3 && atom.type[i] < 300)
        {
            if (icoord == 0)
            {
             ia = atom.iat[i][0];
             ib = atom.iat[i][1];
             ic = atom.iat[i][2];
            } else
            {
              icoord = 0;
              for(j=0; j < MAXIAT; j++)
              {
                 if (atom.iat[i][j] !=0 && atom.bo[i][j] != 9)
                 {
                     itemp[icoord] = atom.iat[i][j];
                     icoord++;
                 }
              }
              ia = itemp[0];
              ib = itemp[1];
              ic = itemp[2];
            }
            angles.i13[angles.nang][0] = ia;
            angles.i13[angles.nang][1] = i;
            angles.i13[angles.nang][2] = ib;
            angles.i13[angles.nang][3] = ic;
            angles.nang++;            
            angles.i13[angles.nang][0] = ia;
            angles.i13[angles.nang][1] = i;
            angles.i13[angles.nang][2] = ic;
            angles.i13[angles.nang][3] = ib;
            angles.nang++;            
            angles.i13[angles.nang][0] = ib;
            angles.i13[angles.nang][1] = i;
            angles.i13[angles.nang][2] = ic;
            angles.i13[angles.nang][3] = ia;
            angles.nang++;            
        } else if (jj == 4 && atom.type[i] < 300) 
        {
            if (icoord == 0)
            {
              ia = atom.iat[i][0];
              ib = atom.iat[i][1];
              ic = atom.iat[i][2];
              id = atom.iat[i][3];
            } else
            {
                icoord = 0;
                for (j=0; j < MAXIAT; j++)
                {
                    if (atom.iat[i][j] != 0 && atom.bo[i][j] != 9)
                    {
                        itemp[icoord] = atom.iat[i][j];
                        icoord++;
                    }
                }
                ia = itemp[0];
                ib = itemp[1];
                ic = itemp[2];
                id = itemp[3];
            }
            angles.i13[angles.nang][0] = ia;
            angles.i13[angles.nang][1] = i;
            angles.i13[angles.nang][2] = ib;
            angles.i13[angles.nang][3] = 0;
            angles.nang++;            
            angles.i13[angles.nang][0] = ia;
            angles.i13[angles.nang][1] = i;
            angles.i13[angles.nang][2] = ic;
            angles.i13[angles.nang][3] = 0;
            angles.nang++;            
            angles.i13[angles.nang][0] = ia;
            angles.i13[angles.nang][1] = i;
            angles.i13[angles.nang][2] = id;
            angles.i13[angles.nang][3] = 0;
            angles.nang++;            
            angles.i13[angles.nang][0] = ib;
            angles.i13[angles.nang][1] = i;
            angles.i13[angles.nang][2] = ic;
            angles.i13[angles.nang][3] = 0;
            angles.nang++;            
            angles.i13[angles.nang][0] = ib;
            angles.i13[angles.nang][1] = i;
            angles.i13[angles.nang][2] = id;
            angles.i13[angles.nang][3] = 0;
            angles.nang++;            
            angles.i13[angles.nang][0] = ic;
            angles.i13[angles.nang][1] = i;
            angles.i13[angles.nang][2] = id;
            angles.i13[angles.nang][3] = 0;
            angles.nang++;                        
        } else if (jj == 5 && atom.type[i] < 300)
        {
            if (icoord == 0)
            {
              ia = atom.iat[i][0];
              ib = atom.iat[i][1];
              ic = atom.iat[i][2];
              id = atom.iat[i][3];
              ie = atom.iat[i][4];
            } else
            {
                icoord = 0;
                for (j=0; j < MAXIAT; j++)
                {
                    if (atom.iat[i][j] != 0 && atom.bo[i][j] != 9)
                    {
                        itemp[icoord] = atom.iat[i][j];
                        icoord++;
                    }
                }
                ia = itemp[0];
                ib = itemp[1];
                ic = itemp[2];
                id = itemp[3];
                ie = itemp[4];
            }
            angles.i13[angles.nang][0] = ia;
            angles.i13[angles.nang][1] = i;
            angles.i13[angles.nang][2] = ib;
            angles.i13[angles.nang][3] = 0;
            angles.nang++;            
            angles.i13[angles.nang][0] = ia;
            angles.i13[angles.nang][1] = i;
            angles.i13[angles.nang][2] = ic;
            angles.i13[angles.nang][3] = 0;
            angles.nang++;            
            angles.i13[angles.nang][0] = ia;
            angles.i13[angles.nang][1] = i;
            angles.i13[angles.nang][2] = id;
            angles.i13[angles.nang][3] = 0;
            angles.nang++;            
            angles.i13[angles.nang][0] = ia;
            angles.i13[angles.nang][1] = i;
            angles.i13[angles.nang][2] = ie;
            angles.i13[angles.nang][3] = 0;
            angles.nang++;            
            angles.i13[angles.nang][0] = ib;
            angles.i13[angles.nang][1] = i;
            angles.i13[angles.nang][2] = ic;
            angles.i13[angles.nang][3] = 0;
            angles.nang++;            
            angles.i13[angles.nang][0] = ib;
            angles.i13[angles.nang][1] = i;
            angles.i13[angles.nang][2] = id;
            angles.i13[angles.nang][3] = 0;
            angles.nang++;            
            angles.i13[angles.nang][0] = ib;
            angles.i13[angles.nang][1] = i;
            angles.i13[angles.nang][2] = ie;
            angles.i13[angles.nang][3] = 0;
            angles.nang++;            
            angles.i13[angles.nang][0] = ic;
            angles.i13[angles.nang][1] = i;
            angles.i13[angles.nang][2] = id;
            angles.i13[angles.nang][3] = 0;
            angles.nang++;                        
            angles.i13[angles.nang][0] = ic;
            angles.i13[angles.nang][1] = i;
            angles.i13[angles.nang][2] = ie;
            angles.i13[angles.nang][3] = 0;
            angles.nang++;                        
            angles.i13[angles.nang][0] = id;
            angles.i13[angles.nang][1] = i;
            angles.i13[angles.nang][2] = ie;
            angles.i13[angles.nang][3] = 0;
            angles.nang++;                        
        }

        if (angles.nang > MAXANG)
        {
            message_alert("Too many angles. PCMODEL will now exit","Angle Setup");
            fprintf(pcmlogfile,"Too many angles. PCMODEL will now exit\n");
            exit(0);
        }
    }
       
    for (i=0; i < angles.nang; i++)
       angles.angin[i] = FALSE;
}

