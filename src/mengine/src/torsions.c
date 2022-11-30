#define EXTERN extern

#include "pcwin.h"
#include "torsions.h"

int isbond(int, int);
void message_alert(char *, char *);
int is_linear(int,int *);
int is_allene(int,int,int,int,int *);
void max_torsions(int natom,int *type,int **iat,int **bo);
void get_torsions(int natom,int *type,int **iat,int **bo);

struct t_allene {
    int nallene, ntor[10];
     } allene;

int is_allene(int i, int j, int k, int l,int *type)
{
    if (type[i] == 2 && type[j] == 4 && type[k] == 2)
       return FALSE;
    if (type[j] == 2 && type[k] == 4 && type[l] == 2)
       return FALSE;
    return TRUE;
}
// ==================================================
int is_linear(int ia,int *type)
{
    if (type[ia] == 4)
       return TRUE;
    if (type[ia] == 53)
       return TRUE;
    return FALSE;
}       
//   ==================================         
void get_torsions(int natom,int *type,int **iat,int **bo)
{
    int i, j, k, l,m;
    int ia1, ia2, nRc;
    int katm, latm;
    
    l = k = m = 0;
    torsions.ntor = 0;
    allene.nallene = 0;
    for (i=1; i <=natom; i++)
    {
        ia1 = i;
        for(j=0; j< MAXIAT; j++)
        {
            if (iat[i][j] != 0 && isbond(i,iat[i][j]) && type[iat[i][j]] < 300 )
            {
                ia2 = iat[i][j];
                if (is_linear(ia2,type) == FALSE)
                {
                  for(k=0; k<MAXIAT; k++)
                  {
                    if (iat[ia2][k] != 0 && isbond(iat[i][j],iat[iat[i][j]][k]) && type[iat[i][j]] < 300 )
                    {
                        katm = iat[iat[i][j]][k];
                    if (isbond(ia2,katm) && type[katm] < 300 )
                    {
                        if (i != katm)
                        {
                            for(l=0; l<MAXIAT; l++)
                            {
			      if (isbond(katm,iat[katm][l]) && !is_linear(katm,type) )
                                {
                                  latm = iat[katm][l];
                                  nRc = is_allene(i,iat[i][j],katm,latm,type);
                                  if ( (latm != iat[i][j]) && (latm >i) && nRc)
                                  {
                                      torsions.i14[torsions.ntor][0] = i;
                                      torsions.i14[torsions.ntor][1] = iat[i][j];
                                      torsions.i14[torsions.ntor][2] = katm;
                                      torsions.i14[torsions.ntor][3] = latm;
                                      torsions.ntor++;
                                      if (torsions.ntor > MAXTOR)
                                      {
                                          message_alert("Too many torsions. PCMODEL will now exit","Torsion Setup");
                                          fprintf(pcmlogfile,"Too many torsions in torsion setup.\n");
                                          exit(0);
                                      }
                                  }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
        }
    }
    // look for allenes
    for (i=1; i <= natom; i++)
    {
        if (type[i] == 4)
        {
            if (bo[i][0] == 2 && bo[i][1] == 2) // got an cummulene
            {
                if (type[iat[i][0]] == 2 && type[iat[i][1]] == 2) // got an allene
                {
                    katm = iat[i][0];
                    latm = iat[i][1];
                    if (iat[katm][0] == i)
                    {
                        j = iat[katm][1];
                        k = iat[katm][2];
                    }else if (iat[katm][1] == i)
                    {
                        j = iat[katm][0];
                        k = iat[katm][2];
                    }else if (iat[katm][2] == i)
                    {
                        j = iat[katm][1];
                        k = iat[katm][2];
                    }
                    
                    if (iat[latm][0] == i)
                    {
                        l = iat[latm][1];
                        m = iat[latm][2];
                    }else if (iat[katm][1] == i)
                    {
                        l = iat[latm][0];
                        m = iat[latm][2];
                    }else if (iat[katm][2] == i)
                    {
                        l = iat[latm][1];
                        m = iat[latm][2];
                    }
                    torsions.i14[torsions.ntor][0] = j;
                    torsions.i14[torsions.ntor][1] = katm;
                    torsions.i14[torsions.ntor][2] = latm;
                    torsions.i14[torsions.ntor][3] = l;
                    allene.ntor[allene.nallene] = torsions.ntor;
                    allene.nallene++;
                    torsions.ntor++;
                    
                    torsions.i14[torsions.ntor][0] = j;
                    torsions.i14[torsions.ntor][1] = katm;
                    torsions.i14[torsions.ntor][2] = latm;
                    torsions.i14[torsions.ntor][3] = m;
                    allene.ntor[allene.nallene] = torsions.ntor;
                    allene.nallene++;
                    torsions.ntor++;
                    
                    torsions.i14[torsions.ntor][0] = k;
                    torsions.i14[torsions.ntor][1] = katm;
                    torsions.i14[torsions.ntor][2] = latm;
                    torsions.i14[torsions.ntor][3] = l;
                    allene.ntor[allene.nallene] = torsions.ntor;
                    allene.nallene++;
                    torsions.ntor++;
                        
                    torsions.i14[torsions.ntor][0] = k;
                    torsions.i14[torsions.ntor][1] = katm;
                    torsions.i14[torsions.ntor][2] = latm;
                    torsions.i14[torsions.ntor][3] = m;
                    allene.ntor[allene.nallene] = torsions.ntor;
                    allene.nallene++;
                    torsions.ntor++;
                    if (torsions.ntor > MAXTOR)
                    {
                        message_alert("Too many allene torsions. PCMODEL will now exit","Torsion Setup");
                        fprintf(pcmlogfile,"Too many allene torsions in torsion setup.\n");
                    }
                    if (allene.nallene > 10)
                    {
                        message_alert("Too many allene torsions. PCMODEL will now exit","Torsion Setup");
                        fprintf(pcmlogfile,"Too many allene torsions in torsion setup.\n");
                    }
                }
            }
        }
    }                  
}
// ======================================================================
void max_torsions(int natom,int *type,int **iat,int **bo)
{
    int i, j, k, l;
    int ia1, ia2;
    int katm, latm;
    
    torsions.ntor = 0;
    for (i=1; i <=natom; i++)
    {
        ia1 = i;
        for(j=0; j< MAXIAT; j++)
        {
            if (isbond(i,iat[i][j]) )
            {
                ia2 = iat[i][j];
                for(k=0; k<MAXIAT; k++)
                {
                    if (isbond(iat[i][j],iat[iat[i][j]][k]) )
                    {
                        katm = iat[iat[i][j]][k];
                        if (i != katm)
                        {
                            for(l=0; l<MAXIAT; l++)
                            {
                                if (isbond(katm,iat[katm][l]) )
                                {
                                  latm = iat[katm][l];
                                  if ( (latm != iat[i][j]) && (latm >i) )
                                  {
                                      torsions.ntor++;
                                  }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    // find allenes
    for (i=0; i <= natom; i++)
    {
        if (type[i] == 4)
        {
            if (bo[i][0] == 2 && bo[i][1] == 2)
            {
                if (type[iat[i][0]] == 2 && type[iat[i][1]] == 2)
                {
                    torsions.ntor += 4;
                }
            }
        }
    }                   
}


