/* NOTICE: this source code file has been modified for use with FreeMOL */
#define EXTERN extern
#include "pcwin.h"
#include "pcmod.h"
#include "utility.h"

int read_mol2(FILE *);
void write_mol2(void);
void make_bond(int,int,int);
int make_atom(int,float,float,float,char *);
FILE * fopen_path ( char * , char * , char * ) ;
int FetchRecord(FILE *, char *);
void mopaco(int numbonds,int *ib1,int *ib2);
void get_molecule_memory(int);
void get_rings(int natom,long int *flags,int *atomnum,int **iat,int **bo);
void set_connectivity(void);
void assign_aromaticity(void);
void quick_type(void);

EXTERN struct t_files {
        int nfiles, append, batch, icurrent;
        int istrselect;
        }       files;
EXTERN struct ElementType { char symbol[3];
                             int atomnum;
                             float weight, extweight, covradius, vdwradius;
                              int s,p,d,f,type ;
                            } Elements[];

EXTERN char *Struct_Title;
EXTERN char *sybname[];
                    
int read_mol2(FILE *infile)
{
    int niatom,nibond,icount;
    int i,j,newatom, ia1,ia2,ibo;
    int numbonds,*ib1,*ib2;
    float xtmp,ytmp,ztmp;
    char line[180],dumm1[20],dumm2[20],dumm3[20];
    char atname[7],atype[10],atsym[3];
    
//  either start read or continue reading file
    icount = 0;
    numbonds = 0;
    ib1 = NULL;
    ib2 = NULL;
    while ( FetchRecord(infile,line))
    {
        if (strcasecmp(line,"@<TRIPOS>MOLECULE") == 0)
        {
             FetchRecord(infile,line);  // title
	     //             strncpy(Struct_Title,line,sizeof(Struct_Title));
             FetchRecord(infile,line);  // natoms and nbonds
             sscanf(line,"%d %d",&niatom,&nibond);
	     // get memory for molecule
	     get_molecule_memory(niatom);
	     // allocate space for aromatic bonds
	     numbonds = 0;
	     ib1 = ivector(0,nibond);
	     ib2 = ivector(0,nibond);
	     for (i=0; i < numbonds; i++)
	       {
		 ib1[i] = 0;
		 ib2[i] = 0;
	       }
         } else if (strcasecmp(line,"@<TRIPOS>ATOM") == 0)
         {            
             for (i=1; i <= niatom; i++)
             {
                 FetchRecord(infile,line);
		 strcpy(atname,"");
                 sscanf(line,"%d %s %s %s %s %s",&ibo,atname,dumm1,dumm2,dumm3,atype);
                 xtmp = atof(dumm1);
                 ytmp = atof(dumm2);
                 ztmp = atof(dumm3);
                 
                 for(j=0; j < 3; j++)
                 {
                     if (isdigit(atname[j]) )
                     {
                         atsym[j] = '\0';
                         break;
                     } else
                         atsym[j] = atname[j];
                 }
		 if (isupper(atsym[1])) atsym[1] = '\0';
				 
                 newatom = make_atom(0,xtmp,ytmp,ztmp,atsym);
                 if (newatom == -1)
                   return(FALSE);
		 // set formal charges
		 if (strcmp(atype,"N.4") == 0)
		   atom.formal_charge[newatom] = 1;
		 else if (strcmp(atype,"N.pl3") == 0)
		   atom.formal_charge[newatom] = 1;
		 else if (strcmp(atype,"O.co2") == 0)
		   atom.formal_charge[newatom] = -1;

             }
         } else if (strcasecmp(line,"@<TRIPOS>BOND") == 0)
         {
	   numbonds = 0;
             for (i=0; i < nibond; i++)
             {
                 FetchRecord(infile,line);
                 sscanf(line,"%d %d %d %s",&j,&ia1,&ia2,dumm1);
                 if (strcmp(dumm1,"am") == 0 || strcmp(dumm1,"AM") == 0)
                     ibo = 1;
                 else if (strcmp(dumm1,"ar") == 0 || strcmp(dumm1,"AR") == 0)
		   {
                     ibo = 1;
		     ib1[numbonds] = ia1;
		     ib2[numbonds] = ia2;
		     numbonds++;
		   }
                 else
                    ibo = atoi(dumm1);
                 make_bond(ia1,ia2,ibo);
             }
         }
    }
// 
    quick_type();
     get_rings(natom,atom.flags,atom.atomnum,atom.iat,atom.bo);
     if (numbonds > 0)
       mopaco(numbonds,ib1,ib2);
     
     set_connectivity();
     // could check for formal charges here if we assume complete structure
     assign_aromaticity();
     quick_type();
     free_ivector(ib1,0,nibond);
     free_ivector(ib2,0,nibond);
     return TRUE;  
}
// ========================
void write_mol2()
{
    FILE *wfile = pcmoutfile;
    int i,  j, nbond;
    char number[4],atname[7],atomtype[4];

    nbond = 0;
    /*     **  calculate the number of bonds in the molecule ** */
    for( j = 1; j <= natom; j++ )
    {
      for( i = 0; i < atom.nconnect[j]; i++ )
      {
         if( atom.iat[j][i] != 0 )
         {
            if( atom.iat[j][i] < j )
               nbond = nbond + 1;
         }
       }
     }

    /*     *** write the input file for sybyl *** */
    fprintf(wfile, "@<TRIPOS>MOLECULE\n");
    fprintf(wfile,"%s\n",Struct_Title);
    fprintf(wfile, " %d %d 0 0 0\n", natom, nbond);
    fprintf(wfile, "SMALL\n");
    // charge type : None, MMFF94, Gasteiger, USER (for MMX)
    fprintf(wfile, "USER_CHARGES\n");
    fprintf(wfile, "\n");

    fprintf(wfile, "@<TRIPOS>ATOM\n");
//
    for( i = 1; i <= natom; i++ )
    {
        sprintf(number,"%d",i);
        strcpy(atname,sybname[atom.atomnum[i]-1]);
                
       fprintf(wfile,"    %3d %-5s      %10.4f%10.4f%10.4f %-5s  1  LIG  %7.3f \n",
             i,atname,atom.x[i], atom.y[i], atom.z[i],atomtype,atom.charge[i] );
    }

    fprintf(wfile, "@<TRIPOS>BOND\n");

     nbond = 0;
     for( j = 1; j <= natom; j++ )
     {
       for( i = 0; i < MAXIAT; i++ )
       {
         if( atom.iat[j][i] != 0 )
         {
            if( atom.iat[j][i] < j )
            {
               nbond = nbond + 1;
               fprintf(wfile,"%d %d %d %d \n",nbond,j,atom.iat[j][i],atom.bo[j][i]);
             }
          }
        }
     }
     return;
} 
