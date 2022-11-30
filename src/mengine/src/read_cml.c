#define EXTERN extern

#include "pcwin.h"
#include "pcmod.h"

EXTERN struct t_files {
        int nfiles, append, batch, icurrent, ibatno;
        }       files;
EXTERN struct ElementType { char symbol[3];
                             int atomnum;
                             float weight, extweight, covradius, vdwradius;
                              int s,p,d,f,type ;
                            } Elements[];

EXTERN char *Struct_Title;

void write_cml(void);
int make_atom(int,float,float,float,char *);
void initialize(void);
void make_bond(int,int,int);
void message_alert(char *, char *);
FILE * fopen_path ( char * , char * , char * ) ;
int FetchRecord(FILE *, char *);
int read_cml(char *);
void get_molecule_memory(int);
char *get_structure_title(void);

/* ------------------------------- */
int tokenizer(int n,char *string,char *attr, char *value)
{
    int i,j;
    int icount, jcount;

    icount = jcount = 0;
    strcpy(attr,"");
    strcpy(value,"");
    
    for (i = n; i < strlen(string); i++)
    {
        if (string[i] == '>')
        {
            attr[icount] = '\0';
            value[jcount] = '\0';
            return (i+1);
        }else if (string[i] == '<')
          continue;
        else if (string[i] == ' ')
        {
            attr[icount] = '\0';
            value[jcount] = '\0';
            return (i+1);
        }else if (string[i] == '=')
        {
            attr[icount] = '\0';
            i += 2;
            for (j=i; j < strlen(string); j++)
            {
                if (string[j] == '"')
                {
                    value[jcount] = '\0';
                    return (j+2);
                } else
                {
                    value[jcount] = string[j];
                    jcount++;
                }
            }
        } else
        {
            attr[icount] = string[i];
            icount++;
        }
    }
    return(1000);      
}   
// ====================================
int read_cml(char *infilename)
{
    int i,n, pos;
    int newatom,niatm,at1,at2,ibo;
    int Makeatom,Makebond;
    char atomID[1000][5], as1[5],as2[5];
    char line[255];
    char dummy[30],dummy1[90],atomchar[3];
    float x1,y1,z1;
    FILE *infile;
    
    infile = fopen(infilename,"r");
    if (infile == NULL)
    {
        message_alert("Error opening cml file","Error");
        return -1;
    }
// check for valid XML file should be <?xml
    FetchRecord(infile,line);
    *(line+ strlen(line) ) = '\0';
    sscanf(line,"%s %s",dummy,dummy1);
    if (strcmp(dummy,"<?xml") != 0)
    {
        message_alert("This does not appear to be an XML file","Error");
        fclose(infile);
        return -1;
    }
    // count number of atoms - then close and reopen
    niatm = 0;
    ibo = 0;
 L_20:
    FetchRecord(infile,line);
    *(line + strlen(line) ) = '\0';
    if (feof(infile))
                goto L_DONE;
    n = 0;
    pos = 0;
 
    while (TRUE)
      {
	n = pos;
	pos = tokenizer(n,line,dummy,dummy1);
	if (pos == 1000) goto L_20;
	if (strcmp(dummy,"atomArray") == 0) // start of atoms
	  {
L_05:
	    n = pos = 0;
	    FetchRecord(infile,line);
L_10:
	    n = pos;
	    pos = tokenizer(n,line,dummy,dummy1);
	    if (pos == 1000) goto L_05;
	    if (strcmp(dummy,"/atomArray") == 0)
	      {
		fclose(infile);
		get_molecule_memory(niatm);
		infile = fopen(infilename,"r");
		goto L_30;
	      } else if (strcmp(dummy,"atom") == 0)
	      {
		niatm++;
		goto L_05;
	      } else 
	       goto L_10;
	  }
      }

L_30:
    if (niatm == 0)
      {
	printf("no atoms found\n");
	return -1;
      }
    FetchRecord(infile,line);
    *(line + strlen(line) ) = '\0';
    if (feof(infile))
                goto L_DONE;
    n = 0;
    pos = 0;
  // start of token read
    
L_40:
    n = pos;
    pos = tokenizer(n,line,dummy,dummy1);
    if (pos == 1000) goto L_30;
    if (strcmp(dummy,"molecule") == 0)
    {
        n = pos;
        pos = tokenizer(n,line,dummy,dummy1);
        if (strcmp(dummy,"title") == 0)
        {
            strcpy(Struct_Title,dummy1);
        }
        goto L_40;
    } else if (strcmp(dummy,"/molecule") == 0)
    {
        goto L_DONE;
    } else if (strcmp(dummy,"atomArray") == 0)
    {
        niatm = 0;
L_50:
       n = 0;
       pos = 0;
       FetchRecord(infile,line);
       *(line + strlen(line) ) = '\0';
       if (feof(infile))
           goto L_DONE;
L_51:
       n = pos;
       pos = tokenizer(n,line,dummy,dummy1);
       if (pos == 1000) goto L_50;
       if (strcmp(dummy,"/atomArray") == 0)
           goto L_30;
       else if (strcmp(dummy,"atom") == 0)
       {
          x1 = y1 = z1 = 0.0;
          Makeatom = FALSE;
L_60:
          n = pos;
          pos = tokenizer(n,line,dummy,dummy1);
          if (pos == 1000) // end of line
          {
              if (Makeatom == FALSE)
              {
                 newatom = make_atom(0,x1,y1,z1,atomchar);
                 Makeatom = TRUE;
              }
              goto L_50;
          }
          if (strcmp(dummy,"id") == 0)
          {
              strcpy(atomID[niatm],dummy1);
              niatm++;
              goto L_60;
          } else if (strcmp(dummy,"elementType") == 0)
          {
              strcpy(atomchar,dummy1);
              goto L_60;
          } else if (strcmp(dummy,"formalCharge") == 0)
          {
              goto L_60;
          } else if (strcmp(dummy,"x2") == 0)
          {
              sscanf(dummy1,"%f",&x1);
              goto L_60;
          } else if (strcmp(dummy,"x3") == 0)
          {
              sscanf(dummy1,"%f",&x1);
              goto L_60;
          } else if (strcmp(dummy,"y2") == 0)
          {
              sscanf(dummy1,"%f",&y1);
              goto L_60;
          } else if (strcmp(dummy,"y3") == 0)
          {
              sscanf(dummy1,"%f",&y1);
              goto L_60;
          } else if (strcmp(dummy,"z3") == 0)
          {
              sscanf(dummy1,"%f",&z1);
              goto L_60;
          } else if (strcmp(dummy,"/atom") == 0)
          {
              newatom = make_atom(0,x1,y1,z1,atomchar);
              goto L_50;
          } else if (strcmp(dummy,"/") == 0)
          {
              newatom = make_atom(0,x1,y1,z1,atomchar);
              Makeatom = TRUE;
              goto L_50;
          }
          else
             goto L_60;
       }
        goto L_51;  // unrecognized token
    } else if (strcmp(dummy,"bondArray") == 0)
    {
L_70:
       n = 0;
       pos = 0;
       FetchRecord(infile,line);
       *(line + strlen(line) ) = '\0';
       if (feof(infile))
           goto L_DONE;
L_71:
       n = pos;
       pos = tokenizer(n,line,dummy,dummy1);
       if (strcmp(dummy,"/bondArray") == 0)
           goto L_30;
       else if (strcmp(dummy,"bond") == 0)
       {
           at1 = at2 = 0;
           Makebond = FALSE;
L_80:
           n = pos;
           pos = tokenizer(n,line,dummy,dummy1);
           if (pos == 1000 && Makebond == FALSE)
           {
               Makebond = TRUE;
               if (at1 > 0 && at2 > 0)
                 make_bond(at1,at2,ibo);
               goto L_70;
           }
           if (strcmp(dummy,"atomRefs2") == 0)
           {
               sscanf(dummy1,"%s %s",as1,as2);
               for (i=0; i < niatm; i++)
               {
                   if (strcmp(as1,atomID[i]) == 0)
                   {
                      at1 = i+1;
                      break;
                   }
               }
               for (i=0; i < niatm; i++)
               {
                   if (strcmp(as2,atomID[i]) == 0)
                   {
                      at2 = i+1;
                      break;
                   }
               }
               goto L_80;
           } else if (strcmp(dummy,"atomRef") == 0)
           {
               sscanf(dummy1,"%s",as1);
               for (i=0; i < niatm; i++)
               {
                   if (strcmp(as1,atomID[i]) == 0)
                   {
                       if (at1 == 0)
                          at1 = i+1;
                       else
                          at2 = i+1;  
                      break;
                   }
               }
               goto L_80;               
           } else if (strcmp(dummy,"atomRef1") == 0)
           {
               sscanf(dummy1,"%s",as1);
               for (i=0; i < niatm; i++)
               {
                   if (strcmp(as1,atomID[i]) == 0)
                   {
                       at1 = i+1;
                       break;
                   }
               }
               goto L_80;
           } else if (strcmp(dummy,"atomRef2") == 0)
           {
               sscanf(dummy1,"%s",as2);
               for (i=0; i < niatm; i++)
               {
                   if (strcmp(as2,atomID[i]) == 0)
                   {
                       at2 = i+1;
                       break;
                   }
               }
               goto L_80;
           } else if (strcmp(dummy,"order") == 0)
           {
               sscanf(dummy1,"%s",as1);
               if (strcmp(as1,"S") == 0)
                 ibo = 1;
               else if (strcmp(as1,"D") == 0)
                 ibo = 2;
               else if (strcmp(as1,"A") == 0)  // aromatic bond 
                 ibo = 1;
               else
                 ibo = atoi(as1);                 
               goto L_80;
           } else if (strcmp(dummy,"/bond") == 0)
           {
               Makebond = TRUE;
               if (at1 > 0 && at2 > 0)
                 make_bond(at1,at2,ibo);
               goto L_70;
           } else
               goto L_80;
       }
       goto L_71;
    }
    goto L_40;
L_DONE:
    fclose(infile); 
    return TRUE;     
}
//  ==============================================
void write_cml()
{
    int i,j,nbond;
    int at1,at2,ibo;
    FILE *wfile = pcmoutfile;

//
    /*     **  calculate the number of bonds in the molecule ** */
    nbond = 0;
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
//
    fprintf(wfile,"<molecule title=\"%s\">\n",Struct_Title);
    fprintf(wfile,"<atomArray>\n");
    for (i=1; i<= natom; i++)
    {
        fprintf(wfile,"<atom id=\"%d\" elementType=\"%s\" x3=\"%f\" y3=\"%f\" z3=\"%f\"></atom>\n",
           i,Elements[atom.atomnum[i]-1].symbol,atom.x[i],atom.y[i],atom.z[i]);
    }
    fprintf(wfile,"</atomArray>\n");
//
    fprintf(wfile,"<bondArray>\n");
    for( j = 1; j <= natom; j++ )
    {
      for( i = 0; i < atom.nconnect[j]; i++ )
      {
         if( atom.iat[j][i] != 0 )
         {
             if (atom.iat[j][i] > j)
             {
                at1 = j;
                at2 = atom.iat[j][i];
                ibo = atom.bo[j][i];
                fprintf(wfile,"<bond atomRef1=\"%d\" atomRef2=\"%d\" order=\"%d\"></bond>\n",at1,at2,ibo);
             }
         }
       }
     }
    fprintf(wfile,"</bondArray>\n");
    fprintf(wfile,"</molecule>\n");

    fclose(wfile);
}

