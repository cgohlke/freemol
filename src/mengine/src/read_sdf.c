/* NOTICE: this source code file has been modified for use with FreeMOL */
#define EXTERN extern

#include "pcwin.h"
#include "pcmod.h"

#include "fix.h"
#include "draw.h"
#include "job_control.h"
#include "vibrate.h"
#include "utility.h"
#include "rings.h"

#define TABULATOR_INCLUDE_IMPLEMENTATION
#include "tabulator.h"

EXTERN struct t_files {
        int nfiles, append, batch, icurrent, ibatno;
        }       files;
EXTERN struct ElementType { char symbol[3];
                             int atomnum;
                             float weight, covradius, vdwradius;
                              int s,p,d,f,type ;
                            } Elements[];
EXTERN struct t_logp {
        float logp;
        } logp_calc;
EXTERN struct t_solvent {
    int type;
    double EPSin, EPSsolv;
    double doffset, p1,p2,p3,p4,p5;
    double *shct,*asolv,*rsolv,*vsolv,*gpol,*rborn;
    } solvent;

char Struct_Title[100];

// ============================
void quick_type(void);
int read_sdf(int,int);
int rd_sdf(FILE *);
void write_sdf(int);
void hdel(int);
FILE * fopen_path ( char * , char * , char * ) ;
int FetchRecord(FILE *, char *);
void avgleg(void);
void get_tag(char *,char *);
void mopaco(int numbonds,int *ib1,int *ib2);
void get_molecule_memory(int);
void get_rings(int natom,long int *flags,int *atomnum,int **iat,int **bo);
double get_dipole_moment(void);
double get_total_energy(void);
char *get_structure_title(void);
void set_connectivity(void);
void compute_valency(int *,int *);
void assign_aromaticity(void);
int aromatic_5(int *array,long int *flags,int *atomnum,int *nconnect,int *input_charge,int **,int **);
int aromatic_6(int *,long int *flags,int *atomnum,int *nconnect,int *input_charge,int **,int **);
void unmark_arobond(int,int,int,int *,int *);
// ==================================
char *get_structure_title()
{
  if (strlen(Struct_Title) > 0)
    return Struct_Title;
  else 
    return NULL;
}
// ==================================
/*
 * flags - indicates whether dipole, xlogp or vibrational calculations were done
 * and if so write them to the output file. Look at the definitions in pcmod.h
 *
 */
void write_sdf(int flags)
{
    int i,j,lptest,nbond,junk;
    FILE *wfile = pcmoutfile;
    
    lptest = 0;
    for( i = 1; i <= natom; i++ )
    {
        if( atom.mmx_type[i] == 20 )
        {
            lptest = 1;
            hdel( lptest );
            break;
        }
    }
    nbond = 0;
    /*     **  calculate the number of bonds in the molecule ** */
    for( j = 1; j <= natom; j++ )
    {
      for( i = 0; i < MAXIAT; i++ )
      {
         if( atom.iat[j][i] != 0 )
         {
            if( atom.iat[j][i] < j )
               nbond = nbond + 1;
         }
       }
     }
    /*     now write the concord file */
/*
    if (files.append)
        wfile = fopen_path(Savebox.path,Savebox.fname,"a");
    else
        wfile = fopen_path(Savebox.path,Savebox.fname,"w");
*/

    j = strlen(Struct_Title);
    for (i=0; i < j; i++)
    {
        if (Struct_Title[i] == '\n')
        {
            Struct_Title[i] = '\0';
            break;
        }
    }
    fprintf(wfile,"%s\n",Struct_Title);
    fprintf(wfile," MENGINE  v1.0   1.00000     0.00000\n");       
    fprintf(wfile,"\n");       

    if (natom > 999) 
      {
	junk = 0;
	fprintf(wfile,"%3d%3d  0  0  0  0              1 V2000\n",junk,nbond);
	goto WRITE_TAB;
      } else
	fprintf(wfile,"%3d%3d  0  0  0  0              1 V2000\n",natom,nbond);


    for (i=1; i <= natom; i++)
    {
        junk = 0;
	if (atom.input_charge[i] == 1)
	  junk = 3;
	else if (atom.input_charge[i] == 2)
	  junk = 2;
	else if (atom.input_charge[i] == -1)
	  junk = 5;
	else if (atom.input_charge[i] == -2)
	  junk = 6;

        fprintf(wfile,"%10.4f%10.4f%10.4f %-3s 0  %d  0  0  0  0 \n",atom.x[i], atom.y[i],
          atom.z[i],Elements[atom.atomnum[i]-1].symbol,junk);
    }
    for (i=1; i <= natom; i++)
    {
        for (j=0; j < MAXIAT; j++)
        {
            if (atom.iat[i][j] != 0 && i < atom.iat[i][j])
               fprintf(wfile,"%3d%3d%3d  0\n",i, atom.iat[i][j], atom.bo[i][j]);
        }
    }
    //
 WRITE_TAB:
    if (natom > 999)
      {
	fprintf(wfile,"> <ATOMS>\n");
	fprintf(wfile,"+   X  Y  Z  SYMBOL FORMAL_CHARGE\n");
	for (i=1; i <= natom; i++)
	  fprintf(wfile,"|  %f  %f  %f  %s  %d\n",atom.x[i],atom.y[i],atom.z[i],Elements[atom.atomnum[i]-1].symbol,atom.input_charge[i]);
      fprintf(wfile,"\n"); 

	fprintf(wfile,"> <BONDS> \n");
	fprintf(wfile,"+  ATOM1  ATOM2 ORDER \n");
	for (i=1; i <= natom; i++)
	  {
	    for (j=0; j < atom.nconnect[i]; j++)
	      {
		if (i < atom.iat[i][j])
		  fprintf(wfile,"| %d  %d  %d\n",i,atom.iat[i][j],atom.bo[i][j]);
	      }
	  }
	fprintf(wfile,"\n");
      }

    fprintf(wfile,"M  END\n\n");
    fprintf(wfile,"> <title>\n");
    fprintf(wfile,"%s\n\n",Struct_Title);
    
    fprintf(wfile,"> <MMFF94 energy>\n");
    fprintf(wfile,"%f\n\n",get_total_energy());

    /*   fprintf(wfile,"> <MMFF94_Atomtypes>\n");
    fprintf(wfile,"+   Atom   Atomtype\n");
    for (i=1; i <=natom; i++)
      fprintf(wfile,"| %d %d\n",i,atom.type[i]);
      fprintf(wfile,"\n"); */

    if (flags & DO_DIPOLE) {
      fprintf(wfile,"> <dipole moment>\n");
      fprintf(wfile,"%f\n\n",get_dipole_moment());
    }

    if (flags & DO_XLOGP) {
      fprintf(wfile,"> <xLogP>\n");
      fprintf(wfile,"%f\n\n",logp_calc.logp);
    }

    if (flags & DO_VIBRATION) {
      fprintf(wfile,"> <Point Group>\n");
      fprintf(wfile,"%s\n\n",vibdata.ptgrp);
    
      fprintf(wfile,"> <Moments of Inertia>\n");
      fprintf(wfile,"%f8.3  %f8.3  %f8.3 \n\n",vibdata.mom_ix,vibdata.mom_iy,vibdata.mom_iz);

      fprintf(wfile,"> <thermodynamics dE, dH, S, dG, CP>\n");
      fprintf(wfile,"%f8.3  %f8.3  %f8.3  %f8.3  %f8.3  \n\n",vibdata.etot, vibdata.htot,
              vibdata.stot,vibdata.gtot,vibdata.cptot);
    }

    if (job_control.monitor)
      fprintf(wfile,"> <Intermediate_Structure>\n\n");

    fprintf(wfile,"\n");
    fprintf(wfile,"$$$$\n");
//    fclose(wfile);
}
/* ===================================  */
// fast read - assume file is open and positioned
// read structure and down to end $$$$
// and return
int rd_sdf(FILE *rfile)
{
    int i, j, niatom, nibond, ia1, ia2, ia3, ia4, ibond, itype, newatom;
    int junk,junk1,got_title,istereo;
    int n_row;
    int has_Aromatic, xPlus,yPlus,zPlus, planar;
    int Ret_Val;
    int numbonds, *ib1, *ib2;
    float xtmp, ytmp, ztmp;
    float fconst,min,max;
    char  c1[4],c2[4];
    char  atomchar[3];
    char  inputline[150];
    char  tag[30];
    char ***tab;

    numbonds = 0;
    ib1 = NULL;
    ib2 = NULL;
    has_Aromatic = FALSE;
     Ret_Val = TRUE;
     got_title = FALSE;
     xPlus = yPlus = zPlus = FALSE;
     if ( 0 == FetchRecord(rfile,inputline) )return -1;
     strncpy(Struct_Title, inputline, sizeof(Struct_Title));
     got_title = TRUE;
     FetchRecord(rfile,inputline); // blank line
     FetchRecord(rfile,inputline); // blank line
     
     FetchRecord(rfile,inputline); // natom and nbond
     for (i=0; i <4; i++)
     {
        c1[i] = inputline[i];
        c2[i] = inputline[i+3];
     }
     c1[3] = '\0';c2[3] = '\0';
     niatom = atoi(c1);
     nibond = atoi(c2); 
     if (niatom == 0)  // check for tabulator input firs
       goto READ_TAB;
     // get memory for molecule
     get_molecule_memory(niatom);
     //
     for (i=0; i < niatom; i++)
     {
        FetchRecord(rfile,inputline);
        sscanf(inputline,"%f %f %f %s %d %d",&xtmp, &ytmp, &ztmp, atomchar,&junk,&junk1);
        itype = 0;
        if (xtmp != 0.0) xPlus= TRUE;
        if (ytmp != 0.0) yPlus= TRUE;
        if (ztmp != 0.0) zPlus= TRUE;

	if (strlen(atomchar) == 2)  // two character atom symbol
	  atomchar[1] = tolower(atomchar[1]);

        newatom = make_atom(itype,xtmp,ytmp,ztmp,atomchar);
        if (newatom == -1)
        {
          fprintf(pcmlogfile, "Atom %d %s %s not recognized structure skipped\n",i,atomchar,Struct_Title);
	  natom = niatom; // allocated space for niatom need to deallocate upon return
          // read to end of structure then return
	  while (FetchRecord(rfile,inputline))
	    {
	      if (strncasecmp(inputline,"$$$$",4) == 0)
		{
		  return -2;
		}
	    }
        }
	if (junk1 == 3)
	  {
	    atom.formal_charge[newatom] = 1.0;
	    atom.input_charge[newatom] = 1;
	  } else if (junk1 == 5)
	  {
	    atom.formal_charge[newatom] = -1.0;
	    atom.input_charge[newatom] = -1;
	  } else if (junk1 == 2)
	  {
	    atom.formal_charge[newatom] = 2.0;
	    atom.input_charge[newatom] = 2;
	  } else if (junk1 == 1)
	  {
	    atom.formal_charge[newatom] = 3.0;
	    atom.input_charge[newatom] = 3;
	  } else if (junk1 == 6)
	  {
	    atom.formal_charge[newatom] = -2.0;
	    atom.input_charge[newatom] = -2;
	  }

     }
     has_Aromatic = FALSE;
     
     if (natom > niatom)
       {
	 printf("Error in read_sdf natom %d != niatom %d\n",natom,niatom);
	 printf("Structure: %s\n",get_structure_title());
	 exit(0);
       }
     planar = TRUE;
     if (xPlus && yPlus && zPlus) planar = FALSE;

     // allocate space for aromatic bonds
     numbonds = 0;
     ib1 = ivector(0,nibond);
     ib2 = ivector(0,nibond);
     
     for (i=0; i < nibond; i++)
     {
        FetchRecord(rfile,inputline);
        for (j=0; j <4; j++)
        {
           c1[j] = inputline[j];
           c2[j] = inputline[j+3];
         }  
        c1[3] = '\0';c2[3] = '\0';
        ia1 = atoi(c1);
        ia2 = atoi(c2);
        istereo = 0;
	if (ia1 > natom || ia2 > natom)
	  {
	    printf("error reading bonds in read sdf: natom %d  ia1 %d ia2 %d\n",natom,ia1,ia2);
	  }
        sscanf(inputline,"%3d%3d%3d%3d",&junk, &junk1, &ibond,&istereo);
        if (ibond >= 4)
        {
// use bonds array as temporary storage for aromatic bonds
           make_bond(ia1,ia2,1);
           ib1[numbonds] = ia1;
           ib2[numbonds] = ia2;
           numbonds++;
           has_Aromatic = TRUE;
        } else
           make_bond(ia1, ia2, ibond);
        if (istereo == 1 && planar)
        {
            atom.z[ia2] += 0.7;
            if (atom.atomnum[ia2] == 1)
              atom.type[ia2] = 60;
        }
        if (istereo == 6 && planar)
        {
            atom.z[ia2] -= 0.7;
            if (atom.atomnum[ia2] == 1)
             atom.type[ia2] = 60;
        }
     }
// read to end of structure
// parse any tags here
// agreed tags
//  <FIXED_ATOMS>
//  <RESTRAINED_ATOMS>
//  <RESTRAINED_DISTANCES>
//  <RESTRAINED_ANGLES>
//  <RESTRAINED_DIHEDRALS>
//  <ELECTROSTATICS>
//  <PARTIAL_CHARGES>
// <MONITOR>
// <ATOMS>
 READ_TAB:
     while (FetchRecord(rfile,inputline))
     {
         if (strncasecmp(inputline,"$$$$",4) == 0)
         {
            goto L_DONE;
         } else if (inputline[0] == '>')
         {
           get_tag(inputline,tag);
           if (strcasecmp(tag,"fixed_atoms") == 0)
             {
               tab = tabulator_new_from_file_using_header(rfile,"ATOM",0);
               if (tab)
                 {
                   n_row = tabulator_height(tab);
                   for (i=0; i < n_row; i++)
                     {
                       ia1 = atoi(tab[i][0]);
                       if (ia1 > 0 && ia1 <= natom)
                         {
                           fx_atom.katom_fix[fx_atom.natom_fix] = ia1;
                           fx_atom.natom_fix++;
                         }
                     }
                 }
	     } else if (strcasecmp(tag,"atoms") == 0)
	     {
               tab = tabulator_new_from_file_using_header(rfile,"X Y Z SYMBOL FORMAL_CHARGE",0);
               if (tab)
                 {
		   itype = 0;
                   n_row = tabulator_height(tab);  // n_row should be the number of entries and thus the number of atoms
		   get_molecule_memory(n_row);
		   for (i=0; i < n_row; i++)
		     {
		       xtmp = atof(tab[i][0]);
		       ytmp = atof(tab[i][1]);
		       ztmp = atof(tab[i][2]);
		       strcpy(atomchar,tab[i][3]);
		       newatom = make_atom(itype,xtmp,ytmp,ztmp,atomchar);
		       atom.formal_charge[newatom] = atoi(tab[i][4]);
		     }
		 }
	     } else if (strcasecmp(tag,"bonds") == 0)
	     {
               tab = tabulator_new_from_file_using_header(rfile,"ATOM1 ATOM2 ORDER",0);
               if (tab)
                 {
                   n_row = tabulator_height(tab);  // n_row should be the number of entries and thus the number of bonds
		   nibond = n_row;
		   numbonds = 0;
		   ib1 = ivector(0,nibond);
		   ib2 = ivector(0,nibond);
		   for (i=0; i < n_row; i++)
		     {
		       ia1 = atoi(tab[i][0]);
		       ia2 = atoi(tab[i][1]);
		       ibond = atoi(tab[i][2]);
		       if (ibond >= 4)
			 {
			   // use bonds array as temporary storage for aromatic bonds
			   make_bond(ia1,ia2,1);
			   ib1[numbonds] = ia1;
			   ib2[numbonds] = ia2;
			   numbonds++;
			   has_Aromatic = TRUE;
			 } else
			 make_bond(ia1, ia2, ibond);
		     }
		 }
             } else if (strcasecmp(tag,"restrained_atoms") == 0)
             {
               tab = tabulator_new_from_file_using_header(rfile,"ATOM MAX F_CONST X Y Z",0);
               if (tab)
                 {
                   n_row = tabulator_height(tab);
                   for (i=0; i < n_row; i++)
                     {
                       ia1 = atoi(tab[i][0]);
                       max = atof(tab[i][1]);
                       fconst = atof(tab[i][2]);
                       if (ia1 > 0 && ia1 <= natom)
                         {
                           restrain_atom.katom_restrain[restrain_atom.natom_restrain] = ia1;
                           restrain_atom.restrain_const[restrain_atom.natom_restrain] = fconst;
                           restrain_atom.restrain_max[restrain_atom.natom_restrain] = max;
                           restrain_atom.restrain_position[restrain_atom.natom_restrain][0] = atom.x[ia1]; // default positions
                           restrain_atom.restrain_position[restrain_atom.natom_restrain][1] = atom.y[ia1];
                           restrain_atom.restrain_position[restrain_atom.natom_restrain][2] = atom.z[ia1];
                           if (strlen(tab[i][3]) > 0)  // x postion
                             {
                               xtmp = atof(tab[i][3]);
                               restrain_atom.restrain_position[restrain_atom.natom_restrain][0] = xtmp;
                             }
                           if (strlen(tab[i][4]) > 0)  // y postion
                             {
                               xtmp = atof(tab[i][4]);
                               restrain_atom.restrain_position[restrain_atom.natom_restrain][1] = xtmp;
                             }
                           if (strlen(tab[i][5]) > 0)  // y postion
                             {
                               xtmp = atof(tab[i][5]);
                               restrain_atom.restrain_position[restrain_atom.natom_restrain][2] = xtmp;
                             }
                           restrain_atom.natom_restrain++;
                         }
                     }
                 }
             } else if (strcasecmp(tag,"restrained_distances") == 0)
             {
               tab = tabulator_new_from_file_using_header(rfile,"ATOM1 ATOM2 MIN MAX F_CONST",0);
               if (tab)
                 {
                   n_row = tabulator_height(tab);
                   for (i=0; i < n_row; i++)
                     {
                       ia1 = atoi(tab[i][0]);
                       ia2 = atoi(tab[i][1]);
                       min = atof(tab[i][2]);
                       max = atof(tab[i][3]);
                       fconst = atof(tab[i][4]);
                       if ((ia1 > 0 && ia1 <= natom) && (ia2 > 0 && ia2 < natom))
                         {
                           fx_dist.kdfix[fx_dist.ndfix][0] = ia1;
                           fx_dist.kdfix[fx_dist.ndfix][1] = ia2;
                           fx_dist.fdconst[fx_dist.ndfix] = fconst;
                           fx_dist.min_dist[fx_dist.ndfix] = min;
                           fx_dist.max_dist[fx_dist.ndfix] = max;
                           fx_dist.ndfix++;
                         }
                     }
                 }
             } else if (strcasecmp(tag,"restrained_angles") == 0)
             {
               tab = tabulator_new_from_file_using_header(rfile,"ATOM1 ATOM2 ATOM3 MIN MAX F_CONST",0);
               if (tab)
                 {
                   n_row = tabulator_height(tab);
                   for (i=0; i < n_row; i++)
                     {
                       ia1 = atoi(tab[i][0]);
                       ia2 = atoi(tab[i][1]);
                       ia3 = atoi(tab[i][2]);
                       min = atof(tab[i][3]);
                       max = atof(tab[i][4]);
                       fconst = atof(tab[i][5]);
                       if ((ia1 > 0 && ia1 <= natom) && (ia2 > 0 && ia2 < natom) && (ia3 > 0 && ia3 <= natom))
                         {
                           fx_angle.kafix[fx_angle.nafix][0] = ia1;
                           fx_angle.kafix[fx_angle.nafix][1] = ia2;
                           fx_angle.kafix[fx_angle.nafix][2] = ia3;
                           fx_angle.faconst[fx_angle.nafix] = fconst;
                           fx_angle.min_ang[fx_angle.nafix] = min;
                           fx_angle.max_ang[fx_angle.nafix] = max;
                           fx_angle.nafix++;
                         }
                     }
                 }
             } else if (strcasecmp(tag,"restrained_dihedrals") == 0)
             {
               tab = tabulator_new_from_file_using_header(rfile,"ATOM1 ATOM2 ATOM3 ATOM4 MIN MAX F_CONST",0);
               if (tab)
                 {
                   n_row = tabulator_height(tab);
                   for (i=0; i < n_row; i++)
                     {
                       ia1 = atoi(tab[i][0]);
                       ia2 = atoi(tab[i][1]);
                       ia3 = atoi(tab[i][2]);
                       ia4 = atoi(tab[i][3]);
                       min = atof(tab[i][4]);
                       max = atof(tab[i][5]);
                       fconst = atof(tab[i][6]);
                       if ((ia1 > 0 && ia1 <= natom) && (ia2 > 0 && ia2 < natom) && (ia3 > 0 && ia3 <= natom) && (ia4 > 0 && ia4 <= natom))
                         {
                           fx_torsion.ktfix[fx_torsion.ntfix][0] = ia1;
                           fx_torsion.ktfix[fx_torsion.ntfix][1] = ia2;
                           fx_torsion.ktfix[fx_torsion.ntfix][2] = ia3;
                           fx_torsion.ktfix[fx_torsion.ntfix][3] = ia4;
                           fx_torsion.ftconst[fx_torsion.ntfix] = fconst;
                           fx_torsion.min_tor[fx_torsion.ntfix] = min;
                           fx_torsion.max_tor[fx_torsion.ntfix] = max;
                           fx_torsion.ntfix++;
                         }
                     }
                 }
             } else if (strcasecmp(tag,"electrostatics") == 0)
	     {
               tab = tabulator_new_from_file_using_header(rfile,"TREATMENT DIELECTRIC SCALE_FACTOR GBSA_INTERNAL GBSA_EXTERNAL",0);
               if (tabulator_height(tab))
		 {
		   if (strcmp(tab[0][0],"NONE") == 0)
		     {
		       job_control.use_charge = TRUE;
		       job_control.scale = 0.0;
		     } else if (strcmp(tab[0][0],"NORMAL") == 0)
		     {
		       if(tab[0][1])
			 units.dielec = atof(tab[0][1]);
		       else 
			 units.dielec = 1.0F;
		     } else if (strcmp(tab[0][0],"SCALED") == 0)
		     {
		       job_control.use_charge = TRUE;
		       if(tab[0][2]) 
			 job_control.scale = atof(tab[0][2]);
		       else
			 job_control.scale = 1.0F;
		     } else if (strcmp(tab[0][0],"GBSA") == 0)
		     {
		       job_control.use_gbsa = TRUE;
		       solvent.type = 1;  // STILL
		       if(tab[0][3]) 
			 solvent.EPSin = atof(tab[0][3]);
		       else
			 solvent.EPSin = 1.0F;                 
		       if(tab[0][4]) 
			 solvent.EPSsolv = atof(tab[0][4]);
		       else
			 solvent.EPSsolv = 78.30F;                 
		     }
		 }
	     } else if (strcasecmp(tag,"monitor") == 0)
	     {
               tab = tabulator_new_from_file_using_header(rfile,"UPDATE_INTERVAL",0);
               if (tab)
		 {
		   fconst = atof(tab[0][0]);
		   if (fconst > 0.0 && fconst < 10.0)
		     {
		       job_control.monitor = TRUE;
		       job_control.interval = fconst; // assume update interval is in 0.1 secs
		       if (natom < 100)  // don't do small jobs
			 job_control.monitor = FALSE;
		     }
		 }
	     }
	 }
     }
     if (natom == 0) // no atoms read
       return FALSE;
     // if (has_aromatic)
     // need to deal with aromatic bonds
 L_DONE:
     set_connectivity();
     quick_type();
     get_rings(natom,atom.flags,atom.atomnum,atom.iat,atom.bo);
     if (has_Aromatic)
     {
       mopaco(numbonds,ib1,ib2);
     }
     assign_aromaticity();
     quick_type();
     free_ivector(ib1,0,nibond);
     free_ivector(ib2,0,nibond);
     return TRUE;
}
// ===========================
void assign_aromaticity()
{
  int i,j,aromatic;
  int array[6];
  long int aromatic_mask;

  aromatic_mask = (1L << AROMATIC_MASK);
  for (i=0; i < rings.nring5; i++)
    {
      for (j=0; j < 5; j++)
	array[j] = rings.r15[i][j];
      aromatic = aromatic_5(array,atom.flags,atom.atomnum,atom.nconnect,atom.input_charge,atom.iat,atom.bo);
      if (aromatic == TRUE)
      {
         for (j = 0; j < 5; j++)
             atom.flags[array[j]] |= aromatic_mask;
      }
    }
  for (i=0; i < rings.nring6; i++)
    {
      for (j=0; j < 6; j++)
	array[j] = rings.r16[i][j];
      aromatic = aromatic_6(array,atom.flags,atom.atomnum,atom.nconnect,atom.input_charge,atom.iat,atom.bo);
      if (aromatic == TRUE)
      {
         for (j = 0; j < 6; j++)
             atom.flags[array[j]] |= aromatic_mask;
      }
    }
}
// =====================
void get_tag(char *line,char *tag)
{
  int i,j,icount;
  icount = 0;
  strcpy(tag,"");
  for (i=0; i < strlen(line); i++)
    {
      if (line[i] == '<')
        {
          for (j=i+1; j < strlen(line); j++)
            {
              if (line[j] == '>')
                {
                  tag[icount] = '\0';
                  return;
                } else
                {
                  tag[icount] = line[j];
                  icount++;
                }
            }
        }
    }
}
/* ------------------------------------*/      
// new mopaco version - input list of aromatic bonds
// test bond if part of ring, get ring, test if all atoms of ring can be aromatic
// get fused rings
// add bonds and check valency of all aromatic atoms
void mopaco(int nbond,int *ia1,int *ia2)
{
  int i,j,ia,ib;
  int ibi,ibj;
  int npi,done;
  int array[6],tbo[6];

  //
  for (i=0; i < rings.nring6; i++)
    {
      for (j=0; j < 6; j++)
	array[j] = rings.r16[i][j];
      compute_valency(array,tbo);
      npi = 0;
      for (j=0; j < 6; j++)
	{
	  if (atom.atomnum[array[j]] == 6 && tbo[j] == 3)
	    npi++;
	  else if (atom.atomnum[array[j]] == 7 && tbo[j] == 2)
	    npi++;
	}
      if (npi == 6) // no double bonds add 3
	{
	  make_bond(array[0],array[1],1);
	  unmark_arobond(array[0],array[1],nbond,ia1,ia2);
	  make_bond(array[2],array[3],1);
	  unmark_arobond(array[2],array[3],nbond,ia1,ia2);
	  make_bond(array[4],array[5],1);
	  unmark_arobond(array[4],array[5],nbond,ia1,ia2);
	} else if (npi == 4) // add two bonds
	{
	  // find four centers 
	  done = FALSE;
	  if (tbo[0] == 3 && tbo[1] == 3 && tbo[2] == 3 && tbo[3] == 3)
	    {
	      make_bond(array[0],array[1],1);
	      unmark_arobond(array[0],array[1],nbond,ia1,ia2);
	      make_bond(array[2],array[3],1);
	      unmark_arobond(array[2],array[3],nbond,ia1,ia2);
	      done = TRUE;
	    } else if (tbo[1] == 3 && tbo[2] == 3 && tbo[3] == 3 && tbo[4] == 3)
	    {
	      make_bond(array[1],array[2],1);
	      unmark_arobond(array[1],array[2],nbond,ia1,ia2);
	      make_bond(array[3],array[4],1);
	      unmark_arobond(array[3],array[4],nbond,ia1,ia2);
	      done = TRUE;
	    } else if (tbo[2] == 3 && tbo[3] == 3 && tbo[4] == 3 && tbo[5] == 3)
	    {
	      make_bond(array[2],array[3],1);
	      unmark_arobond(array[2],array[3],nbond,ia1,ia2);
	      make_bond(array[4],array[5],1);
	      unmark_arobond(array[4],array[5],nbond,ia1,ia2);
	      done = TRUE;
	    } else if (tbo[3] == 3 && tbo[4] == 3 && tbo[5] == 3 && tbo[0] == 3)
	    {
	      make_bond(array[3],array[4],1);
	      unmark_arobond(array[3],array[4],nbond,ia1,ia2);
	      make_bond(array[5],array[0],1);
	      unmark_arobond(array[5],array[0],nbond,ia1,ia2);
	      done = TRUE;
	    } else if (tbo[4] == 3 && tbo[5] == 3 && tbo[0] == 3 && tbo[1] == 3)
	    {
	      make_bond(array[0],array[1],1);
	      unmark_arobond(array[0],array[1],nbond,ia1,ia2);
	      make_bond(array[4],array[5],1);
	      unmark_arobond(array[4],array[5],nbond,ia1,ia2);
	      done = TRUE;
	    } else if (tbo[5] == 3 && tbo[0] == 3 && tbo[1] == 3 && tbo[2] == 3)
	    {
	      make_bond(array[0],array[5],1);
	      unmark_arobond(array[0],array[5],nbond,ia1,ia2);
	      make_bond(array[1],array[2],1);
	      unmark_arobond(array[1],array[2],nbond,ia1,ia2);
	      done = TRUE;
	    }
	  if (done == FALSE)
	    fprintf(pcmlogfile,"Unable to make bond for npi = 1: %d %d %d %d %d %d = %d %d %d %d %d %d\n",
		   array[0],array[1],array[2],array[3],array[4],array[5],
		   tbo[0],tbo[1],tbo[2],tbo[3],tbo[4],tbo[5]);
	} else if (npi == 2) // add one bond
	{
	  done = FALSE;
	  for (j=0; j < 5; j++)
	    {
	      if (tbo[j] == 3) // found first atoms
		{
		  if (j == 0) 
		    {
		      if (tbo[5] == 3)
			{
			  make_bond(array[0],array[5],1);
			  unmark_arobond(array[0],array[5],nbond,ia1,ia2);
			  done = TRUE;
			} else if (tbo[1] == 0)
			{
			  make_bond(array[0],array[1],1);
			  unmark_arobond(array[0],array[1],nbond,ia1,ia2);
			  done = TRUE;
			}
		    } else if (tbo[j+1] == 3)
		    {
			make_bond(array[j],array[j+1],1);
		        unmark_arobond(array[j],array[j+1],nbond,ia1,ia2);
			done = TRUE;
		    }
		}
	    }
	  if (done == FALSE)
	    fprintf(pcmlogfile,"Unable to make bond for npi = 1: %d %d %d %d %d %d = %d %d %d %d %d %d\n",
		   array[0],array[1],array[2],array[3],array[4],array[5],
		   tbo[0],tbo[1],tbo[2],tbo[3],tbo[4],tbo[5]);
	}
    }
  // make any bonds left
  for (i=0; i < nbond; i++)
    {
      if (ia1[i] != 0 && ia2[i] != 0)
	{
	  ia = ia1[i];
	  ib = ia2[i];
	  ibi = 0;
	  ibj = 0;
	  for (j=0; j < MAXIAT; j++)
	    {
	      ibi += atom.bo[ia][j];
	      ibj += atom.bo[ib][j];
	    }
	  if (atom.atomnum[ia] == 6 && atom.atomnum[ib] == 6)
	    {
	      if (ibi == 3 && ibj == 3)
		make_bond(ia,ib,1);
	    } else if (atom.atomnum[ia] == 6 && atom.atomnum[ib] == 7)
	    {
	      if (ibi == 3 && ibj == 2)
		make_bond(ia,ib,1);
	    } else if (atom.atomnum[ia] == 7 && atom.atomnum[ib] == 6)
	    {
	      if (ibi == 2 && ibj == 3)
		make_bond(ia,ib,1);
	    }
	}
    }
}
// =========================
void unmark_arobond(int ia,int ib,int nbond,int *ia1,int *ia2)
{
  int i;
  for (i=0; i < nbond; i++)
    {
      if (ia == ia1[i] && ib == ia2[i])
	{
	  ia1[i] = 0;
	  ia2[i] = 0;
	  return;
	} else if (ia == ia2[i] && ib == ia1[i])
	{
	  ia1[i] = 0;
	  ia2[i] = 0;
	  return;
	}
    }
}
// =====================
void compute_valency(int *array,int *tbo)
{
  int i;
  int ia,j;

  for (i=0; i < 6; i++)
    tbo[i] = 0;
  for (i=0; i < 6; i++)
    {
      ia = array[i];
      for (j=0; j < MAXIAT; j++)
	tbo[i] += atom.bo[ia][j];
    }
}
// =========================
int all_aromatic(int *array,int naromatic,int *aro_atoms)
{
  int i,j,ia,found;

  for (i=0; i < 6; i++)
    {
      ia = array[i];
      found = FALSE;
      for (j=0; j < naromatic; j++)
	{
	  if (ia == aro_atoms[j])
	    {
	      found = TRUE;
	      break;
	    }
	}
      if (found == FALSE)
	return FALSE;
    }
  return TRUE;
}

