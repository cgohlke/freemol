/* NOTICE: this source code file has been modified for use with FreeMOL */

#define notKEG_DEBUG

#define EXTERN
#include "pcwin.h"
#include "pcmod.h"
#include "angles.h"
#include "attached.h"
#include "torsions.h"
#include "nonbond.h"
#include "bonds_ff.h"
#include "derivs.h"
#include "hess.h"
#include "atom_k.h"
#include "cutoffs.h"
#include "fix.h"
#include "vibrate.h"
#include "job_control.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
//#include <malloc.h>

struct t_minim_values {
        int iprint, ndc, nconst;
        float dielc;
        } minim_values;
        
struct t_minim_control {
        int type, method, field, added_const;
        char added_path[256],added_name[256];
        } minim_control;
EXTERN struct t_files {
        int nfiles, append, batch, icurrent, ibatno;
        }       files;
char Savename[80];
       
/*
void initialize_pcmodel(char*);
void mmp22mod(int,int);
void pcmfin(int, int);
void initialize(void);
void mmxsub(int);
void check_numfile(int, char *);
void mac2mod(int, int);
void bbchk(void);
void type(void);
void initialize_gmmx(void);
void read_gmmxinp(char*, char*, char*, FILE*, int);
void run_gmmx(void);
void search_rings(int);
*/
void initialize_pcmodel(char*);
void initialize(void);
void read_gmmxinp(char*, FILE*, int);
void set_connectivity(void);
void type(void);
int rd_sdf(FILE *);
int read_mol2(FILE *);
int read_cml(char *);
void write_sdf(int);
void zero_data(void);
void read_datafiles(char *);
char *get_structure_title(void);
void free_molecule_memory(void);

void usage(void) {
  fprintf(stderr,"\nUsage: mengine [options]\n Default reads sdf on stdin, writes on stdout\n\n"
         "-h\tThis help page\n"
	 "-i\tInput filename\n"
	 "-fi\tInput filetype\n"
	 "-o\tOutput filename\n"
	 "-fo\tOutput filetype\n"
         "-v\tVerbose output\n"
	 "-t\tFilename for test suite output\n"
	 "-debug\tDump atoms types and parameters for input file\n"
         "-a\tAdd hydrogen atoms; delete existing ones\n"
         "-d\tEvaluate dipole moment\n"
         "-x\tEvaluate XlogP\n"
         "-vi\tEvaluate vibrational data\n\n");
  fprintf(stderr,"Input file formats:\n");
  fprintf(stderr,"MDL standard data  sdf   1\n");
  fprintf(stderr,"Sybyl Mol2         mol2  2\n");
  fprintf(stderr,"Chemical Markup    cml   3\n\n");
  fprintf(stderr,"Output file formats:\n");
  fprintf(stderr,"MDL standard data  sdf   1\n");
  fprintf(stderr,"Sybyl Mol2         mol2  2\n");
  fprintf(stderr,"Chemical Markup    cml   3\n");
  fprintf(stderr,"Protein Databank   pdb   4\n\n");
}

/* ==================================================== */
int main(int argc, char *argv[])
{
     FILE *infile;

    char infilename[50];
    char outfilename[50];
    char *mmff94filename = NULL;
    char *mmxfilename = NULL;                  
    char *logfilename = NULL;

    char inftype[10];  
    char outftype[10];
    int infiletype = 0;
    int outfiletype = 0;  
    int i,nret,icount;
    int flags = 0;

    // defaults
    strcpy(inftype,"sdf");
    strcpy(outftype,"sdf");
    infiletype = 1;
    outfiletype = 1;
    strcpy(infilename,"");
    strcpy(outfilename,"");

    if (argc == 2 && strcmp(argv[1],"-h") == 0)
	{
	  printf("argv %s\n",argv[1]);
	  usage();
	  exit(1);
	}

    for (i=1; i < argc ; i++)
      {

	if (strcmp(argv[i],"-a") == 0)
	  flags = flags | DO_ADDH;
	else if (strcmp(argv[i],"-d") == 0)
	  flags = flags | DO_DIPOLE;
	else if (strcmp(argv[i],"-x") == 0)
	  flags = flags | DO_XLOGP;
	else if (strcmp(argv[i],"-vi") == 0)
	  flags = flags | DO_VIBRATION;
	else if (strcmp(argv[i],"-v") == 0)
	  VERBOSE = 1;
	else if (strcmp(argv[i],"-t") == 0)
	  {
	    strcpy(testfilename,argv[i+1]);
	    flags = flags | DO_TEST;
	  }
	else if (strcmp(argv[i],"-debug") == 0)
	  {
	    flags = flags | DO_DEBUG;
	  }
	else if (strcmp(argv[i],"-h") == 0)
	  {
	    usage();
	    exit(1);
	  }
	else if (strcmp(argv[i],"-i") == 0)
	  {
	    strcpy(infilename,argv[i+1]);
	    i++;
	  }
	else if (strcmp(argv[i],"-fi") == 0)
	  {
	    strcpy(inftype,argv[i+1]);
	    i++;
	  }
	else if (strcmp(argv[i],"-o") == 0)
	  {
	    strcpy(outfilename,argv[i+1]);
	    i++;
	  }
	else if (strcmp(argv[i],"-fo") == 0)
	  {
	    strcpy(outftype,argv[i+1]);
	    i++;
	  }
      }
	
#ifdef KEG_DEBUG
    strcpy(infilename,"Components-pub.sdf");
    strcpy(outfilename,"output.sdf");
#endif

    /* default behavior is to use the hardcoded parameter files */

    if(!mmff94filename) mmff94filename = strdup("|mmff94"); /* hard-coded builtins */
    if(!mmxfilename) mmxfilename = strdup("|mmxconst"); /* hard-coded builtins */

    /* if no input file provide, then use standard input */

    if(strlen(infilename) > 0) 
      {
	infile = fopen(infilename,"rb");
	if ( (strcasecmp(inftype,"sdf") == 0) || (strcmp(inftype,"1") == 0)) 
	  infiletype = 1;
	else if ((strcasecmp(inftype,"mol2") == 0) || (strcmp(inftype,"2") == 0)) 
	  infiletype = 2;
	else if ((strcasecmp(inftype,"cml") == 0) || (strcmp(inftype,"3") == 0)) 
	  infiletype = 3;
      } else 
      {
	infile = stdin;
      }

    /* if no output file provided, then use standard output */

    if(strlen(outfilename) > 0) 
     {
       pcmoutfile = fopen(outfilename,"wb");
       if ((strcasecmp(outftype,"sdf") == 0) || (strcmp(outftype,"1") == 0))
	   outfiletype = 1;
       else if ((strcasecmp(outftype,"mol2") == 0) || (strcmp(outftype,"2") == 0))
	   outfiletype = 2;
       else if ((strcasecmp(outftype,"cml") == 0) || (strcmp(outftype,"3") == 0))
	   outfiletype = 3;
       else if ((strcasecmp(outftype,"pdb") == 0) || (strcmp(outftype,"4") == 0))
	   outfiletype = 4;
     } else 
     {
	pcmoutfile = stdout;
     }

    /* if no log file provided, then use standard error */

    if(logfilename) {
      pcmlogfile = fopen(logfilename,"wb");
    } else {
      pcmlogfile = stderr;
    }
    // end of filename setup

    /* initialize the system using parameter file */

    initialize_pcmodel(mmxfilename);
    minim_values.iprint = 0; //FALSE;
    minim_control.field = MMFF94;
    zero_data();
    read_datafiles(mmff94filename);  
    /* run the job */

    icount = 0;
    while (TRUE)
      {
	initialize();  // clear out atom data
	// read routine reads molecular data and allocates space for new molecule
	if (infiletype == 1)  // SDF
	  {
	    nret = rd_sdf(infile);
	    if (nret == -2) 
	      goto GET_NEXT;
	    else if (nret == -1)
	      goto DONE;
	  } else if (infiletype == 2)  // MOL2
	  {
	    nret = read_mol2(infile);
	    if (nret == -1) break;
	  } else if (infiletype == 3)  // CML
	  {
	    if (strlen(infilename) == 0)
	      {
		printf("Error - cml input requires local file\n");
		exit(0);
	      }
	    fclose(infile);
	    nret = read_cml(infilename);
	    if (nret == -1) break;
	  }
    // run job
	icount++;
	//	printf("Processing %d %s\n",icount,get_structure_title());
	set_connectivity();
	type();
        read_gmmxinp(mmff94filename, infile, flags);
	
    // output new structure
	if (outfiletype == 1) // SDF
	  write_sdf(flags);
	else if (outfiletype == 2) // MOL2
	  continue;
	else if (outfiletype == 3) // CML
	  continue;
	else if (outfiletype == 4) // PDB
	  continue;
	if (infiletype != 1)
	  break;
	// free molecule memory
      GET_NEXT:
	free_molecule_memory();
	continue;
      }
    //    if (infile) {
    //  read_gmmxinp(mmff94filename, infile, flags);
    /* release temporary strings */
 DONE:
    if(mmff94filename) free(mmff94filename);
    if(mmxfilename) free(mmxfilename);
    if(logfilename) free(logfilename);

    /* close file handles */

    if(infile && (infile!=stdin)) fclose(infile);
    if(pcmoutfile && (pcmoutfile!=stdout)) fclose(pcmoutfile);
    if(pcmlogfile && (pcmlogfile!=stderr)) fclose(pcmlogfile);

    exit(0);
}

