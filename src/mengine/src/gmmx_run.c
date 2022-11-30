/* NOTICE: this source code file has been modified for use with FreeMOL */
#define notKEG_DEBUG
#define notMMFF_DEBUG


#define EXTERN extern
#include "pcwin.h"
#include "pcmod.h"
#include "vibrate.h"

#include <errno.h>

struct t_logp {
        float logp;
        } logp_calc;
    
EXTERN struct t_files {
        int nfiles, append, batch, icurrent, ibatno;
        }       files;

EXTERN struct t_minim_values {
        int iprint, ndc, nconst;
        float dielc;
        } minim_values;
EXTERN struct t_minim_control {
        int type, method, field, added_const;
        char added_path[256],added_name[256];
        } minim_control;

int rd_sdf(FILE *);
int setup_calculation(void);
void end_calculation(void);
void minimize(int natom,int *use,double *x,double *y,double *z);
static int close_contact(int);
void initialize(void);
void type(void);
void read_gmmxinp(char*,FILE*,int);
void hdel(int);
void hadd(void);
int FetchRecord(FILE *, char *);
void write_sdf(int);
void zero_data(void);
void type_mmx(void);
void read_datafiles(char *);
float xlogp(int natom,int *atomnum,int **iat,int **bo,long int *flags,int *type);
void print_energy_totals(void);
char *get_structure_title(void);
void vibrate(int natom,double *atomwt,double *x,double *y,double *z,double *charge);
void charge_dipole(int natom,double *x,double *y,double *z,double *atomwt,double *charge);
void set_connectivity(void);
char *get_structure_title(void);
double get_total_energy(void);

// =====================================
void read_gmmxinp(char *paramfile,
		  FILE *infile,
		  int flags)
{
   int i,ncount,nret;
   int icount;
   int ngood, nbparam, nbatom,nbcontact;
   float logp;
   //char *std_file;
   FILE *logfile;
   FILE *efile = NULL;


   if (flags & DO_TEST)
     efile = fopen(testfilename,"w");


// check filetype
    ncount = 0;
    
    // open logfile
    logfile = pcmlogfile;
    
    files.nfiles = 0;
//
    minim_control.method = 3;  // just do first deriv minimization
    icount = 0;
//

    ngood = nbparam = nbatom = nbcontact = 0;

    if (flags & DO_DEBUG)
      {
	for (i=1; i<= natom; i++)
	  printf("Atom: %d  type: %d  flags: %ld\n",i,atom.type[i],atom.flags[i]);
	minim_values.iprint = TRUE;
      }         
    nret = setup_calculation();

    // debug run - exit after printing all interactions and parameters
    if (flags & DO_DEBUG)
      {
	 print_energy_totals();
	 exit(0);
      }

    if (nret == TRUE)
    {
       minimize(natom,atom.use,atom.x,atom.y,atom.z);
       if (flags & DO_ADDH)
	 {
	   end_calculation();
	   hdel(0);
	   hadd();
	   set_connectivity();
	   type();
	   //	   type_mmx();
	   nret = setup_calculation();
	   if (nret == TRUE)
		minimize(natom,atom.use,atom.x,atom.y,atom.z);
	   else 
	     {
	       end_calculation();
	       printf("Failed setup after hydrogen addition\n");
	       return;
	     }
	 }
       // compute xlogp
       if (flags & DO_XLOGP) {
	   if (VERBOSE) fprintf(pcmlogfile, "  Evaluating XlogP\n");
           logp = xlogp(natom,atom.atomnum,atom.iat,atom.bo,atom.flags,atom.type);
           logp_calc.logp = logp;
	}
	// compute dipole moment
	if (flags & DO_DIPOLE) {
	   if (VERBOSE) fprintf(pcmlogfile, "  Evaluating dipole moment\n");
  	   charge_dipole(natom,atom.x,atom.y,atom.z,atom.atomwt,atom.charge);
	}
	// compute vib
	if (flags & DO_VIBRATION) {
	   if (VERBOSE) fprintf(pcmlogfile, "  Evaluating vibrational parameters\n");		      
	   vibrate(natom,atom.atomwt,atom.x,atom.y,atom.z,atom.charge);
	}
	//	write_sdf(flags);
	++ngood;
	++files.nfiles;
	if (flags & DO_TEST)
	   fprintf(efile,"%s %f\n",get_structure_title(),get_total_energy());
    } else
    {
	nbparam++;
	fprintf(logfile,"Strnum: %d  CID: %s fails on missing parameters\n",files.nfiles,get_structure_title());
	if (VERBOSE) 
	  fprintf(stdout, "Missing parameters - no calc - %s\n",get_structure_title());
    }
    end_calculation();
    if (VERBOSE) {
      fprintf(stdout,"Strnum: %d %s natom %d \n",files.nfiles,get_structure_title(), natom);
    }
   if (flags & DO_TEST)
     fclose(efile);
}
// ================================================
int close_contact(int mode)
{
  int i, j;
  double dx,dy,dz;

  for (i=1; i < natom; i++)
  {
    for (j=i+1; j <= natom; j++)
      {
        dx = fabs(atom.x[i] - atom.x[j]);
        dy = fabs(atom.y[i] - atom.y[j]);
        dz = fabs(atom.z[i] - atom.z[j]);
        if ( (dx + dy + dz) < 0.1)
        {
            return FALSE;
        }
      }
  }
    return TRUE;
}

