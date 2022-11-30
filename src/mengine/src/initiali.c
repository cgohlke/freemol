#define EXTERN extern

#include "pcwin.h"
#include "pcmod.h"
#include "cutoffs.h"
#include "fix.h"
#include "atom_k.h"
#include "job_control.h"

void reset_atom_data(void);
void reset_calc_parameters(void);
void zero_data(void);
void read_datafiles(char *);
void initialize_pcmodel(char *);
void set_field(int);
int get_field(void);
void message_alert(char *, char *);


EXTERN struct t_files {
        int nfiles, append, batch, icurrent;
        int ibatno;
        } files;

EXTERN struct t_minim_control {
        int type, method, field, added_const;
        char added_path[256],added_name[256];
        } minim_control;

EXTERN struct t_minim_values {
        int iprint, ndc, nconst;
        float dielc;
        } minim_values;
                
struct t_user {
        int dielec;
        } user;

EXTERN struct t_pcmfile {
        char string[200];
        int head;
        char token[20];
        int state;
        unsigned int nocaps;
        }       pcmfile;


// =======================================
void message_alert(char *astring, char *title)
{
  fprintf(pcmlogfile,"%s\n",astring);
}
// ===================================
void initialize_pcmodel(char *mmxfilename)
{
  int i,field;
    user.dielec = FALSE;
    units.dielec = 1.0;
    zero_data();
    read_datafiles(mmxfilename);
    set_field(MMX);
    // initialize atom definitions
    atom_def.natomtype = atom_k.natomtype;
    for (i=1; i < MAXATOMTYPE; i++)
      {
        atom_def.type[i] = atom_k.type[i];
        atom_def.valency[i] = atom_k.valency[i];
        atom_def.number[i] = atom_k.number[i];
        atom_def.ligands[i] = atom_k.ligands[i];
        atom_def.weight[i] = atom_k.weight[i];
        strcpy(atom_def.symbol[i],atom_k.symbol[i]);
      }
    field = get_field();
    Openbox.ftype = FTYPE_PCM;
    Savebox.ftype = FTYPE_PCM;
    strcpy(Openbox.path,pcwindir);
    strcpy(Savebox.path,pcwindir);
    reset_calc_parameters();
    natom = 0;
    /*  minimizer control */
    minim_control.method = 3;
    minim_control.field = MMX;
    minim_control.added_const = FALSE;
    // printout
    minim_values.iprint = FALSE;
    /*  cutoffs    */
    cutoffs.vdwcut = 10.0;
    cutoffs.pmecut = 9.0;
    cutoffs.chrgcut = 100.0;
    cutoffs.dipcut = 8.0;
    // fixed stuff
    fx_atom.natom_fix = 0;
    fx_dist.ndfix = 0;
    fx_angle.nafix = 0;
    fx_torsion.ntfix = 0;
    restrain_atom.natom_restrain = 0;
    // job control
    job_control.use_charge = FALSE;
    job_control.use_scale_charge = FALSE;
    job_control.use_gbsa = FALSE;
    job_control.scale = 1.0;
    job_control.monitor = FALSE;
    job_control.interval = 0.0;
}

int initialize(void)
{
  natom = 0;
  reset_calc_parameters();
    //    reset_atom_data();
    return(0);
}
/* ============================================== */
void reset_calc_parameters(void)
{

/* default to hydrogen bonding on  */
   minim_values.ndc = 4;
   minim_values.nconst= 0;
}
/* ============================================== */       
void reset_atom_data(void)
{

}
