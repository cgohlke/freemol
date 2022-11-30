/* NOTICE: this source code file has been modified for use with FreeMOL */

// File Information
#define FTYPE_PCM                   109
#define FTYPE_SDF                   123


/*  PCMODEL specific definitions */

/* EXTERN struct t_atom {
  int type[MAXATOM],tclass[MAXATOM],mmx_type[MAXATOM], mm3_type[MAXATOM], mmff_type[MAXATOM], atomnum[MAXATOM], use[MAXATOM];
  int iat[MAXATOM][MAXIAT],bo[MAXATOM][MAXIAT];
  long int flags[MAXATOM];
  double x[MAXATOM], y[MAXATOM], z[MAXATOM], atomwt[MAXATOM];
  double  charge[MAXATOM], formal_charge[MAXATOM], sigma_charge[MAXATOM], radius[MAXATOM];
  char name[MAXATOM][3];
  } atom;  */

typedef char LABEL[3];

EXTERN struct t_atom {
  int *type, *tclass, *mmx_type, *mmff_type, *gaff_type, *atomnum, *use;
  int *nconnect,*tbo;
  int *input_charge; // same as formal charge on input
  int **iat, **bo;
  long int *flags;
  double *x, *y, *z, *atomwt;
  double *charge, *formal_charge, *sigma_charge, *radius;
  LABEL *name;
  } atom;


EXTERN int              natom;
EXTERN FILE             *pcmoutfile;
EXTERN char             pcwindir[80];
EXTERN int              **skip;

/* ##define MAXATOM         1000
##define MAXBND          6*MAXATOM/5
##define MAXANG          12*MAXATOM/5
##define MAXTOR          4*MAXATOM */
