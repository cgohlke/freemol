/* NOTICE: this source code file has been modified for use with FreeMOL */
#define EXTERN extern

#include "pcwin.h"
#include "pcmod.h"
#include "atom_k.h"

FILE * fopen_path ( char * , char * , char * ) ;
void zero_data(void);
void readpid(char *);
void read_datafiles(char *);
void readmmxdata(char *);
void torphase(int, float *, int *);
void read_parameterfile(char **);
void numeral(int, char *, int);
void message_alert(char *, char *);
int FetchRecord(FILE *, char *);
void read_added_const(void);
void get_added_const(void);
int check_dipoleconst(int,int,float);
int check_torsionconst(int,int,int,int,float,float,float,int,int,int);
int check_torsion4const(int,int,int,int,float,float,float,int,int,int);
int check_torsion5const(int,int,int,int,float,float,float,int,int,int);
int check_angleconst(int,int,int,float,float,float,float);
int check_angle5const(int,int,int,float,float,float,float);
int check_angle4const(int,int,int,float,float,float,float);
int check_angle3const(int,int,int,float,float,float,float);
int check_bondconst(int,int,float,float);
int check_bond5const(int,int,float,float);
int check_bond4const(int,int,float,float);
int check_bond3const(int,int,float,float);
int check_opbend(int,int,int,int,float);
int check_strbnd(int,int,int,float,float,float);
int check_vdwpr(int,int,float,float);
void set_field(int type);
void set_field_name(char *);
void set_field_bondunit(float ftemp);
void set_field_bondcubic(float ftemp);
void set_field_bondquartic(float ftemp);
void set_field_angleunit(float ftemp);
void set_field_anglecubic(float ftemp);
void set_field_anglequartic(float ftemp);
void set_field_anglepentic(float ftemp);
void set_field_anglesextic(float ftemp);
void set_field_strbndunit(float ftemp);
void set_field_angangunit(float ftemp);
void set_field_strtorunit(float ftemp);
void set_field_torsionunit(float ftemp);
void set_field_vdwtype(char *name);
void set_field_radiustype(char *name);
void set_field_radiussize(char *name);
void set_field_radiusrule(char *name);
void set_field_epsrule(char *name);
void set_field_aterm(float ftemp);
void set_field_bterm(float ftemp);
void set_field_cterm(float ftemp);
void set_field_vdwscale(float ftemp);
void set_field_chrgscale(float ftemp);
void set_field_dielectric(float ftemp);

EXTERN struct t_minim_control {
        int type, method, field, added_const;
        char added_path[256],added_name[256];
        } minim_control;
EXTERN struct t_minim_values {
        int iprint, ndc, nconst;
        float dielc;
        } minim_values;

struct  t_bondk1 {
        int use_ring3, use_ring4, use_ring5;
        int nbnd, nbnd3, nbnd4, nbnd5, ndeloc;
        char kb[MAXBONDCONST][7], kb3[MAXBOND3CONST][7],kb4[MAXBOND4CONST][7],kb5[MAXBOND5CONST][7];
        char kbdel[MAXBONDDELOC][7];
        float s[MAXBONDCONST], t[MAXBONDCONST];
        float s3[MAXBOND3CONST], t3[MAXBOND3CONST];
        float s4[MAXBOND4CONST], t4[MAXBOND4CONST];
        float s5[MAXBOND5CONST], t5[MAXBOND5CONST];
        float sdel[MAXBONDDELOC], tdel[MAXBONDDELOC];
        }  bondk1;

struct t_electroneg {
        int nelecti, nelectj;
        int itype[200],ibond[200],iattach[200];
        int jtype[200],jbond[200],jattach[200];              
        float icorr[200],jcorr[200];
        }  electroneg;
                
struct t_angk1 {
        int use_ang3, use_ang4, use_ang5;
        int nang, nang3, nang4, nang5;
        int ndel, ndel3, ndel4;
        char  ktype[MAXANGCONST][10],ktype3[MAXANG3CONST][10],ktype4[MAXANG4CONST][10],ktype5[MAXANG5CONST][10];
        char  kdel[MAXANGDEL][10],kdel3[MAXANG3DEL][10],kdel4[MAXANG4DEL][10];
        int index[MAXANGCONST],index3[MAXANG3CONST],index4[MAXANG4CONST],index5[MAXANG5CONST];
        int indexdel[MAXANGDEL],indexdel3[MAXANG3DEL],indexdel4[MAXANG4DEL];
        float con[MAXANGCONST], ang[MAXANGCONST][3];
        float con3[MAXANG3CONST], ang3[MAXANG3CONST][3];        
        float con4[MAXANG4CONST], ang4[MAXANG4CONST][3];        
        float con5[MAXANG5CONST], ang5[MAXANG5CONST][3];
        float condel[MAXANGDEL],  angdel[MAXANGDEL][3];      
        float condel3[MAXANG3DEL], angdel3[MAXANG3DEL][3];      
        float condel4[MAXANG4DEL], angdel4[MAXANG4DEL][3];      
         } angk1;
        
struct t_angf {
        int use_angf, nfang;
        char kftype[MAXANGCONST][10];
        float fcon[MAXANGCONST],fc0[MAXANGCONST],fc1[MAXANGCONST],fc2[MAXANGCONST];
        float fc3[MAXANGCONST],fc4[MAXANGCONST],fc5[MAXANGCONST],fc6[MAXANGCONST];
        } angf;
        
struct t_ureybrad_k {
         int nurey_brad;
         char kubang[MAXUREY][10];
         float ubconst[MAXUREY],ubdist[MAXUREY];
         } ureybrad_k;
         
struct t_crossterm_k {
        int nangang, nstrbnd, nstrtor;
        int ang_ang[MAXAA], stbnindex[MAXSTBN] ;
        char str_tor[MAXSTRTOR][7], stbn[MAXSTBN][10];
        float aacon[MAXAA][3], stbncon[MAXSTBN][3], str_torcon[MAXSTRTOR];
        } crossterm_k;
        
struct t_ooplane_k {
         int nopbend;
         char iopb[MAXOOP][13];
         float copb[MAXOOP];
        } ooplane_k;

struct t_improptor_k {
        int nimptor, nimprop;
        char  kv[MAXIMP][13];
        float cimptor[MAXIMP], tdi[MAXIMP];
        float v1[MAXIMP], v2[MAXIMP], v3[MAXIMP];
        int   ph1[MAXIMP], ph2[MAXIMP], ph3[MAXIMP];
        } improptor_k;

struct t_torkn1 {
        int  use_tor4, use_tor5;
        int  ntor, ntor4, ntor5, ntordel, torindex[MAXTORDEL];
        char  kv[MAXTORCONST][13], kv4[MAXTOR4CONST][13], kv5[MAXTOR5CONST][13];
        char  kvdel[MAXTORDEL][13];
        float tv1[MAXTORCONST], tv2[MAXTORCONST], tv3[MAXTORCONST];
        float tv4[MAXTORCONST], tv5[MAXTORCONST], tv6[MAXTORCONST];
        int   phase1[MAXTORCONST], phase2[MAXTORCONST], phase3[MAXTORCONST];
        int   phase4[MAXTORCONST], phase5[MAXTORCONST], phase6[MAXTORCONST];
        float tv41[MAXTOR4CONST], tv42[MAXTOR4CONST], tv43[MAXTOR4CONST];
        int   phase41[MAXTOR4CONST], phase42[MAXTOR4CONST], phase43[MAXTOR4CONST];
        float tv51[MAXTOR5CONST], tv52[MAXTOR5CONST], tv53[MAXTOR5CONST];
        int   phase51[MAXTOR5CONST], phase52[MAXTOR5CONST], phase53[MAXTOR5CONST];
        float tvdel1[MAXTORDEL], tvdel2[MAXTORDEL], tvdel3[MAXTORDEL];
        int   phasedel1[MAXTORDEL], phasedel2[MAXTORDEL], phasedel3[MAXTORDEL];
        } torkn1;

struct t_vdw1 {
        int  nvdw;
        float rad[MAXVDWCONST], eps[MAXVDWCONST];
        int lpd[MAXVDWCONST], ihtyp[MAXVDWCONST], ihdon[MAXVDWCONST];
        float alpha[MAXVDWCONST],n[MAXVDWCONST],a[MAXVDWCONST],g[MAXVDWCONST];
        char da[MAXVDWCONST][2];
        } vdw1;

struct t_vdwpr_k {
        int nvdwpr;
        int  ia1[MAXBONDCONST],ia2[MAXBONDCONST];
        char kv[MAXBONDCONST][7];
        float radius[MAXBONDCONST], eps[MAXBONDCONST];
        } vdwpr_k;
        
struct t_charge_k {
         int ncharge, nbndchrg, nbndchrgdel;
         int type[MAXATOMTYPE], btype[MAXBONDCONST], btypedel[MAXBONDCONST];
         float charge[MAXATOMTYPE], bcharge[MAXBONDCONST], formchrg[MAXATOMTYPE], bchargedel[MAXBONDCONST];
         float typechrg[MAXATOMTYPE];
         } charge_k;
         
struct  t_dipole_k {
        int ndipole,ndipole3,ndipole4,ndipole5;
        char   kb[MAXBONDCONST][7], kb3[MAXBOND3CONST][7], kb4[MAXBOND4CONST][7], kb5[MAXBOND5CONST][7];
        float bmom[MAXBONDCONST],bmom3[MAXBOND3CONST],bmom4[MAXBOND4CONST],bmom5[MAXBOND5CONST];
         } dipole_k;

struct t_metaldata {
        int nmetal, type[200];
        char name[200][3];
        float radius[200], eps[200];
        }  metaldata;

EXTERN struct t_user {
        int dielec;
        } user;

void zero_data()
{
      int ii,j;

      bondk1.use_ring3 = FALSE;
      bondk1.use_ring4 = FALSE;
      bondk1.use_ring5 = FALSE;

      angk1.use_ang3 = FALSE;
      angk1.use_ang4 = FALSE;
      angk1.use_ang5 = FALSE;
      angf.use_angf = FALSE;

      torkn1.use_tor4 = FALSE;
      torkn1.use_tor5 = FALSE;

      units.bndunit = 1.0;
      units.cbnd = 0.0;
      units.qbnd = 0.0;
      units.angunit = 1.0 / (radian*radian);
      units.cang = 0.0;
      units.qang = 0.0;
      units.pang = 0.0;
      units.sang = 0.0;
      units.aaunit = 1.0 / (radian*radian);
      units.stbnunit = 1.0;
      units.ureyunit = 1.0;
      units.torsunit = 1.0;
      units.storunit = 1.0;
      units.v14scale = 1.0;
      units.aterm = 0.0;
      units.bterm = 0.0;
      units.cterm = 0.0;
//      units.dielec = 1.0;
      units.chgscale = 1.0;
//  atom type constants
      atom_k.natomtype = 0;
      charge_k.ncharge = 0;
      charge_k.nbndchrg = 0;
      charge_k.nbndchrgdel = 0;
      for (ii = 0; ii < MAXATOMTYPE; ii ++)
      {
          atom_k.type[ii]= 0;
          atom_k.valency[ii]= 0;
          atom_k.tclass[ii] = 0;
          atom_k.tclass1[ii] = 0;
          atom_k.tclass2[ii] = 0;
          atom_k.number[ii] = 0;
          atom_k.ligands[ii] = 0;
          atom_k.weight[ii] = 0.00F;
          strcpy(atom_k.symbol[ii],"    ");
          strcpy(atom_k.description[ii],"                  ");
          charge_k.type[ii] = 0;
          charge_k.charge[ii] = 0.00F;
          charge_k.formchrg[ii] = 0.00F;
          charge_k.typechrg[ii] = 0.00F;
      }
//  metals
      metaldata.nmetal = 0;
      for (ii = 0; ii < 200; ii++)
      {
          metaldata.type[ii] = 0;
          metaldata.radius[ii] = 0.0;
          metaldata.eps[ii] = 0.0;
          strcpy(metaldata.name[ii],"  ");
      }
//  bond constant data
      bondk1.nbnd = 0;
      vdwpr_k.nvdwpr = 0;
      dipole_k.ndipole = 0;
      electroneg.nelecti = 0;
      electroneg.nelectj = 0;       

      for( ii = 0; ii < MAXBONDCONST; ii++ )
      {
          strcpy(bondk1.kb[ii],"      ");
          bondk1.s[ii] = 0.;
          bondk1.t[ii] = 0.;
          strcpy(vdwpr_k.kv[ii],"      ");
          vdwpr_k.ia1[ii] = 0;
          vdwpr_k.ia2[ii] = 0;
          vdwpr_k.radius[ii] = 0.0;
          vdwpr_k.eps[ii] = 0.0;
          strcpy(dipole_k.kb[ii],"      ");
          dipole_k.bmom[ii] = 0.;
          charge_k.btype[ii] = 0;
          charge_k.bcharge[ii] = 0.0;
          charge_k.btypedel[ii] = 0;
          charge_k.bchargedel[ii] = 0.0;
      }
//  bond 3 constant
      bondk1.nbnd3 = 0;
      dipole_k.ndipole3 = 0;
      for( ii = 0; ii < MAXBOND3CONST; ii++ )
      {
          strcpy(bondk1.kb3[ii],"      ");
          bondk1.s3[ii] = 0.;
          bondk1.t3[ii] = 0.;
          strcpy(dipole_k.kb3[ii],"      ");
          dipole_k.bmom3[ii] = 0.;
      }
// bond 4 constant
      bondk1.nbnd4 = 0;
      dipole_k.ndipole4 = 0;
      for( ii = 0; ii < MAXBOND4CONST; ii++ )
      {
          strcpy(bondk1.kb4[ii],"      ");
          bondk1.s4[ii] = 0.;
          bondk1.t4[ii] = 0.;
          strcpy(dipole_k.kb4[ii],"      ");
          dipole_k.bmom4[ii] = 0.;
      }
// bond 5 constant
      bondk1.nbnd5 = 0;
      dipole_k.ndipole5 = 0;
      for( ii = 0; ii < MAXBOND5CONST; ii++ )
      {
          strcpy(bondk1.kb5[ii],"      ");
          bondk1.s5[ii] = 0.;
          bondk1.t5[ii] = 0.;
          strcpy(dipole_k.kb5[ii],"      ");
          dipole_k.bmom5[ii] = 0.;
      }
// deloc bond  constant
      bondk1.ndeloc = 0;
      for( ii = 0; ii < MAXBONDDELOC; ii++ )
      {
          strcpy(bondk1.kbdel[ii],"      ");
          bondk1.sdel[ii] = 0.;
          bondk1.tdel[ii] = 0.;
      }
// VDW constants
      vdw1.nvdw = 0;
      for( ii = 0; ii < MAXVDWCONST; ii++ )
      {
          vdw1.rad[ii] = 0.;
          vdw1.eps[ii] = 0.;
          vdw1.lpd[ii] = 0;
          vdw1.ihtyp[ii] = 0;
          vdw1.ihdon[ii] = 0;
          vdw1.alpha[ii] = 0.0;
          vdw1.n[ii] = 0.0;
          vdw1.a[ii] = 0.0;
          vdw1.g[ii] = 0.0;
          strcpy(vdw1.da[ii]," ");
      }
//  torsion
      torkn1.ntor = 0;
      for( ii = 0; ii < MAXTORCONST; ii++ )
      {
          strcpy(torkn1.kv[ii],"            ");
          torkn1.tv1[ii] = 0.;
          torkn1.tv2[ii] = 0.;
          torkn1.tv3[ii] = 0.;
          torkn1.tv4[ii] = 0.;
          torkn1.tv5[ii] = 0.;
          torkn1.tv6[ii] = 0.;
          torkn1.phase1[ii] = torkn1.phase2[ii] = torkn1.phase3[ii] = 0;
          torkn1.phase4[ii] = torkn1.phase5[ii] = torkn1.phase6[ii] = 0;
      }
//  torsion 4
      torkn1.ntor4 = 0;
      for( ii = 0; ii < MAXTOR4CONST; ii++ )
      {
          strcpy(torkn1.kv4[ii],"            ");
          torkn1.tv41[ii] = 0.;
          torkn1.tv42[ii] = 0.;
          torkn1.tv43[ii] = 0.;
          torkn1.phase41[ii] = torkn1.phase42[ii] = torkn1.phase43[ii] = 0;
      }
//  torsion 5
      torkn1.ntor5 = 0;
      for( ii = 0; ii < MAXTOR5CONST; ii++ )
      {
          strcpy(torkn1.kv5[ii],"            ");
          torkn1.tv51[ii] = 0.;
          torkn1.tv52[ii] = 0.;
          torkn1.tv53[ii] = 0.;
          torkn1.phase51[ii] = torkn1.phase52[ii] = torkn1.phase53[ii] = 0;
      }
//  torsion delocalized
      torkn1.ntordel = 0;
      for( ii = 0; ii < MAXTORDEL; ii++ )
      {
          strcpy(torkn1.kvdel[ii],"            ");
          torkn1.torindex[ii] = 0;
          torkn1.tvdel1[ii] = 0.;
          torkn1.tvdel2[ii] = 0.;
          torkn1.tvdel3[ii] = 0.;
          torkn1.phasedel1[ii] = torkn1.phasedel2[ii] = torkn1.phasedel3[ii] = 0;
      }
// fourier angles
      angf.nfang = 0;
      for (ii=0; ii < MAXANGCONST; ii++ )
      {
          strcpy(angf.kftype[ii],"         ");
          angf.fcon[ii] = 0.0;
          angf.fc0[ii] = 0.0;
          angf.fc1[ii] = 0.0;
          angf.fc2[ii] = 0.0;
          angf.fc3[ii] = 0.0;
          angf.fc4[ii] = 0.0;
          angf.fc5[ii] = 0.0;
          angf.fc6[ii] = 0.0;
      }
//  angle
      angk1.nang = 0;
      for( ii = 0; ii < MAXANGCONST; ii++ )
      {
          strcpy(angk1.ktype[ii],"         ");
	  angk1.index[ii] = 0;
          angk1.con[ii] = 0.;
          for( j = 0; j < 3; j++ )
              angk1.ang[ii][j] = 0.;
      }
//  angle 3
      angk1.nang3 = 0;
      for( ii = 0; ii < MAXANG3CONST; ii++ )
      {
          strcpy(angk1.ktype3[ii],"         ");
	  angk1.index3[ii] = 0;
          angk1.con3[ii] = 0.;
          for( j = 0; j < 3; j++ )
              angk1.ang3[ii][j] = 0.;
      }
//  angle 4
      angk1.nang4 = 0;
      for( ii = 0; ii < MAXANG4CONST; ii++ )
      {
          strcpy(angk1.ktype4[ii],"         ");
	  angk1.index4[ii] = 0;
          angk1.con4[ii] = 0.;
          for( j = 0; j < 3; j++ )
              angk1.ang4[ii][j] = 0.;
      }
//  angle 5
      angk1.nang5 = 0;
      for( ii = 0; ii < MAXANG5CONST; ii++ )
      {
          strcpy(angk1.ktype5[ii],"         ");
	  angk1.index5[ii] = 0;
          angk1.con5[ii] = 0.;
          for( j = 0; j < 3; j++ )
              angk1.ang5[ii][j] = 0.;
      }
//  angle delocalized
      angk1.ndel = 0;
      for( ii = 0; ii < MAXANGDEL; ii++ )
      {
          strcpy(angk1.kdel[ii],"         ");
	  angk1.indexdel[ii] = 0;
          angk1.condel[ii] = 0.;
          for( j = 0; j < 3; j++ )
              angk1.angdel[ii][j] = 0.;
      }
//  angle3 delocalized
      angk1.ndel3 = 0;
      for( ii = 0; ii < MAXANG3DEL; ii++ )
      {
          strcpy(angk1.kdel3[ii],"         ");
	  angk1.indexdel3[ii] = 0;
          angk1.condel3[ii] = 0.;
          for( j = 0; j < 3; j++ )
              angk1.angdel3[ii][j] = 0.;
      }
//  angle4 delocalized
      angk1.ndel4 = 0;
      for( ii = 0; ii < MAXANG4DEL; ii++ )
      {
          strcpy(angk1.kdel4[ii],"         ");
	  angk1.indexdel4[ii] = 0;
          angk1.condel4[ii] = 0.;
          for( j = 0; j < 3; j++ )
              angk1.angdel4[ii][j] = 0.;
      }
// urey bradley
      ureybrad_k.nurey_brad = 0;
      for (ii = 0; ii < MAXUREY; ii++)
      {
          strcpy(ureybrad_k.kubang[ii],"         ");
          ureybrad_k.ubconst[ii] = 0.0;
          ureybrad_k.ubdist[ii] = 0.0;
      }
//  out of plane
      ooplane_k.nopbend = 0;
      for( ii = 0; ii < MAXOOP; ii++ )
      {
          strcpy(ooplane_k.iopb[ii],"            ");
          ooplane_k.copb[ii] = 0.;
      }
// improper torsions
      improptor_k.nimptor = 0;
      improptor_k.nimprop = 0;
      for( ii = 0; ii < MAXIMP; ii++ )
      {
          strcpy(improptor_k.kv[ii],"            ");
          improptor_k.cimptor[ii] = 0.;
          improptor_k.tdi[ii] = 0.;
          improptor_k.v1[ii] = 0.0;
          improptor_k.v2[ii] = 0.0;
          improptor_k.v3[ii] = 0.0;
          improptor_k.ph1[ii] = 0;
          improptor_k.ph2[ii] = 0;
          improptor_k.ph3[ii] = 0;
      }
// cross terms
        crossterm_k.nangang = 0;
        crossterm_k.nstrbnd = 0;
        crossterm_k.nstrtor = 0;
        for (ii = 0; ii < MAXAA; ii++)
        {
            crossterm_k.ang_ang[ii] = 0;
            crossterm_k.aacon[ii][0] = 0.0;
            crossterm_k.aacon[ii][1] = 0.0;
            crossterm_k.aacon[ii][2] = 0.0;
        }
        for (ii = 0; ii < MAXSTBN; ii++)
        {
            strcpy(crossterm_k.stbn[ii],"         ");
            crossterm_k.stbnindex[ii] = 0;
            crossterm_k.stbncon[ii][0] = 0.0;
            crossterm_k.stbncon[ii][1] = 0.0;
            crossterm_k.stbncon[ii][2] = 0.0;
        }
        for (ii = 0; ii < MAXSTRTOR; ii++)
        {
            strcpy(crossterm_k.str_tor[ii],"      ");
            crossterm_k.str_torcon[ii] = 0.0;
        }        
}
/* =========================================== */
void get_added_const()
{
    int ia,ib,ic,id, ip1,ip2,ip3,iz,if3,if4,if5;
    float f1,f2,f3,f4;
    float v1[6];
    int se1[6];
    char pa[3],pb[3],pc[3],pd[3],pt[13];
    char line[151], iptemp[21], dummy[10];
    FILE *wfile;

    wfile = fopen_path(minim_control.added_path,minim_control.added_name,"r");
    if (wfile == NULL)
    {
        message_alert("Error opening added constants file","Error");
        fprintf(pcmlogfile,"Error reading added constants file: %s\n",minim_control.added_name);
        return;
    }

    fprintf(pcmlogfile,"\nThe Following Added Parameters were read:\n\n");
    while ( FetchRecord(wfile,line) )
    {
        sscanf(line,"%s",iptemp);

        if (strcmp(iptemp,"bond") == 0)
        {
            sscanf(line,"%s %d %d %f %f",dummy,&ia,&ib,&f1,&f2);
            fprintf(pcmlogfile,"Bond: AtomTypes %d %d   Blen   %f    Bk %f \n",ia,ib,f1,f2);
            if (ia < MAXATOMTYPE && ib < MAXATOMTYPE)
            {
              iz = check_bondconst(ia,ib,f1,f2);
              if (iz == FALSE)
              {
                bondk1.s[bondk1.nbnd] = f1;
                bondk1.t[bondk1.nbnd] = f2;
                numeral(ia,pa,3);
                numeral(ib,pb,3);
                strcpy(pt,pa);
                strcat(pt,pb);
                strcpy(bondk1.kb[bondk1.nbnd],pt);
                bondk1.nbnd++;
              }
            } else
            {
                message_alert("Error reading bond atom types in added constants file","Error");
            }
        } else if (strcmp(iptemp,"bond5") == 0)
        {
            iz = sscanf(line,"%s %d %d %f %f",dummy,&ia,&ib,&f1,&f2);
            fprintf(pcmlogfile,"Bond5: AtomTypes %d %d   Blen   %f    Bk %f \n",ia,ib,f1,f2);
            if (ia < MAXATOMTYPE && ib < MAXATOMTYPE && iz > 3)
            {
              iz = check_bond5const(ia,ib,f1,f2);
              if (iz == FALSE)
              {
                bondk1.s5[bondk1.nbnd5] = f1;
                bondk1.t5[bondk1.nbnd5] = f2;
                numeral(ia,pa,3);
                numeral(ib,pb,3);
                strcpy(pt,pa);
                strcat(pt,pb);
                strcpy(bondk1.kb5[bondk1.nbnd5],pt);
                bondk1.nbnd5++;
              }
            } else
            {
                message_alert("Error reading bond atom types in added constants file","Error");
            }
        } else if (strcmp(iptemp,"bond4") == 0)
        {
            iz = sscanf(line,"%s %d %d %f %f",dummy,&ia,&ib,&f1,&f2);
            fprintf(pcmlogfile,"Bond4: AtomTypes %d %d   Blen   %f    Bk %f \n",ia,ib,f1,f2);
            if (ia < MAXATOMTYPE && ib < MAXATOMTYPE && iz > 3)
            {
              iz = check_bond4const(ia,ib,f1,f2);
              if (iz == FALSE)
              {
                bondk1.s4[bondk1.nbnd4] = f1;
                bondk1.t4[bondk1.nbnd4] = f2;
                numeral(ia,pa,3);
                numeral(ib,pb,3);
                strcpy(pt,pa);
                strcat(pt,pb);
                strcpy(bondk1.kb4[bondk1.nbnd4],pt);
                bondk1.nbnd4++;
              }
            } else
            {
                message_alert("Error reading bond atom types in added constants file","Error");
            }
        } else if (strcmp(iptemp,"bond3") == 0)
        {
            sscanf(line,"%s %d %d %f %f",dummy,&ia,&ib,&f1,&f2);
            fprintf(pcmlogfile,"Bond3: AtomTypes %d %d   Blen   %f    Bk %f \n",ia,ib,f1,f2);
            if (ia < MAXATOMTYPE && ib < MAXATOMTYPE)
            {
              iz = check_bond3const(ia,ib,f1,f2);
              if (iz == FALSE)
              {
                bondk1.s3[bondk1.nbnd3] = f1;
                bondk1.t3[bondk1.nbnd3] = f2;
                numeral(ia,pa,3);
                numeral(ib,pb,3);
                strcpy(pt,pa);
                strcat(pt,pb);
                strcpy(bondk1.kb3[bondk1.nbnd3],pt);
                bondk1.nbnd3++;
              }
            } else
            {
                message_alert("Error reading bond atom types in added constants file","Error");
            }
        } else if (strcmp(iptemp,"angle") == 0)
        {
            f1 =f2=f3=f4=0.0;
            sscanf(line,"%s %d %d %d %f %f %f %f",dummy,&ia,&ib,&ic,&f1,&f2,&f3,&f4);
            fprintf(pcmlogfile,"Angle: AtomTypes %d %d %d  Anat   %f    Acon %f \n",ia,ib,ic,f2,f1);
            if (ia < MAXATOMTYPE && ib < MAXATOMTYPE && ic < MAXATOMTYPE)
            {
              iz = check_angleconst(ia,ib,ic,f1,f2,f3,f4);
              if (iz == FALSE)
              {
               angk1.con[angk1.nang] = f1;
               angk1.ang[angk1.nang][0] = f2; 
               angk1.ang[angk1.nang][1] = f3; 
               angk1.ang[angk1.nang][2] = f4; 
               numeral(ia,pa,3);
               numeral(ib,pb,3);
               numeral(ic,pc,3);
               strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc);
               strcpy(angk1.ktype[angk1.nang],pt);
               angk1.nang++;
              }
            } else
            {
                message_alert("Error reading Angle atom types in added constants file","Error");
            }
        } else if (strcmp(iptemp,"angle5") == 0)
        {
            f1 =f2=f3=f4=0.0;
            sscanf(line,"%s %d %d %d %f %f %f %f",dummy,&ia,&ib,&ic,&f1,&f2,&f3,&f4);
            fprintf(pcmlogfile,"Angle5: AtomTypes %d %d %d  Anat   %f    Acon %f \n",ia,ib,ic,f2,f1);
            if (ia < MAXATOMTYPE && ib < MAXATOMTYPE && ic < MAXATOMTYPE)
            {
              iz = check_angle5const(ia,ib,ic,f1,f2,f3,f4);
              if (iz == FALSE)
              {
               angk1.con5[angk1.nang5] = f1;
               angk1.ang5[angk1.nang5][0] = f2; 
               angk1.ang5[angk1.nang5][1] = f3; 
               angk1.ang5[angk1.nang5][2] = f4; 
               numeral(ia,pa,3);
               numeral(ib,pb,3);
               numeral(ic,pc,3);
               strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc);
               strcpy(angk1.ktype5[angk1.nang5],pt);
               angk1.nang5++;
              }
            } else
            {
                message_alert("Error reading Angle atom types in added constants file","Error");
            }
        } else if (strcmp(iptemp,"angle4") == 0)
        {
            f1 =f2=f3=f4=0.0;
            sscanf(line,"%s %d %d %d %f %f %f %f",dummy,&ia,&ib,&ic,&f1,&f2,&f3,&f4);
            fprintf(pcmlogfile,"Angle4: AtomTypes %d %d %d  Anat   %f    Acon %f \n",ia,ib,ic,f2,f1);
            if (ia < MAXATOMTYPE && ib < MAXATOMTYPE && ic < MAXATOMTYPE)
            {
              iz = check_angle4const(ia,ib,ic,f1,f2,f3,f4);
              if (iz == FALSE)
              {
               angk1.con4[angk1.nang4] = f1;
               angk1.ang4[angk1.nang4][0] = f2; 
               angk1.ang4[angk1.nang4][1] = f3; 
               angk1.ang4[angk1.nang4][2] = f4; 
               numeral(ia,pa,3);
               numeral(ib,pb,3);
               numeral(ic,pc,3);
               strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc);
               strcpy(angk1.ktype4[angk1.nang4],pt);
               angk1.nang4++;
              }
            } else
            {
                message_alert("Error reading Angle atom types in added constants file","Error");
            }
        } else if (strcmp(iptemp,"angle3") == 0)
        {
            f1 =f2=f3=f4=0.0;
            sscanf(line,"%s %d %d %d %f %f %f %f",dummy,&ia,&ib,&ic,&f1,&f2,&f3,&f4);
            fprintf(pcmlogfile,"Angle3: AtomTypes %d %d %d  Anat   %f    Acon %f \n",ia,ib,ic,f2,f1);
            if (ia < MAXATOMTYPE && ib < MAXATOMTYPE && ic < MAXATOMTYPE)
            {
              iz = check_angle3const(ia,ib,ic,f1,f2,f3,f4);
              if (iz == FALSE)
              {
               angk1.con3[angk1.nang3] = f1;
               angk1.ang3[angk1.nang3][0] = f2; 
               angk1.ang3[angk1.nang3][1] = f3; 
               angk1.ang3[angk1.nang3][2] = f4; 
               numeral(ia,pa,3);
               numeral(ib,pb,3);
               numeral(ic,pc,3);
               strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc);
               strcpy(angk1.ktype3[angk1.nang3],pt);
               angk1.nang3++;
              }
            } else
            {
                message_alert("Error reading Angle atom types in added constants file","Error");
            }
        } else if (strcmp(iptemp,"torsion") == 0)
        {
            sscanf(line,"%s %d %d %d %d %f %d %f %d %f %d",dummy,&ia,&ib,&ic,&id,&f1,&ip1,&f2,
                &ip2,&f3,&ip3);
            fprintf(pcmlogfile,"Torsion: AtomTypes %d %d %d %d  V1   %f    V2 %f    V3 %f\n",ia,ib,ic,id,f1,f2,f3);
            if (ia < MAXATOMTYPE && ib < MAXATOMTYPE && ic < MAXATOMTYPE && id < MAXATOMTYPE)
            {
              iz = check_torsionconst(ia,ib,ic,id,f1,f2,f3,ip1,ip2,ip3);
              if (iz == FALSE)
              {
                v1[0] = f1; v1[1] = f2; v1[2] = f3;
                v1[3] = v1[4] = v1[5] = 0.0;
                se1[0] = ip1; se1[1] = ip2; se1[2] = ip3;
                se1[3] = se1[4] = se1[5] = 0;
                torphase(6, v1, se1);
                torkn1.tv1[torkn1.ntor] = v1[0];
                torkn1.phase1[torkn1.ntor] = se1[0];
                torkn1.tv2[torkn1.ntor] = v1[1];
                torkn1.phase2[torkn1.ntor] = se1[1];
                torkn1.tv3[torkn1.ntor] = v1[2];
                torkn1.phase3[torkn1.ntor] = se1[2];
                torkn1.tv4[torkn1.ntor] = v1[3];
                torkn1.phase4[torkn1.ntor] = se1[3];
                torkn1.tv5[torkn1.ntor] = v1[4];
                torkn1.phase5[torkn1.ntor] = se1[4];
                torkn1.tv6[torkn1.ntor] = v1[5];
                torkn1.phase6[torkn1.ntor] = se1[5];
                 numeral(ia,pa,3);
                 numeral(ib,pb,3);
                 numeral(ic,pc,3);
                 numeral(id,pd,3);
                 strcpy(pt,pa);
                 strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                 strcpy(torkn1.kv[torkn1.ntor],pt);
                 torkn1.ntor++; 
              }
            } else
            {
                message_alert("Error reading Torsion atom types in added constants file","Error");
            }
        } else if (strcmp(iptemp,"torsion4") == 0)
        {
            sscanf(line,"%s %d %d %d %d %f %d %f %d %f %d",dummy,&ia,&ib,&ic,&id,&f1,&ip1,&f2,
                &ip2,&f3,&ip3);
            fprintf(pcmlogfile,"Torsion4: AtomTypes %d %d %d %d  V1   %f    V2 %f    V3 %f\n",ia,ib,ic,id,f1,f2,f3);
            if (ia < MAXATOMTYPE && ib < MAXATOMTYPE && ic < MAXATOMTYPE && id < MAXATOMTYPE)
            {
              iz = check_torsion4const(ia,ib,ic,id,f1,f2,f3,ip1,ip2,ip3);
              if (iz == FALSE)
              {
                v1[0] = f1; v1[1] = f2; v1[2] = f3;
                v1[3] = v1[4] = v1[5] = 0.0;
                se1[0] = ip1; se1[1] = ip2; se1[2] = ip3;
                se1[3] = se1[4] = se1[5] = 0;
                torphase(6, v1, se1);
                torkn1.tv41[torkn1.ntor4] = v1[0];
                torkn1.phase41[torkn1.ntor4] = se1[0];
                torkn1.tv42[torkn1.ntor4] = v1[1];
                torkn1.phase42[torkn1.ntor4] = se1[1];
                torkn1.tv43[torkn1.ntor4] = v1[2];
                torkn1.phase43[torkn1.ntor4] = se1[2];
                 numeral(ia,pa,3);
                 numeral(ib,pb,3);
                 numeral(ic,pc,3);
                 numeral(id,pd,3);
                 strcpy(pt,pa);
                 strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                 strcpy(torkn1.kv4[torkn1.ntor4],pt);
                 torkn1.ntor4++; 
              }
            } else
            {
                message_alert("Error reading Torsion atom types in added constants file","Error");
            }
        } else if (strcmp(iptemp,"torsion5") == 0)
        {
            sscanf(line,"%s %d %d %d %d %f %d %f %d %f %d",dummy,&ia,&ib,&ic,&id,&f1,&ip1,&f2,
                &ip2,&f3,&ip3);
            fprintf(pcmlogfile,"Torsion5: AtomTypes %d %d %d %d  V1   %f    V2 %f    V3 %f\n",ia,ib,ic,id,f1,f2,f3);
            if (ia < MAXATOMTYPE && ib < MAXATOMTYPE && ic < MAXATOMTYPE && id < MAXATOMTYPE)
            {
              iz = check_torsion5const(ia,ib,ic,id,f1,f2,f3,ip1,ip2,ip3);
              if (iz == FALSE)
              {
                v1[0] = f1; v1[1] = f2; v1[2] = f3;
                v1[3] = v1[4] = v1[5] = 0.0;
                se1[0] = ip1; se1[1] = ip2; se1[2] = ip3;
                se1[3] = se1[4] = se1[5] = 0;
                torphase(6, v1, se1);
                torkn1.tv51[torkn1.ntor5] = v1[0];
                torkn1.phase51[torkn1.ntor5] = se1[0];
                torkn1.tv52[torkn1.ntor5] = v1[1];
                torkn1.phase52[torkn1.ntor5] = se1[1];
                torkn1.tv53[torkn1.ntor5] = v1[2];
                torkn1.phase53[torkn1.ntor5] = se1[2];
                 numeral(ia,pa,3);
                 numeral(ib,pb,3);
                 numeral(ic,pc,3);
                 numeral(id,pd,3);
                 strcpy(pt,pa);
                 strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                 strcpy(torkn1.kv5[torkn1.ntor5],pt);
                 torkn1.ntor5++; 
              }
            } else
            {
                message_alert("Error reading Torsion atom types in added constants file","Error");
            }
        } else if (strcmp(iptemp,"vdw") == 0)
        {
            sscanf(line,"%s %d %f %f %d %d %d",dummy,&ia,&f1,&f2,&if3,&if4,&if5);
            fprintf(pcmlogfile,"VDW: AtomTypes %d   Rad   %f    Eps %f \n",ia,f1,f2);
            vdw1.rad[ia] = f1;
            vdw1.eps[ia] = f2;
            vdw1.lpd[ia] = if3;
            vdw1.ihdon[ia] = if4;
            vdw1.ihtyp[ia] = if5;
        } else if (strcmp(iptemp,"metal") == 0)
        {
            sscanf(line,"%s %d %f %f",dummy,&ia,&f1,&f2);
            fprintf(pcmlogfile,"Metal: AtomTypes %d   Rad   %f    Eps %f \n",ia,f1,f2);
            metaldata.radius[ia] = f1;
            metaldata.eps[ia] = f2;
        } else if (strcmp(iptemp,"vdwpr") == 0)
        {
            sscanf(line,"%s %d %d %f %f",dummy, &ia, &ib, &f1,&f2);
            fprintf(pcmlogfile,"VDWPR: AtomTypes %d %d  Rad   %f    Eps %f \n",ia,ib,f1,f2);
            if (ia < MAXATOMTYPE && ib < MAXATOMTYPE)
            {
              iz = check_vdwpr(ia,ib,f1,f2);
              if (iz == FALSE)
              {
                numeral(ia,pa,3);
                numeral(ib,pb,3);
                strcpy(pt,pa); strcat(pt,pb);
                strcpy(vdwpr_k.kv[vdwpr_k .nvdwpr],pt);
                vdwpr_k.ia1[vdwpr_k .nvdwpr] = ia;
                vdwpr_k.ia2[vdwpr_k .nvdwpr] = ib;
                vdwpr_k.radius[vdwpr_k .nvdwpr] = f1;
                vdwpr_k.eps[vdwpr_k .nvdwpr] = f2;
                vdwpr_k .nvdwpr++;
              }
            } else
            {
                message_alert("Error reading VDWPR atom types in added constants file","Error");
            }
        } else if (strcmp(iptemp,"opbend") == 0)
        {
            sscanf( line, "%s %d %d %d %d %f", dummy, &ia, &ib, &ic,&id,&f1 );
            fprintf(pcmlogfile,"OPBEND: AtomTypes %d %d %d %d  Const  %f \n",ia,ib,ic,id,f1);
            if (ia < MAXATOMTYPE && ib < MAXATOMTYPE && ic < MAXATOMTYPE)
            {
              iz = check_opbend(ia,ib,ic,id,f1);
              if (iz == FALSE)
              {
               numeral(ia,pa,3);
               numeral(ib,pb,3);
               numeral(ic,pc,3);
               numeral(id,pd,3);
               strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
               strcpy(ooplane_k.iopb[ooplane_k.nopbend],pt);
               ooplane_k.nopbend++;
              }
            } else
            {
                message_alert("Error reading Opbend atom types in added constants file","Error");
            }
        } else if (strcmp(iptemp,"strbnd") == 0)
        {
            id = 0;
            sscanf(line,"%s %d %d %d %f %f %f",dummy, &ia, &ib, &ic,
                &f1,&f2,&f3);
            fprintf(pcmlogfile,"STRBEND: AtomTypes %d %d Const  %f %f %f\n",ia,ib,f1,f2,f3);
            if (ia < MAXATOMTYPE && ib < MAXATOMTYPE && ic < MAXATOMTYPE)
            {
              iz = check_strbnd(ia,ib,ic,f1,f2,f3);
              if (iz == FALSE)
              {
                numeral(ia,pa,3);
                numeral(ib,pb,3);
                numeral(ic,pc,3);
                strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc);
                strcpy(crossterm_k.stbn[crossterm_k.nstrbnd],pt);
                crossterm_k.stbncon[crossterm_k.nstrbnd][0] = f1;
                crossterm_k.stbncon[crossterm_k.nstrbnd][1] = f2;
                crossterm_k.stbncon[crossterm_k.nstrbnd][2] = f3;
                crossterm_k.stbnindex[crossterm_k.nstrbnd] = 0;
                crossterm_k.nstrbnd++;
             }
            } else
            {
                message_alert("Error reading Strbend atom types in added constants file","Error");
            }                        
        } else if (strcmp(iptemp,"dipole") == 0)
        {
            sscanf(line,"%s %d %d %f",dummy,&ia,&ib,&f1);
            fprintf(pcmlogfile,"DIPOLE: AtomTypes %d %d  Bond Moment  %f \n",ia,ib,f1);
            if (ia < MAXATOMTYPE && ib < MAXATOMTYPE)
            {
              iz = check_dipoleconst(ia,ib,f1);
              if (iz == FALSE)
              {
                 dipole_k.bmom[dipole_k.ndipole] = f1;
                 numeral(ia,pa,3);
                 numeral(ib,pb,3);
                 strcpy(pt,pa); strcat(pt,pb);
                 strcpy(dipole_k.kb[dipole_k.ndipole],pt);
                 dipole_k.ndipole++;
              }
            } else
            {
                message_alert("Error reading Dipole atom types in added constants file","Error");
            }                        
        }
    }
    fclose(wfile);  
}
/* ========================================== */
int check_vdwpr(int ia, int ib,float f1, float f2)
{
    int i;
    char pa[3],pb[3],pt[13];

    numeral(ia,pa,3);
    numeral(ib,pb,3);
    strcpy(pt,pa);
    strcat(pt,pb);
    for (i=0; i < vdwpr_k .nvdwpr; i++)
    {
        if (strcmp(pt,vdwpr_k.kv[i]) == 0)
        {
            vdwpr_k.radius[i] = f1;
            vdwpr_k.eps[i] = f2;
            return TRUE;
        }
    }
    return FALSE;
}
/* ========================================== */
int check_strbnd(int ia, int ib, int ic,float f1, float f2, float f3)
{
    int i;
    char pa[3],pb[3],pc[3],pt[13];

    numeral(ia,pa,3);
    numeral(ib,pb,3);
    numeral(ic,pc,3);
    strcpy(pt,pa);
    strcat(pt,pb);
    strcat(pt,pc);
    for (i=0; i < crossterm_k.nstrbnd; i++)
    {
        if (strcmp(pt,crossterm_k.stbn[i]) == 0)
        {
            crossterm_k.stbncon[i][0] = f1;
            crossterm_k.stbncon[i][1] = f2;
            crossterm_k.stbncon[i][2] = f3;
            return TRUE;
        }
    }
    return FALSE;
}
/* ========================================== */
int check_opbend(int ia,int ib,int ic,int id, float ftemp)
{
    int i;
    char pa[3],pb[3],pc[3],pd[3],pt[13];
    numeral(ia,pa,3);
    numeral(ib,pb,3);
    numeral(ic,pc,3);
    numeral(id,pd,3);
    strcpy(pt,pa);
    strcat(pt,pb);
    strcat(pt,pc);
    strcat(pt,pd);
    for (i=0; i < ooplane_k.nopbend; i++)
    {
        if (strcmp(pt,ooplane_k.iopb[i]) == 0)
        {
            ooplane_k.copb[i] = ftemp;
            return TRUE;
        }
    }
    return FALSE;
}
/* ========================================== */
int check_dipoleconst(int ia,int ib,float bmom)
{
    int i;
    char pa[3],pb[3],pt[13];
    
    numeral(ia,pa,3);
    numeral(ib,pb,3);
    strcpy(pt,pa); strcat(pt,pb);
    for (i=0; i < dipole_k.ndipole; i++)
    {
        if (strcmp(pt,dipole_k.kb[i]) == 0)
        {        
           dipole_k.bmom[i] = bmom;
           return TRUE;
        }
    }
    return FALSE;
}  
/* ========================================== */
int check_torsionconst(int ia,int ib,int ic,int id,float v1,float v2, float v3,
        int ip1, int ip2, int ip3)
{
    int i;
    char pa[3],pb[3],pc[3],pd[3],pt[13];

    numeral(ia,pa,3);
    numeral(ib,pb,3);
    numeral(ic,pc,3);
    numeral(id,pd,3);
    strcpy(pt,pa);
    strcat(pt,pb);
    strcat(pt,pc);
    strcat(pt,pd);
    for (i= 0; i < torkn1.ntor; i++)
    {
        if (strcmp(pt,torkn1.kv[i]) == 0)
        {
            torkn1.tv1[i] = v1;
            torkn1.tv2[i] = v2;
            torkn1.tv3[i] = v3;
            return TRUE;
        }
    }
    return FALSE;
}
/* ========================================== */
int check_torsion4const(int ia,int ib,int ic,int id,float v1,float v2, float v3,
        int ip1, int ip2, int ip3)
{
    int i;
    char pa[3],pb[3],pc[3],pd[3],pt[13];

    numeral(ia,pa,3);
    numeral(ib,pb,3);
    numeral(ic,pc,3);
    numeral(id,pd,3);
    strcpy(pt,pa);
    strcat(pt,pb);
    strcat(pt,pc);
    strcat(pt,pd);
    for (i= 0; i < torkn1.ntor4; i++)
    {
        if (strcmp(pt,torkn1.kv4[i]) == 0)
        {
            torkn1.tv41[i] = v1;
            torkn1.tv42[i] = v2;
            torkn1.tv43[i] = v3;
            return TRUE;
        }
    }
    return FALSE;
}
/* ========================================== */
int check_torsion5const(int ia,int ib,int ic,int id,float v1,float v2, float v3,
        int ip1, int ip2, int ip3)
{
    int i;
    char pa[3],pb[3],pc[3],pd[3],pt[13];

    numeral(ia,pa,3);
    numeral(ib,pb,3);
    numeral(ic,pc,3);
    numeral(id,pd,3);
    strcpy(pt,pa);
    strcat(pt,pb);
    strcat(pt,pc);
    strcat(pt,pd);
    for (i= 0; i < torkn1.ntor5; i++)
    {
        if (strcmp(pt,torkn1.kv5[i]) == 0)
        {
            torkn1.tv51[i] = v1;
            torkn1.tv52[i] = v2;
            torkn1.tv53[i] = v3;
            return TRUE;
        }
    }
    return FALSE;
}
/* ========================================== */    
int check_angleconst(int ia,int ib,int ic,float f1, float f2,float f3,float f4)
{
    int i;
    char pa[3],pb[3],pc[3],pt[13];
    
    numeral(ia,pa,3);
    numeral(ib,pb,3);
    numeral(ic,pc,3);
    strcpy(pt,pa);
    strcat(pt,pb);
    strcat(pt,pc);
    for (i=0; i < angk1.nang; i++)
    {
        if (strcmp(pt,angk1.ktype[i]) == 0)
        {
            angk1.con[i] = f1;
            angk1.ang[i][0] = f2;
            angk1.ang[i][1] = f3;
            angk1.ang[i][2] = f4;
            return TRUE;
        }
    }
    return FALSE;
}         
/* ========================================== */    
int check_angle5const(int ia,int ib,int ic,float f1, float f2,float f3,float f4)
{
    int i;
    char pa[3],pb[3],pc[3],pt[13];
    
    numeral(ia,pa,3);
    numeral(ib,pb,3);
    numeral(ic,pc,3);
    strcpy(pt,pa);
    strcat(pt,pb);
    strcat(pt,pc);
    for (i=0; i < angk1.nang5; i++)
    {
        if (strcmp(pt,angk1.ktype5[i]) == 0)
        {
            angk1.con5[i] = f1;
            angk1.ang5[i][0] = f2;
            angk1.ang5[i][1] = f3;
            angk1.ang5[i][2] = f4;
            return TRUE;
        }
    }
    return FALSE;
}         
/* ========================================== */    
int check_angle4const(int ia,int ib,int ic,float f1, float f2,float f3,float f4)
{
    int i;
    char pa[3],pb[3],pc[3],pt[13];
    
    numeral(ia,pa,3);
    numeral(ib,pb,3);
    numeral(ic,pc,3);
    strcpy(pt,pa);
    strcat(pt,pb);
    strcat(pt,pc);
    for (i=0; i < angk1.nang4; i++)
    {
        if (strcmp(pt,angk1.ktype4[i]) == 0)
        {
            angk1.con4[i] = f1;
            angk1.ang4[i][0] = f2;
            angk1.ang4[i][1] = f3;
            angk1.ang4[i][2] = f4;
            return TRUE;
        }
    }
    return FALSE;
}         
/* ========================================== */    
int check_angle3const(int ia,int ib,int ic,float f1, float f2,float f3,float f4)
{
    int i;
    char pa[3],pb[3],pc[3],pt[13];
    
    numeral(ia,pa,3);
    numeral(ib,pb,3);
    numeral(ic,pc,3);
    strcpy(pt,pa);
    strcat(pt,pb);
    strcat(pt,pc);
    for (i=0; i < angk1.nang3; i++)
    {
        if (strcmp(pt,angk1.ktype3[i]) == 0)
        {
            angk1.con3[i] = f1;
            angk1.ang3[i][0] = f2;
            angk1.ang3[i][1] = f3;
            angk1.ang3[i][2] = f4;
            return TRUE;
        }
    }
    return FALSE;
}         
/* ========================================== */    
int check_bondconst(int ia,int ib,float f1, float f2)
{
    int i;
    char pa[3],pb[3],pt[13];
    
    numeral(ia,pa,3);
    numeral(ib,pb,3);
    strcpy(pt,pa);
    strcat(pt,pb);
    for (i=0; i < bondk1.nbnd; i++)
    {
        if (strcmp(pt,bondk1.kb[i]) == 0)
        {
            bondk1.s[i] = f1;
            bondk1.t[i] = f2;
            return TRUE;
        }
    }
    return FALSE;
}         
/* ========================================== */    
int check_bond5const(int ia,int ib,float f1, float f2)
{
    int i;
    char pa[3],pb[3],pt[13];
    
    numeral(ia,pa,3);
    numeral(ib,pb,3);
    strcpy(pt,pa);
    strcat(pt,pb);
    for (i=0; i < bondk1.nbnd5; i++)
    {
        if (strcmp(pt,bondk1.kb5[i]) == 0)
        {
            bondk1.s5[i] = f1;
            bondk1.t5[i] = f2;
            return TRUE;
        }
    }
    return FALSE;
}         
/* ========================================== */    
int check_bond4const(int ia,int ib,float f1, float f2)
{
    int i;
    char pa[3],pb[3],pt[13];
    
    numeral(ia,pa,3);
    numeral(ib,pb,3);
    strcpy(pt,pa);
    strcat(pt,pb);
    for (i=0; i < bondk1.nbnd4; i++)
    {
        if (strcmp(pt,bondk1.kb4[i]) == 0)
        {
            bondk1.s4[i] = f1;
            bondk1.t4[i] = f2;
            return TRUE;
        }
    }
    return FALSE;
}         
/* ========================================== */    
int check_bond3const(int ia,int ib,float f1, float f2)
{
    int i;
    char pa[3],pb[3],pt[13];
    
    numeral(ia,pa,3);
    numeral(ib,pb,3);
    strcpy(pt,pa);
    strcat(pt,pb);
    for (i=0; i < bondk1.nbnd3; i++)
    {
        if (strcmp(pt,bondk1.kb3[i]) == 0)
        {
            bondk1.s3[i] = f1;
            bondk1.t3[i] = f2;
            return TRUE;
        }
    }
    return FALSE;
}         

#include "mmxconst.h"
#include "mmff94.h"

/* =============================================  */
void read_datafiles(char *pfile)
{
  if(pfile) {
    if( pfile[0]!='|') { /* vertical bar signals use of in-memory parameter files */

      /* TO IMPLEMENT: read file into array a string for passing intoread_parameter */
#if 0
      char filename[256];
      
      char *ev = NULL;
      
      if ( (ev=getenv("ENGINE_DIR")) != NULL)
        {
          strcpy(filename,ev);
          strcat(filename,"/");
          strcat(filename,datafile);
        } else
        strcpy(filename,datafile);
      read_parameterfile(filename);
#endif
    } else {
      if (!strncmp(pfile, "|mmxconst",9)) {
        read_parameterfile(mmxconst);
      } else if (!strncmp(pfile, "|mmff94",7))  {
        read_parameterfile(mmff94);
      }
    }
    return;
  }
}
// =================================================
void read_parameterfile(char **parray)
{
        int iz;
        int ihtype,lpde,ihdonor;
        int i, ii, kk, kp1, kp2, kp3, kp4, kt1, kt2, kt3;
	int irec=0;
        int   se1[6];
        int iser, inumber, iligand, iclass, iclass1, iclass2, ivalence;
        float radius,epsilon,alpha,ntmp,atmp,gtmp,ftemp;
        float v1[6];
        float wtt;
        float charge, fchrg,tchrg;
        char symbol[3], descript[25];
        char dumm[21], iptemp[21],name[21];
        char pa[4],pb[4],pt[7];
        char pc[4],pd[4],pang[10],ptor[13];
	char *line;

/*
        FILE *datafile;
        datafile = fopen(string, "rt");

        if (datafile == NULL)
        {
           char message[80];
           sprintf(message,"Unable to open the data file %s.\nPlease check that this file was installed\n and is in the same directory as PCWIN",string);
           message_alert(message,"Read Parameter");
           fprintf(pcmlogfile,"Unable to open the data file %s.\nPlease check that this file was installed\n and is in the same directory as PCWIN\n",string);
           fclose(datafile);
           return;
         }

         while ( FetchRecord(datafile,line)) 
*/
	while (strlen(parray[irec]))
         {
	    line = parray[irec];
	    ++irec;
            sscanf(line,"%s",iptemp);
        /*  force field descriptors   */
            if (strcmp(iptemp,"forcefield") == 0)
            {
                sscanf(line,"%s %s",dumm, name);
		set_field_name(name);
                if (strcmp(name,"MMX") == 0)
		  set_field(MMX);
                else if (strcmp(name,"GAFF") == 0 )
		  set_field(GAFF);
                else if (strcmp(name,"MMFF94") == 0)
		  set_field(MMFF94);
                else
		  set_field(MMX);
            } else if (strcmp(iptemp,"bondunit") == 0)
            {
                sscanf(line,"%s %f",dumm, &ftemp);
		set_field_bondunit(ftemp);
                units.bndunit = ftemp;
            } else if (strcmp(iptemp,"bond-cubic") == 0)
            {
                sscanf(line,"%s %f",dumm, &ftemp);
		set_field_bondcubic(ftemp);
                units.cbnd = ftemp;
            } else if (strcmp(iptemp,"bond-quartic") == 0)
            {
                sscanf(line,"%s %f",dumm, &ftemp);
		set_field_bondquartic(ftemp);
                units.qbnd = ftemp;
            } else if (strcmp(iptemp,"angleunit") == 0)
            {
                sscanf(line,"%s %f",dumm, &ftemp);
		set_field_angleunit(ftemp);
                units.angunit = ftemp;
            } else if (strcmp(iptemp,"angle-cubic") == 0)
            {
                sscanf(line,"%s %f",dumm, &ftemp);
		set_field_anglecubic(ftemp);
                  units.cang = ftemp;
            } else if (strcmp(iptemp,"angle-quartic") == 0)
            {
                sscanf(line,"%s %f",dumm, &ftemp);
		set_field_anglequartic(ftemp);
                  units.qang = ftemp;
            } else if (strcmp(iptemp,"angle-pentic") == 0)
            {
                sscanf(line,"%s %f",dumm, &ftemp);
		set_field_anglepentic(ftemp);
                units.pang = ftemp;
            } else if (strcmp(iptemp,"angle-sextic") == 0)
            {
                sscanf(line,"%s %f",dumm, &ftemp);
		set_field_anglesextic(ftemp);
                units.sang = ftemp;
            } else if (strcmp(iptemp,"str-bndunit") == 0)
            {
                sscanf(line,"%s %f",dumm, &ftemp);
		set_field_strbndunit(ftemp);
                units.stbnunit = ftemp;
            } else if (strcmp(iptemp,"ang-angunit") == 0)
            {
                sscanf(line,"%s %f",dumm, &ftemp);
		set_field_angangunit(ftemp);
                units.aaunit = ftemp;
            } else if (strcmp(iptemp,"torsionunit") == 0)
            {
                sscanf(line,"%s %f",dumm, &ftemp);
		set_field_torsionunit(ftemp);
                units.torsunit = ftemp;
            } else if (strcmp(iptemp,"str-torunit") == 0)
            {
                sscanf(line,"%s %f",dumm, &ftemp);
		set_field_strtorunit(ftemp);
                units.storunit = ftemp;
            } else if (strcmp(iptemp,"vdwtype") == 0)
            {
                  sscanf(line,"%s %s",dumm, name);
		  set_field_vdwtype(name);
            } else if (strcmp(iptemp,"radiustype") == 0)
            {
                  sscanf(line,"%s %s",dumm, name);
		  set_field_radiustype(name);
            } else if (strcmp(iptemp,"radiussize") == 0)
            {
                  sscanf(line,"%s %s",dumm, name);
		  set_field_radiussize(name);
            } else if (strcmp(iptemp,"radiusrule") == 0)
            {
                  sscanf(line,"%s %s",dumm, name);
		  set_field_radiusrule(name);
            } else if (strcmp(iptemp,"epsilonrule") == 0)
            {
                  sscanf(line,"%s %s",dumm, name);
		  set_field_epsrule(name);
            } else if (strcmp(iptemp,"a-expterm") == 0)
            {
                  sscanf(line,"%s %f",dumm, &ftemp);
		  set_field_aterm(ftemp);
                  units.aterm = ftemp;
            } else if (strcmp(iptemp,"b-expterm") == 0)
            {
                  sscanf(line,"%s %f",dumm, &ftemp);
		  set_field_bterm(ftemp);
                  units.bterm = ftemp;
            } else if (strcmp(iptemp,"c-expterm") == 0)
            {
                  sscanf(line,"%s %f",dumm, &ftemp);
		  set_field_cterm(ftemp);
                  units.cterm = ftemp;
            } else if (strcmp(iptemp,"vdw-14-scale") == 0)
            {
                  sscanf(line,"%s %f",dumm, &ftemp);
		  set_field_vdwscale(ftemp);
            } else if (strcmp(iptemp,"dielectric") == 0)
            {
                  sscanf(line,"%s %f",dumm, &ftemp);
		  set_field_dielectric(ftemp);
                 if (user.dielec == FALSE)
                     units.dielec = ftemp;
            } else if (strcmp(iptemp,"chg-14-scale") == 0)
            {
                  sscanf(line,"%s %f",dumm, &ftemp);
		  set_field_chrgscale(ftemp);
                  units.chgscale = ftemp;
            } else if (strcmp(iptemp,"atom") == 0)    /* atom type descriptors */
            {
                if (atom_k.natomtype >= MAXATOMTYPE )
                {
                     message_alert("Error - Too many atom types","Read Parameter");
                     fprintf(pcmlogfile,"Maximum number of atom types exceeded. Others will be ignored\n");
                     //fclose(datafile);
                     return;
                }
                
                sscanf(line, "%s %d %s %22c %d %f %d %d %d %d %d",dumm,&iser, symbol, descript, &inumber, &wtt,
                              &iligand, &ivalence, &iclass, &iclass1, &iclass2);

                 atom_k.type[iser] = iser;
                 strcpy(atom_k.symbol[iser], symbol);
                 strcpy(atom_k.description[iser], descript);
                 atom_k.description[iser][20] = '\0';
                 atom_k.number[iser]  = inumber;
                 atom_k.weight[iser]  = wtt;
                 atom_k.ligands[iser] = iligand;
                 atom_k.tclass[iser]   = iclass;
                 atom_k.valency[iser] = ivalence;
                 atom_k.tclass1[iser]  = iclass1;
                 atom_k.tclass2[iser]  = iclass2;
                 atom_k.natomtype++;
            } else if ( strcmp(iptemp,"metal") == 0) /* metal data */
            {
                if (metaldata.nmetal > 200)
                {
                    message_alert("Error-Too many metals in Parameter file","Read Parameter");
                    //fclose(datafile);
                    return;
                }
                sscanf(line,"%s %d %s %f %f",dumm,&metaldata.type[metaldata.nmetal], metaldata.name[metaldata.nmetal],&metaldata.radius[metaldata.nmetal],
                    &metaldata.eps[metaldata.nmetal]);
                metaldata.nmetal++;
            } else if( strcmp(iptemp,"bond") == 0 )    /* bond constants  */
            {
                if (bondk1.nbnd >= MAXBONDCONST)
                {
                     message_alert("Error - Too many bond constants","Read Parameter");
                     fprintf(pcmlogfile,"Maximum number of bond constants exceeded. Others will be ignored\n");
                     //fclose(datafile);
                     return;
                }
                        
                sscanf( line, "%s %d %d %f %f", dumm, &ii, &kk, 
                         &bondk1.s[bondk1.nbnd], &bondk1.t[bondk1.nbnd]);
                numeral(ii,pa,3);
                numeral(kk,pb,3);
                strcpy(pt,pa);
                strcat(pt,pb);
                strcpy(bondk1.kb[bondk1.nbnd],pt);
                bondk1.nbnd++;
            } else if( strcmp(iptemp,"bond3") == 0 )    /* bond constants  */
            {
                if (bondk1.nbnd3 >= MAXBOND3CONST)
                {
                     message_alert("Error - Too many bond3 constants","Read Parameter");
                     fprintf(pcmlogfile,"Maximum number of bond3 constants exceeded. Others will be ignored\n");
                     //fclose(datafile);
                     return;
                }
                        
                sscanf( line, "%s %d %d %f %f", dumm, &ii, &kk, 
                         &bondk1.s3[bondk1.nbnd3], &bondk1.t3[bondk1.nbnd3]);
                numeral(ii,pa,3);
                numeral(kk,pb,3);
                strcpy(pt,pa);
                strcat(pt,pb);
                strcpy(bondk1.kb3[bondk1.nbnd3],pt);
//                bondk1.kb3[bondk1.nbnd3] = ii*100 + kk;
                bondk1.nbnd3++;
            } else if( strcmp(iptemp,"bond4") == 0 )    /* bond constants  */
            {
                if (bondk1.nbnd4 >= MAXBOND4CONST)
                {
                     message_alert("Error - Too many bond4 constants","Read Parameter");
                     fprintf(pcmlogfile,"Maximum number of bond4 constants exceeded. Others will be ignored\n");
                     //fclose(datafile);
                     return;
                }
                        
                sscanf( line, "%s %d %d %f %f", dumm, &ii, &kk, 
                         &bondk1.s4[bondk1.nbnd4], &bondk1.t4[bondk1.nbnd4]);
                numeral(ii,pa,3);
                numeral(kk,pb,3);
                strcpy(pt,pa);
                strcat(pt,pb);
                strcpy(bondk1.kb4[bondk1.nbnd4],pt);
//                bondk1.kb4[bondk1.nbnd4] = ii*100 + kk;
                bondk1.nbnd4++;
            } else if( strcmp(iptemp,"bond5") == 0 )    /* bond constants  */
            {
                if (bondk1.nbnd5>= MAXBOND5CONST)
                {
                     message_alert("Error - Too many bond5 constants","Read Parameter");
                     fprintf(pcmlogfile,"Maximum number of bond5 constants exceeded. Others will be ignored\n");
                     //fclose(datafile);
                     return;
                }
                        
                sscanf( line, "%s %d %d %f %f", dumm, &ii, &kk, 
                         &bondk1.s5[bondk1.nbnd5], &bondk1.t5[bondk1.nbnd5]);
                numeral(ii,pa,3);
                numeral(kk,pb,3);
                strcpy(pt,pa);
                strcat(pt,pb);
                strcpy(bondk1.kb5[bondk1.nbnd5],pt);
                bondk1.nbnd5++;
            } else if( strcmp(iptemp,"bonddel") == 0 )    /* bond constants  */
            {
                if (bondk1.ndeloc>= MAXBONDDELOC)
                {
                     message_alert("Error - Too many delocalized bond constants","Read Parameter");
                     fprintf(pcmlogfile,"Maximum number of delocalized bond constants exceeded. Others will be ignored\n");
                     //fclose(datafile);
                     return;
                }
                        
                sscanf( line, "%s %d %d %f %f", dumm, &ii, &kk, 
                         &bondk1.sdel[bondk1.ndeloc], &bondk1.tdel[bondk1.ndeloc]);
                numeral(ii,pa,3);
                numeral(kk,pb,3);
                strcpy(pt,pa);
                strcat(pt,pb);
                strcpy(bondk1.kbdel[bondk1.ndeloc],pt);
                bondk1.ndeloc++;
            } else if (strcmp(iptemp, "electroi") == 0)
            {
                sscanf(line,"%s %d %d %d %f",dumm,&electroneg.itype[electroneg.nelecti],
                       &electroneg.ibond[electroneg.nelecti],
                       &electroneg.iattach[electroneg.nelecti],
                       &electroneg.icorr[electroneg.nelecti]);
                electroneg.nelecti++;
            } else if (strcmp(iptemp, "electroj") == 0)
            {
                sscanf(line,"%s %d %d %d %f",dumm,&electroneg.jtype[electroneg.nelectj],
                       &electroneg.jbond[electroneg.nelectj],
                       &electroneg.jattach[electroneg.nelectj],
                       &electroneg.jcorr[electroneg.nelectj]);
                electroneg.nelectj++;
            } else if (strcmp(iptemp,"dipole") == 0)
            {
                if (dipole_k.ndipole > MAXBONDCONST)
                {
                    message_alert("Error - Too many bond dipoles","Read Parameter");
                    fprintf(pcmlogfile,"Maximum number of bond dipole constants exceeded. Others will be ignored\n");
                    //fclose(datafile);
                    return;
                }
                sscanf( line,"%s %d %d %f",dumm, &ii, &kk, &radius);
                numeral(ii,pa,3);
                numeral(kk,pb,3);
                strcpy(pt,pa); strcat(pt,pb);
                strcpy(dipole_k.kb[dipole_k.ndipole],pt);
                dipole_k.bmom[dipole_k.ndipole] = radius;
                dipole_k.ndipole++;
            } else if (strcmp(iptemp, "charge") == 0)
            {
                if (charge_k.ncharge > MAXATOMTYPE)
                {
                    message_alert("Error - Too many atom charges","Read Parameter");
                    fprintf(pcmlogfile,"Maximum number of atom charge constants exceeded. Others will be ignored\n");
                    //fclose(datafile);
                    return;
                }
                sscanf( line,"%s %d %f", dumm, &charge_k.type[charge_k.ncharge], &charge_k.charge[charge_k.ncharge]);
                charge_k.ncharge++;
            } else if (strcmp(iptemp, "mmffchrg") == 0)
            {
                if (charge_k.ncharge > MAXATOMTYPE)
                {
                    message_alert("Error - Too many atom charges","Read Parameter");
                   fprintf(pcmlogfile,"Maximum number of atom charge constants exceeded. Others will be ignored\n");
                    //fclose(datafile);
                    return;
                }
                sscanf( line,"%s %d %f %f %f", dumm, &i, &charge, &fchrg, &tchrg);
                 charge_k.type[i] = i;
                 charge_k.charge[i] = charge;
                 charge_k.formchrg[i] = fchrg;
                 charge_k.typechrg[i] = tchrg;
            } else if (strcmp(iptemp,"bndchrgdel") == 0)
            {
                if (charge_k.nbndchrgdel > MAXBONDCONST)
                {
                    message_alert("Error - Too many bond charges","Read Parameter");
                   fprintf(pcmlogfile,"Maximum number of bond charge constants exceeded. Others will be ignored\n");
                    //fclose(datafile);
                    return;
                }
                sscanf(line,"%s %d %d %f",dumm,&kp1,&kp2,&charge_k.bchargedel[charge_k.nbndchrgdel]);
                charge_k.btypedel[charge_k.nbndchrgdel] = kp1*100+kp2;
                charge_k.nbndchrgdel++;
            } else if (strcmp(iptemp,"bndchrg") == 0)
            {
                if (charge_k.nbndchrg > MAXBONDCONST)
                {
                    message_alert("Error - Too many bond charges","Read Parameter");
                   fprintf(pcmlogfile,"Maximum number of bond charge constants exceeded. Others will be ignored\n");
                    //fclose(datafile);
                    return;
                }
                sscanf(line,"%s %d %d %f",dumm,&kp1,&kp2,&charge_k.bcharge[charge_k.nbndchrg]);
                charge_k.btype[charge_k.nbndchrg] = kp1*100+kp2;
                charge_k.nbndchrg++;
            }  else if (strcmp(iptemp,"vdwmmff") == 0)
            {
                if (vdw1.nvdw > MAXVDWCONST)
                {
                    message_alert("Error - Too many vdw constants","Read Parameter");
                   fprintf(pcmlogfile,"Maximum number of vdw constants exceeded. Others will be ignored\n");
                    //fclose(datafile);
                    return;
                }
                sscanf( line, "%s %d %f %f %f %f %s", dumm, &i, &alpha, &ntmp, &atmp, &gtmp, symbol);
                vdw1.alpha[i] = alpha;
                vdw1.n[i] = ntmp;
                vdw1.a[i] = atmp;
                vdw1.g[i] = gtmp;
                strcpy(vdw1.da[i],symbol);
                vdw1.nvdw++;               
            }  else if( strcmp(iptemp,"vdw")  == 0 )
            {
                if (vdw1.nvdw > MAXVDWCONST)
                {
                    message_alert("Error - Too many vdw constants","Read Parameter");
                   fprintf(pcmlogfile,"Maximum number of vdw constants exceeded. Others will be ignored\n");
                    //fclose(datafile);
                    return;
                }
                sscanf( line, "%s %d %f %f %d %d %d", dumm, &i, &radius, &epsilon, &lpde, &ihtype, &ihdonor);

                vdw1.rad[i] = radius;
                vdw1.eps[i] = epsilon;
                vdw1.lpd[i] = lpde;
                vdw1.ihdon[i] = ihdonor;
                vdw1.ihtyp[i] = ihtype;
                vdw1.nvdw++;
             }
        /*   torsion read */
             else if( strcmp(iptemp,"torsion4") == 0 )
             {
                 if (torkn1.ntor4 >= MAXTOR4CONST)
                 {
                      message_alert("Error - Too many torsion 4 constants","Read Parameter");
                     fprintf(pcmlogfile,"Maximum number of torsion 4 constants exceeded. Others will be ignored\n");
                      //fclose(datafile);
                      return;
                 }
                for(i=0; i < 6; i++)
                   v1[i] = se1[i] = 0.0;
                 iz = sscanf(line, "%s  %d %d %d %d %f %d %f %d %f %d %f %d %f %d %f %d", dumm, &kp1,
                         &kp2, &kp3, &kp4 ,
                         &v1[0], &se1[0],&v1[1], &se1[1],&v1[2], &se1[2],
                         &v1[3], &se1[3],&v1[4], &se1[4],&v1[5], &se1[5]);
                 torphase(6, v1, se1);
                 torkn1.tv41[torkn1.ntor4] = v1[0];
                 torkn1.tv42[torkn1.ntor4] = v1[1];
                 torkn1.tv43[torkn1.ntor4] = v1[2];
                 torkn1.phase41[torkn1.ntor4] = se1[0];
                 torkn1.phase42[torkn1.ntor4] = se1[1];
                 torkn1.phase43[torkn1.ntor4] = se1[2];

                 numeral(kp1,pa,3);
                 numeral(kp2,pb,3);
                 numeral(kp3,pc,3);
                 numeral(kp4,pd,3);
                 strcpy(ptor,pa);
                 strcat(ptor,pb); strcat(ptor,pc); strcat(ptor,pd);
                 strcpy(torkn1.kv4[torkn1.ntor4],ptor);
                 torkn1.ntor4++;
             }else if( strcmp(iptemp,"torsion5") == 0 )
             {
                 if (torkn1.ntor5 >= MAXTOR5CONST)
                 {
                      message_alert("Error - Too many torsion 5 constants","Read Parameter");
                     fprintf(pcmlogfile,"Maximum number of torsion 5 constants exceeded. Others will be ignored\n");
                      //fclose(datafile);
                      return;
                 }
                for(i=0; i < 6; i++)
                   v1[i] = se1[i] = 0.0;
                 iz = sscanf(line, "%s  %d %d %d %d %f %d %f %d %f %d %f %d %f %d %f %d", dumm, &kp1,
                         &kp2, &kp3, &kp4 ,
                         &v1[0], &se1[0],&v1[1], &se1[1],&v1[2], &se1[2],
                         &v1[3], &se1[3],&v1[4], &se1[4],&v1[5], &se1[5]);
                 torphase(6, v1, se1);
                 torkn1.tv51[torkn1.ntor5] = v1[0];
                 torkn1.tv52[torkn1.ntor5] = v1[1];
                 torkn1.tv53[torkn1.ntor5] = v1[2];
                 torkn1.phase51[torkn1.ntor5] = se1[0];
                 torkn1.phase52[torkn1.ntor5] = se1[1];
                 torkn1.phase53[torkn1.ntor5] = se1[2];
                 numeral(kp1,pa,3);
                 numeral(kp2,pb,3);
                 numeral(kp3,pc,3);
                 numeral(kp4,pd,3);
                 strcpy(ptor,pa);
                 strcat(ptor,pb); strcat(ptor,pc); strcat(ptor,pd);
                 strcpy(torkn1.kv5[torkn1.ntor5],ptor);
                 torkn1.ntor5++;
             }else if( strcmp(iptemp,"torsiondel") == 0 )
             {
                 if (torkn1.ntordel >= MAXTORDEL)
                 {
                      message_alert("Error - Too many delocalized torsion  constants","Read Parameter");
                     fprintf(pcmlogfile,"Maximum number of delocalized torsion  constants exceeded. Others will be ignored\n");
                      //fclose(datafile);
                      return;
                 }
                for(i=0; i < 6; i++)
                   v1[i] = se1[i] = 0.0;
                 iz = sscanf(line, "%s  %d %d %d %d %d %f %d %f %d %f %d %f %d %f %d %f %d", dumm, &torkn1.torindex[torkn1.ntordel],
                         &kp1, &kp2, &kp3, &kp4 ,
                         &v1[0], &se1[0],&v1[1], &se1[1],&v1[2], &se1[2],
                         &v1[3], &se1[3],&v1[4], &se1[4],&v1[5], &se1[5]);
                 torphase(6, v1, se1);
                 torkn1.tvdel1[torkn1.ntordel] = v1[0];
                 torkn1.tvdel2[torkn1.ntordel] = v1[1];
                 torkn1.tvdel3[torkn1.ntordel] = v1[2];
                 torkn1.phasedel1[torkn1.ntordel] = se1[0];
                 torkn1.phasedel2[torkn1.ntordel] = se1[1];
                 torkn1.phasedel3[torkn1.ntordel] = se1[2];
                 numeral(kp1,pa,3);
                 numeral(kp2,pb,3);
                 numeral(kp3,pc,3);
                 numeral(kp4,pd,3);
                 strcpy(ptor,pa);
                 strcat(ptor,pb); strcat(ptor,pc); strcat(ptor,pd);
                 strcpy(torkn1.kvdel[torkn1.ntordel],ptor);
                 torkn1.ntordel++;
             }else if( strcmp(iptemp, "torsion") == 0 )
             {
                if (torkn1.ntor >= MAXTORCONST)
                {
                    message_alert("Error - Too many torsion constants","Read Parameter");
                   fprintf(pcmlogfile,"Maximum number of torsion constants exceeded. Others will be ignored\n");
                    //fclose(datafile);
                    return;
                }
                for(i=0; i < 6; i++)
                   v1[i] = se1[i] = 0.0;
                   
                sscanf( line, "%s %d %d %d %d %f %d %f %d %f %d %f %d %f %d %f %d", dumm,
                         &kp1, &kp2, &kp3, &kp4,
                         &v1[0], &se1[0],&v1[1], &se1[1],&v1[2], &se1[2],
                         &v1[3], &se1[3],&v1[4], &se1[4],&v1[5], &se1[5]);
                torphase(6, v1, se1);
                torkn1.tv1[torkn1.ntor] = v1[0];
                torkn1.phase1[torkn1.ntor] = se1[0];
                torkn1.tv2[torkn1.ntor] = v1[1];
                torkn1.phase2[torkn1.ntor] = se1[1];
                torkn1.tv3[torkn1.ntor] = v1[2];
                torkn1.phase3[torkn1.ntor] = se1[2];
                torkn1.tv4[torkn1.ntor] = v1[3];
                torkn1.phase4[torkn1.ntor] = se1[3];
                torkn1.tv5[torkn1.ntor] = v1[4];
                torkn1.phase5[torkn1.ntor] = se1[4];
                torkn1.tv6[torkn1.ntor] = v1[5];
                torkn1.phase6[torkn1.ntor] = se1[5];
                 numeral(kp1,pa,3);
                 numeral(kp2,pb,3);
                 numeral(kp3,pc,3);
                 numeral(kp4,pd,3);
                 strcpy(ptor,pa);
                 strcat(ptor,pb); strcat(ptor,pc); strcat(ptor,pd);
                 strcpy(torkn1.kv[torkn1.ntor],ptor);

//                torkn1.kv[torkn1.ntor] = ((kp1*100 + kp2)*100 + kp3)*100 + kp4;
                torkn1.ntor++;
             } else if( strcmp(iptemp, "imptors") == 0)
             {
                 if (improptor_k.nimptor >= MAXIMP)
                {
                    message_alert("Error - Too many improper torsion constants","Read Parameter");
                   fprintf(pcmlogfile,"Maximum number of improper torsion constants exceeded. Others will be ignored\n");
                    //fclose(datafile);
                    return;
                }
                sscanf( line, "%s %d %d %d %d %f %d %f %d %f %d", dumm,
                         &kp1, &kp2, &kp3, &kp4,
                         &v1[0], &se1[0], &v1[1], &se1[1], &v1[2], &se1[2]);
                torphase(3,v1,se1);
                improptor_k.v1[improptor_k.nimptor] = v1[0];
                improptor_k.v2[improptor_k.nimptor] = v1[1];
                improptor_k.v3[improptor_k.nimptor] = v1[2];
                improptor_k.ph1[improptor_k.nimptor] = se1[0];
                improptor_k.ph2[improptor_k.nimptor] = se1[1];
                improptor_k.ph3[improptor_k.nimptor] = se1[2];
                 numeral(kp1,pa,3);
                 numeral(kp2,pb,3);
                 numeral(kp3,pc,3);
                 numeral(kp4,pd,3);
                 strcpy(ptor,pa);
                 strcat(ptor,pb); strcat(ptor,pc); strcat(ptor,pd);
                 strcpy(improptor_k.kv[improptor_k.nimptor],ptor);
//                improptor_k.kv[improptor_k.nimptor] = ((kp1*100 + kp2)*100 + kp3)*100 + kp4;
                improptor_k.nimptor++;
                
             } else if( strcmp(iptemp, "improper") == 0)
             {
                 if (improptor_k.nimptor >= MAXIMP)
                {
                    message_alert("Error - Too many improper torsion constants","Read Parameter");
                   fprintf(pcmlogfile,"Maximum number of improper torsion constants exceeded. Others will be ignored\n");
                    //fclose(datafile);
                    return;
                }
                sscanf( line, "%s %d %d %d %d %f %f", dumm,
                         &kp1, &kp2, &kp3, &kp4,
                         &improptor_k.cimptor[improptor_k.nimptor], &improptor_k.tdi[improptor_k.nimptor]);
                 numeral(kp1,pa,3);
                 numeral(kp2,pb,3);
                 numeral(kp3,pc,3);
                 numeral(kp4,pd,3);
                 strcpy(ptor,pa);
                 strcat(ptor,pb); strcat(ptor,pc); strcat(ptor,pd);
                 strcpy(improptor_k.kv[improptor_k.nimptor],ptor);
                 improptor_k.nimptor++;
              } else if( strcmp(iptemp, "ureybrad") == 0)
              {
                  if (ureybrad_k.nurey_brad > MAXUREY)
                  {
                     message_alert("Error - Too many UreyBrad constants","Read Parameter");
                    fprintf(pcmlogfile,"Maximum number of ureybradley constants exceeded. Others will be ignored\n");
                     //fclose(datafile);
                     return;
                  }
                  sscanf(line, "%s %d %d %d %f %f",dumm,
                         &kp1, &kp2, &kp3, &ureybrad_k.ubconst[ureybrad_k.nurey_brad],
                         &ureybrad_k.ubdist[ureybrad_k.nurey_brad]);
                  numeral(kp1,pa,3);
                  numeral(kp2,pb,3);
                  numeral(kp3,pc,3);
                  strcpy(pang,pa); strcat(pang,pb); strcat(pang,pc);
                  strcpy(ureybrad_k.kubang[ureybrad_k.nurey_brad],pang);
                  ureybrad_k.nurey_brad++;
              } else if (strcmp(iptemp, "anglef") == 0)
              {
                  sscanf(line,"%s %d %d %d %f %f %f %f %f %f %f %f",dumm, &kt1,&kt2,&kt3, &angf.fcon[angf.nfang],
                    &angf.fc0[angf.nfang],&angf.fc1[angf.nfang],&angf.fc2[angf.nfang],&angf.fc3[angf.nfang],
                    &angf.fc4[angf.nfang],&angf.fc5[angf.nfang],&angf.fc6[angf.nfang]);
                  numeral(kt1,pa,3);
                  numeral(kt2,pb,3);
                  numeral(kt3,pc,3);
                  strcpy(pang,pa); strcat(pang,pb); strcat(pang,pc);
                  strcpy(angf.kftype[angf.nfang],pang);
                  angf.nfang++;
              } else if( strcmp(iptemp, "angle5") == 0 )
              {
                 if (angk1.nang5 >= MAXANG5CONST)
                 {
                     message_alert("Error - Too many angle 5 constants","Read Parameter");
                    fprintf(pcmlogfile,"Maximum number of angle 5 constants exceeded. Others will be ignored\n");
                     //fclose(datafile);
                     return;
                 }
                 sscanf( line, "%s %d %d %d %f %f %d", dumm, 
                         &kt1, &kt2, &kt3, &angk1.con5[angk1.nang5], &angk1.ang5[angk1.nang5][0],
			 &angk1.index5[angk1.nang5] );
                 numeral(kt1,pa,3);
                 numeral(kt2,pb,3);
                 numeral(kt3,pc,3);
                 strcpy(pang,pa); strcat(pang,pb); strcat(pang,pc);
                 strcpy(angk1.ktype5[angk1.nang5],pang);
                 angk1.nang5++;
             } else if( strcmp(iptemp, "angle4") == 0 )
             {
                 if (angk1.nang4 >= MAXANG4CONST)
                 {
                     message_alert("Error - Too many angle 4 constants","Read Parameter");
                    fprintf(pcmlogfile,"Maximum number of angle 4 constants exceeded. Others will be ignored\n");
                     //fclose(datafile);
                     return;
                 }
                 sscanf( line, "%s %d %d %d %f %f %d", dumm, 
                         &kt1, &kt2, &kt3, &angk1.con4[angk1.nang4], &angk1.ang4[angk1.nang4][0], 
                         &angk1.index4[angk1.nang4] );
                 numeral(kt1,pa,3);
                 numeral(kt2,pb,3);
                 numeral(kt3,pc,3);
                 strcpy(pang,pa); strcat(pang,pb); strcat(pang,pc);
                 strcpy(angk1.ktype4[angk1.nang4],pang);
                 angk1.nang4++;
             } else if( strcmp(iptemp, "angle3") == 0 )
             {
                if (angk1.nang3 >= MAXANG3CONST )
                {
                    message_alert("Error - Too many angle 3 constants","Read Parameter");
                   fprintf(pcmlogfile,"Maximum number of angle 3 constants exceeded. Others will be ignored\n");
                    //fclose(datafile);
                    return;
                }
                sscanf( line, "%s  %d %d %d %f %f %d", dumm, 
                        &kt1, &kt2, &kt3, &angk1.con3[angk1.nang3], &angk1.ang3[angk1.nang3][0], 
                         &angk1.index3[angk1.nang3] );
                 numeral(kt1,pa,3);
                 numeral(kt2,pb,3);
                 numeral(kt3,pc,3);
                 strcpy(pang,pa); strcat(pang,pb); strcat(pang,pc);
                 strcpy(angk1.ktype3[angk1.nang3],pang);
                angk1.nang3++;
           } else if( strcmp(iptemp, "angle") == 0 )
           {
               if (angk1.nang >= MAXANGCONST)
               {
                  message_alert("Error - Too many angle constants","Read Parameter");
                 fprintf(pcmlogfile,"Maximum number of angle constants exceeded. Others will be ignored\n");
                  //fclose(datafile);
                  return;
               }
               sscanf( line, "%s  %d %d %d %f %f %d", dumm,
                        &kt1, &kt2, &kt3, &angk1.con[angk1.nang], &angk1.ang[angk1.nang][0], 
                        &angk1.index[angk1.nang] );
                 numeral(kt1,pa,3);
                 numeral(kt2,pb,3);
                 numeral(kt3,pc,3);
                 strcpy(pang,pa); strcat(pang,pb); strcat(pang,pc);
                 strcpy(angk1.ktype[angk1.nang],pang);
               angk1.nang++;
           } else if( strcmp(iptemp, "angdel") == 0 )
           {
               if (angk1.ndel >= MAXANGDEL)
               {
                  message_alert("Error - Too many delocalized angle constants","Read Parameter");
                 fprintf(pcmlogfile,"Maximum number of delocalized angle constants exceeded. Others will be ignored\n");
                  //fclose(datafile);
                  return;
               }
               sscanf( line, "%s  %d %d %d %f %f %d", dumm,
                        &kt1, &kt2, &kt3, &angk1.condel[angk1.ndel], &angk1.angdel[angk1.ndel][0], 
		        &angk1.indexdel[angk1.ndel] );
                 numeral(kt1,pa,3);
                 numeral(kt2,pb,3);
                 numeral(kt3,pc,3);
                 strcpy(pang,pa); strcat(pang,pb); strcat(pang,pc);
                 strcpy(angk1.kdel[angk1.ndel],pang);
               angk1.ndel++;
           } else if( strcmp(iptemp, "ang3del") == 0 )
           {
               if (angk1.ndel3 >= MAXANG3DEL)
               {
                  message_alert("Error - Too many delocalized angle constants","Read Parameter");
                 fprintf(pcmlogfile,"Maximum number of delocalized angle constants exceeded. Others will be ignored\n");
                  //fclose(datafile);
                  return;
               }
               sscanf( line, "%s  %d %d %d %f %f %d", dumm,
                        &kt1, &kt2, &kt3, &angk1.condel3[angk1.ndel3], &angk1.angdel3[angk1.ndel3][0], 
                        &angk1.indexdel3[angk1.ndel3]);
                 numeral(kt1,pa,3);
                 numeral(kt2,pb,3);
                 numeral(kt3,pc,3);
                 strcpy(pang,pa); strcat(pang,pb); strcat(pang,pc);
                 strcpy(angk1.kdel3[angk1.ndel3],pang);
               angk1.ndel3++;
           } else if( strcmp(iptemp, "ang4del") == 0 )
           {
               if (angk1.ndel4 >= MAXANG4DEL)
               {
                  message_alert("Error - Too many delocalized angle constants","Read Parameter");
                 fprintf(pcmlogfile,"Maximum number of delocalized angle constants exceeded. Others will be ignored\n");
                  //fclose(datafile);
                  return;
               }
               sscanf( line, "%s  %d %d %d %f %f %d", dumm,
                        &kt1, &kt2, &kt3, &angk1.condel4[angk1.ndel4], &angk1.angdel4[angk1.ndel4][0], 
                        &angk1.indexdel4[angk1.ndel4] );
                 numeral(kt1,pa,3);
                 numeral(kt2,pb,3);
                 numeral(kt3,pc,3);
                 strcpy(pang,pa); strcat(pang,pb); strcat(pang,pc);
                 strcpy(angk1.kdel4[angk1.ndel4],pang);
               angk1.ndel4++;
          } else if( strcmp(iptemp, "opbend") == 0 )
          {
               if (ooplane_k.nopbend >= MAXOOP)
               {
                  message_alert("Error - Too many out of plane bending constants","Read Parameter");
                 fprintf(pcmlogfile,"Maximum number of out of plane bend constants exceeded. Others will be ignored\n");
                  //fclose(datafile);
                  return;
               }
               sscanf( line, "%s %d %d %d %d %f", dumm, &kp1, &kp2, &kp3,&kp4,
                         &ooplane_k.copb[ooplane_k.nopbend] );
               numeral(kp1,pa,3);
               numeral(kp2,pb,3);
               numeral(kp3,pc,3);
               numeral(kp4,pd,3);
               strcpy(ptor,pa); strcat(ptor,pb); strcat(ptor,pc); strcat(ptor,pd);
               strcpy(ooplane_k.iopb[ooplane_k.nopbend],ptor);
               ooplane_k.nopbend++;
          } else if (strcmp(iptemp, "strbnd") == 0)
          {
              if (crossterm_k.nstrbnd > MAXSTBN)
              {
                  message_alert("Error - Too many str bend constants","Read Parameter");
                 fprintf(pcmlogfile,"Maximum number of str bend constants exceeded. Others will be ignored\n");
                  //fclose(datafile);
                  return;
              }
              kp4 = 0;
              sscanf(line,"%s %d %d %d %f %f %f %d",dumm, &kp1, &kp2, &kp3,
                &crossterm_k.stbncon[crossterm_k.nstrbnd][0],
                &crossterm_k.stbncon[crossterm_k.nstrbnd][1],
                &crossterm_k.stbncon[crossterm_k.nstrbnd][2], &kp4);
               numeral(kp1,pa,3);
               numeral(kp2,pb,3);
               numeral(kp3,pc,3);
               strcpy(ptor,pa); strcat(ptor,pb); strcat(ptor,pc);
               strcpy(crossterm_k.stbn[crossterm_k.nstrbnd],ptor);
               crossterm_k.stbnindex[crossterm_k.nstrbnd] = kp4;
               crossterm_k.nstrbnd++;
          } else if (strcmp(iptemp, "angang") == 0)
          {
              if (crossterm_k.nangang > MAXAA)
              {
                  message_alert("Error - Too many ang ang constants","Read Parameter");
                 fprintf(pcmlogfile,"Maximum number of ang ang constants exceeded. Others will be ignored\n");
                  //fclose(datafile);
                  return;
              }
              sscanf(line,"%s %d %f %f %f",dumm, &crossterm_k.ang_ang[crossterm_k.nangang],
               &crossterm_k.aacon[crossterm_k.nangang][0],
               &crossterm_k.aacon[crossterm_k.nangang][1],
               &crossterm_k.aacon[crossterm_k.nangang][2]);
              crossterm_k.nangang++;
          } else if (strcmp(iptemp, "strtors") == 0)
          {
              if (crossterm_k.nstrtor > MAXSTRTOR)
              {
                  message_alert("Error - Too many str tor constants","Read Parameter");
                 fprintf(pcmlogfile,"Maximum number of str tor constants exceeded. Others will be ignored\n");
                  //fclose(datafile);
                  return;
              }
              sscanf(line,"%s %d %d %f",dumm, &kt1, &kt2, &crossterm_k.str_torcon[crossterm_k.nstrtor]);
              numeral(kt1,pa,3);
              numeral(kt2,pb,3);
              strcpy(pt,pa); strcat(pt,pb);
              strcpy(crossterm_k.str_tor[crossterm_k.nstrtor],pt);
              crossterm_k.nstrtor++;
          } else if (strcmp(iptemp, "vdwpr") == 0)
          {
              if (vdwpr_k .nvdwpr > MAXBONDCONST)
              {
                  message_alert("Error - Too many vdwpr constants","Read Parameter");
                 fprintf(pcmlogfile,"Maximum number of vdwpr constants exceeded. Others will be ignored\n");
                  //fclose(datafile);
                  return;
              }
              sscanf(line,"%s %d %d %f %f",dumm, &ii, &kk, &vdwpr_k.radius[vdwpr_k .nvdwpr],
              &vdwpr_k.eps[vdwpr_k .nvdwpr]);
              vdwpr_k.ia1[vdwpr_k .nvdwpr] = ii;
              vdwpr_k.ia2[vdwpr_k .nvdwpr] = kk;
              numeral(ii,pa,3);
              numeral(kk,pb,3);
              strcpy(pt,pa); strcat(pt,pb);
              strcpy(vdwpr_k.kv[vdwpr_k .nvdwpr],pt);
              vdwpr_k .nvdwpr++;
	  }
	 }
        //fclose(datafile);

/*
        fprintf(pcmlogfile," field : %s\n",field.name);
        fprintf(pcmlogfile," Atom Types: %d\n",atom_k.natomtype);
        fprintf(pcmlogfile," Bonds: %d Bond3: %d Bond4: %d Bond5: %d\n",
         bondk1.nbnd,bondk1.nbnd3, bondk1.nbnd4, bondk1.nbnd5);
        fprintf(pcmlogfile," Angle: %d Angle3: %d Angle4: %d Angle5: %d\n",
         angk1.nang,angk1.nang3, angk1.nang4, angk1.nang5);
        fprintf(pcmlogfile," Torsion: %d  Torsion4: %d Torsion5: %d\n",
         torkn1.ntor, torkn1.ntor4, torkn1.ntor5);
        fprintf(pcmlogfile," Vdw: %d OOP: %d Dipole: %d Charge: %d Improper: %d\n",
         vdw1.nvdw,ooplane_k.nopbend, dipole_k.ndipole, charge_k.ncharge, improptor_k.nimptor);
        fprintf(pcmlogfile," STBN: %d ANGANG: %d STRTOR: %d VDWPR: %d\n",
         crossterm_k.nstrbnd, crossterm_k.nangang, crossterm_k.nstrtor, vdwpr_k.nvdwpr);
*/
         
        return;
}


void torphase(int icount, float v1[6], int se[6])
{
    int i;
    float amp[6], phase[6];
    int fold[6];

    for (i=0; i < 6; i++)
    {
        fold[i] = 0;
        amp[i] = phase[i] = 0.0;
    }

    for (i=0; i < icount; i++)
    {
        amp[i] = v1[i];
        fold[i] = abs( (int)(2.0*se[i]) - (int) se[i]);
        if (se[i] > 0.0)
          phase[i] = 1.0;
        else
          phase[i] = -1.0;
        v1[i] = 0.0;
        se[i] = 0;
    }

    for (i=0; i < icount; i++)
    {
        if (fold[i] != 0 && fold[i] <= icount)
        {
            v1[fold[i]-1] = amp[i];
            se[fold[i]-1] = phase[i];
        }
    }
}

    
