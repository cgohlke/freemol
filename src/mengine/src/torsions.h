EXTERN struct t_torsions {
        int ntor, **i14;
        float *v1,*v2,*v3,*v4,*v5,*v6;
        int   *ph1,*ph2,*ph3,*ph4,*ph5,*ph6;
        } torsions;

/* struct t_improp {
        int nimprop, nimptors, iiprop[MAXTOR][4];
        float cimp[MAXTOR],tdi[MAXTOR];
        float v1[MAXTOR], v2[MAXTOR], v3[MAXTOR];
        int   ph1[MAXTOR], ph2[MAXTOR], ph3[MAXTOR];        
        } improp; */

/*  EXTERN struct t_torsions {
        int ntor, i14[MAXTOR][4];
        float v1[MAXTOR],v2[MAXTOR],v3[MAXTOR],v4[MAXTOR],v5[MAXTOR],v6[MAXTOR];
        int   ph1[MAXTOR],ph2[MAXTOR],ph3[MAXTOR],ph4[MAXTOR],ph5[MAXTOR],ph6[MAXTOR];
        } torsions;  */

