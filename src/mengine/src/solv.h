#define MAXCELL 1500

EXTERN struct t_solv {
    int nmols;
    int nxcell, nycell, nzcell;
    int ncell, icell[MAXCELL][3];
    double density, boxx,boxy,boxz;
    double xcell,ycell,zcell,xcell2,ycell2,zcell2;
    } solv;

