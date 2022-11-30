
// need
// define MAXFIX
#define MAXFIX 1000
// fixed atom data
//    use inactive atoms to fix
//
EXTERN struct fx_atom {
    int natom_fix;
    int katom_fix[MAXFIX];
    } fx_atom;
// fixed distance
//    ndfix - number of fixed distances
//    kdfix = atom number of fixed dist
//    dfix = force constant and distance range
//
//  restrained atoms - use harmonic restraint to fix atoms
EXTERN struct restrain_atom {
    int natom_restrain;
    int katom_restrain[MAXFIX];
  float restrain_const[MAXFIX], restrain_max[MAXFIX],restrain_position[MAXFIX][3];
    } restrain_atom;
//
EXTERN struct fx_dist {
    int ndfix;
    int kdfix[MAXFIX][2];
    float fdconst[MAXFIX],min_dist[MAXFIX],max_dist[MAXFIX];
    } fx_dist;
//
// fixed angle
//    nafix = number of fixed angles
//    kafix = atom numbers
//    afix = force constant and angle range
EXTERN struct fx_angle {
    int nafix;
    int kafix[MAXFIX][3];
    float faconst[MAXFIX],min_ang[MAXFIX],max_ang[MAXFIX];
    } fx_angle;
//
// fixed torsion
//    ntfix = number of fixed torsions
//    ktfix = atom numbers
//    tfix = force constant and torsion range
EXTERN struct fx_torsion {
    int ntfix;
    int ktfix[MAXFIX][4];
    float ftconst[MAXFIX],min_tor[MAXFIX],max_tor[MAXFIX];
    } fx_torsion;

