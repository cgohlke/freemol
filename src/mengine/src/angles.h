#define HARMONIC  1
#define FOURIER   2

EXTERN struct t_angles {
        int nang, **i13, *index;
        int nopb, *angin,*angtype;
        float *acon, *anat, *copb;
        } angles;

void get_angles(void);

