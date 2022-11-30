
/* EXTERN  struct t_hess {
    double hessx[MAXATOM][3], hessy[MAXATOM][3], hessz[MAXATOM][3];
    } hess;  */
    
 EXTERN struct t_hess {
    float **hessx, **hessy, **hessz;
    } hess; 

void hessian(int, double *,int *, int *, int *,double *);
   

