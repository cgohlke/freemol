#define EXTERN extern

#include "pcwin.h"
#include "fix.h"

void set_active(int natom,int *use,int natom_fix,int *katom_fix);

void set_active(int natom,int *use,int natom_fix,int *katom_fix)
{
   int i;
   for (i=1; i <= natom; i++)
          use[i] = TRUE;

// fixed atoms
   for (i=0; i < natom_fix; i++)
       use[katom_fix[i]] = FALSE;
} 
