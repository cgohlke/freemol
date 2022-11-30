#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stddef.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>
//#include <malloc.h>
//#include <direct.h>

#include "pcmsiz.h"

/* These flags allow us to choose whether to do the extra calculations */
#define DO_VIBRATION     2
#define DO_DIPOLE        4
#define DO_XLOGP         8
#define DO_ADDH         16 
#define DO_TEST         32
#define DO_DEBUG        64
//  flags definitions
#define PI_MASK                 0
#define AROMATIC_MASK           2
// type rules
#define NO_RETYPE               8
// Ring Size
#define RING3                   20
#define RING4                   21
#define RING5                   22
#define RING6                   23

// Force Field Names
#define         MMX             1
#define         MM2             2
#define         GAFF            3
#define         MMFF94          7
#define         UNKNOWN         10

/* Global to indicate verbose output or not */
EXTERN int VERBOSE;

#define strnicmp  strncasecmp

#define TRUE  1
#define FALSE 0
#define True  TRUE
#define False FALSE
#define true  TRUE
#define false FALSE

#ifndef PI   /* Avoid Linux Warnings! */
#define PI   3.14159265358979323846
#endif

#define radian 57.29577951308

#define Rad2Deg      (180.0/PI)
#define Deg2Rad      (PI/180.0)

#define AbsFun(a)    (((a)<0)? -(a) : (a))
#define MinFun(a,b)  (((a)<(b))? (a) : (b) )
#define MaxFun(a,b)  (((a)>(b))? (a) : (b) )

EXTERN struct t_units {
        double bndunit, cbnd, qbnd;
        double angunit, cang, qang, pang, sang, aaunit;
        double stbnunit, ureyunit, torsunit, storunit, v14scale;
        double aterm, bterm, cterm, dielec, chgscale;
        } units;


struct FileInfoStruct {
  int  ftype ;
  char   path[255] ;
  char   fname[255] ;
} ;

typedef struct FileInfoStruct Boxstruct;
EXTERN Boxstruct Openbox ,Savebox;

// global routines
void message_alert(char *, char *);
EXTERN FILE             *pcmlogfile;
EXTERN char              testfilename[50];
EXTERN int              MAXATOM, MAXBND,MAXANG,MAXTOR;

