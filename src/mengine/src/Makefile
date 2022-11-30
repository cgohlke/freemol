# NOTICE: this source code file has been modified for use with FreeMOL 
CC = gcc

# For FreeBSD uncomment the following
# CFLAGS = -g -Wall -O -I/usr/X11R6/include -DFREEBSD

# For RedHat Linux; if you want warnings trys "-Wall" and 
# for debugging "-g".  Most of us (users) like it like this: 
#CFLAGS =  -g -O3 -ftree-vectorize
#CFLAGS = -g -O3 -ftree-vectorize -fopenmp -DKEG_OPENMP
#CFLAGS = -g -O3 
CFLAGS = -g -O3 -ansi -std=gnu99 -pedantic -Wall -Wstrict-prototypes
#CFLAGS = -g -Wall -Wstrict-prototypes -pg -O -I/usr/local/LessTif/Motif1.2/include
#CFLAGS = -g -Wall -Wstrict-prototypes -pg -O -I/usr/local/Motif2.0/include

LDFLAGS = -lm
#LDFLAGS =  -lgomp -lm

RM = rm -f
CP = cp -f

# destination dir; set this for your setup
#BIN_DIR = /usr/local/bin

DST = mengine.exe

SRCS = \
active.c	angles.c	attach.c 	asnsym.c	datastat.c 	dipmom.c	draw.c \
eangle.c	ebond.c		ebufchrg.c	egeom.c		ehal.c		eobpw.c		esolv.c \
estrbnd.c	etorsion.c	field.c		get_mem.c	gmmx_run.c	initiali.c	kangle.c \
kbond.c		kcharge.c	kimptors.c	kopend.c	ksolv.c		kstrbnd.c	ktorsion.c \
kvdw.c		mengine.c	minimize.c	numeral.c	pcm7.c		pcmwin1.c	pcmwin2.c \
read.c		read_sdf.c	readprm.c	rings.c		search.c	solve.c \
tncg.c		torsions.c	type_mmx.c	utility.c	vibrate.c	xlogp.c \
elj.c		eimptors.c      echarge.c       read_syb.c	read_cml.c 

OBJS = \
active.o	angles.o	attach.o 	asnsym.o	datastat.o 	dipmom.o	draw.o \
eangle.o	ebond.o		ebufchrg.o	egeom.o		ehal.o		eobpw.o		esolv.o \
estrbnd.o	etorsion.o	field.o		get_mem.o	gmmx_run.o	initiali.o	kangle.o \
kbond.o		kcharge.o	kimptors.o	kopend.o	ksolv.o		kstrbnd.o	ktorsion.o \
kvdw.o		mengine.o	minimize.o	numeral.o	pcm7.o		pcmwin1.o	pcmwin2.o \
read.o		read_sdf.o	readprm.o	rings.o		search.o	solve.o \
tncg.o		torsions.o	type_mmx.o	utility.o	vibrate.o	xlogp.o \
elj.o 		eimptors.o      echarge.o       read_syb.o      read_cml.o	 


       
${DST}: $(OBJS)
	${CC} -g -o ${DST} ${OBJS}  ${LDFLAGS}
