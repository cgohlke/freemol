
CC = gcc

OBJS = mfwddct.o postdct.o huff.o bitio.o mheaders.o \
mpeg.o subsample.o param.o rgbtoycc.o \
	readframe.o combine.o jrevdct.o frame.o fsize.o frametype.o \
	libpnmrw.o specifics.o rate.o opts.o \
    iframe.o pframe.o bframe.o psearch.o bsearch.o block.o \
	noparallel.o nojpeg.o main.o

mpeg_encode: $(OBJS)
	gcc -o mpeg_encode.exe *.o -lm
