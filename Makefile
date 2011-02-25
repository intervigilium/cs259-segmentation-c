
# makefile
#
# Two Phase Active Contours without Edges (Chan-Vese)
#
#
# Igor Yanovsky
#
# Created at:    March 19, 2010
# Last modified: March 19, 2010


# BEGINNING OF MAKEFILE
TARGET  = twophase

OBJDIR  = $(HOMEDIR)
LIBDIR  = $(HOMEDIR)
CURRDIR = $(HOMEDIR)
OUTPUT  = $(CURRDIR)

SRCDIR  = $(CURRDIR):$(OBJDIR)

vpath %.c $(SRCDIR)

# put all the used .cpp file here:
ALLSRC = IORoutines.c  TwoPhase3D_driver.c  TwoPhase3DRoutines.c  timer.c  
ALLOBJ = IORoutines.o    TwoPhase3D_driver.o    TwoPhase3DRoutines.o    timer.o

#ALLOBJ = $(ALLSRC:%.c=$(OUTPUT)/%.o)

#FFTW3 = /home/yanovsky/fftw/current          # for MIPLDEVLINUX

PAPIDIR= /mnt/jc5/CS259/papi
INCLUDEDIR = -I$(CURRDIR) -I$(OBJDIR) -I$(PAPIDIR)
CFLAGS = -g -pg

LDFLAGS = -L/mnt/jc5/CS259/papi/ -lm -lutil_papi -lpapi

CC = gcc

%.o : %.c 
	$(CC) $(CFLAGS) $(INCLUDEDIR) -c $< -o $@

$(TARGET) : $(ALLOBJ)
	$(CC) $(CFLAGS) $(INCLUDEDIR) -o $@ $^ $(LDFLAGS)

depend : $(ALLSRC)
	makedepend $(CFLAGS) $(INCLUDEDIR) $^


.PHONY : clean
clean:
	rm $(TARGET) $(ALLOBJ)


# DO NOT DELETE
