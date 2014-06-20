# Make file for tnl-dynamics
# By Yuan-Chung Cheng <yccheng@berkeley.edu>
#

CC = gcc
CFLAGS = -O3 -finline-functions -Wall -Winline -funroll-loops -fexpensive-optimizations
#CFLAGS = -g -Wall -DDEBUG -DDEBUGDEBUG
#CFLAGS = -pg -fPIC -O2 -finline-functions -Wall -Winline -funroll-loops -fexpensive-optimizations

#LIBDIR = -L/usr/local/lib -L../lib
LIBDIR = -L/usr/local/lib -L/sw/lib -L../lib -L/home/ntucc/platin/opt/lib
EXTRALIBS = -lqdas -lpthread -lgsl -lgslcblas -lm 
#GSL_LAPACK = -lg2c -lblas -llapack
INCLUDEDIR = -I/usr/local/include -I/usr/X11R6/include -I../lib -I/sw/include
#INCLUDEDIR = -I/usr/local/include -I/usr/X11R6/include -I/sw/include -I../lib -I/home/ntucc/platin/opt/include

INSTALL = /usr/bin/install -c
BINDIR = /usr/local/chem/qdas

BINARYS = tnl-dynamics tnl-dm3pes tnl-lineshape
OBJECTS = tnl-kern.o

all: $(OBJECTS) $(BINARYS)

%.o : %.c %.h
	$(CC) -c $(CFLAGS) $(INCLUDEDIR) $< -o $@

tnl-dynamics: tnl-dynamics.c tnl-dynamics.h $(OBJECTS)
	$(CC) $(CFLAGS) $(INCLUDEDIR) $(LIBDIR) -o tnl-dynamics tnl-dynamics.c $(OBJECTS) $(EXTRALIBS)

tnl-dm3pes: tnl-dm3pes.c tnl-dm3pes.h $(OBJECTS)
	$(CC) $(CFLAGS) $(INCLUDEDIR) $(LIBDIR) -o tnl-dm3pes tnl-dm3pes.c $(OBJECTS) $(EXTRALIBS)

tnl-lineshape: tnl-lineshape.c $(OBJECTS)
	$(CC) $(CFLAGS) $(INCLUDEDIR) $(LIBDIR) -o tnl-lineshape tnl-lineshape.c $(OBJECTS) $(EXTRALIBS)

static: tnl-dynamics.c tnl-lineshape.c tnl-dynamics.h tnl-dm3pes.c tnl-dm3pes.h $(OBJECTS)
	$(CC) -static $(CFLAGS) $(INCLUDEDIR) $(LIBDIR) -o tnl-dynamics tnl-dynamics.c $(OBJECTS) $(EXTRALIBS)
	$(CC) -static $(CFLAGS) $(INCLUDEDIR) $(LIBDIR) -o tnl-dm3pes tnl-dm3pes.c $(OBJECTS) $(EXTRALIBS)
	$(CC) -static $(CFLAGS) $(INCLUDEDIR) $(LIBDIR) -o tnl-lineshape tnl-lineshape.c $(OBJECTS) $(EXTRALIBS)

install: static
	mkdir -p $(BINDIR)
	$(INSTALL) -m 755 $(BINARYS) $(BINDIR)

clean:
	rm -f $(BINARYS) *.o *.so *~ *.mo core a.out gmon.out 

#
# $Log$
# Revision 1.11  2007/07/15 04:51:18  platin
#   - add info on RNG.
#
# Revision 1.10  2007/06/18 19:28:16  platin
#
#   - add cvs revision tags.
#
# Revision 1.9  2007/06/01 18:00:02  platin
#
#   - bump to 0.8p1
#   - improved the Gaussian-Hermite module; now users can select between
#     1-20 point rules can be used.
#   - added keyword to select between MC and GH samplings.
#   - static disorder average loop in both tnl-dm3pes and tnl-dynamics
#     is fixed.
#
# Revision 1.8  2007/05/07 23:51:04  platin
#
#   - minor changes.
#
# Revision 1.7  2007/03/06 23:32:01  platin
#
#   - remove useless -markov codes.
#   - use SS initial condition for all bath aux functions.
#
# Revision 1.6  2007/02/20 23:12:34  platin
#
#   - add a program for TNL Markov dynamics.
#
# Revision 1.5  2007/02/06 06:24:48  platin
#
#   - add a prog. that can be used to test 2D ideas.
#
# Revision 1.4  2007/01/19 20:28:15  platin
#
#   - change def. of t in dm3pes according to normal convention.
#
# Revision 1.2  2006/10/30 06:32:35  platin
#
#   - change source of parser lib.
#
# Revision 1.1  2006/10/27 22:24:54  platin
#
#   - commit re-organized tnl-dm3pes and tnl-dynamics code.
#
# Revision 1.2  2006/08/23 22:30:26  platin
#
#   - add support for two-exciton states in codes regarding H and
#     transition dipoles. The "Assign" format for transition dipole
#     input can also be used to include effects of excited state absorption.
#
#   - basic support for TESLIST keyword.
#
# Revision 1.1  2006/07/11 16:46:32  platin
#
#   - import the tnl-dynamics bit; still requires modification to take
#     initial states etc.
#
#
# 
