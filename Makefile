# Top-level makefile for the qdas package
# By Yuan-Chung Cheng <yuanchung@ntu.edu.tw>
#

CC = gcc
CFLAGS = 

all:
	make -C lib all
	make -C src all
	make -C tools all

clean:
	rm -f *.a *.o *.so *~ *.mo core a.out gmon.out
	make -C lib clean
	make -C src clean
	make -C tools clean

install:
	make -C src install
	make -C tools install

#
# $Log$
# Revision 1.3  2006/10/30 06:33:01  platin
#   - change source of parser lib.
#   - bump version to 0.5.
#
# Revision 1.2  2006/05/26 23:13:54  platin
#
#   - bug fix and minor changes in makefiles.
#
# Revision 1.1.1.1  2006/05/24 00:42:18  platin
#
#   - initial import of the qdas package.
#   - qdas stands for quantum dynamics and spectroscopy.
#   - basic code inherited from the "lineshape" package.
#
#
#
