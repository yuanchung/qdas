# Make file for qas
# By Yuan-Chung Cheng <yccheng@berkeley.edu>
#

CC = gcc
CFLAGS = 

OBJECTS = 
BINARIES = 

INSTALL = /usr/bin/install -c
BINDIR = /usr/local/chem/qdas
INSTALL_BINARIES = qs_archive.sh qs_cleanall.sh qs_tnl-dm3pes_para.sh qs_tnl-dm3pes_para_nocache.sh

install: 
	mkdir -p $(BINDIR)
	$(INSTALL) -m 755 $(INSTALL_BINARIES) $(BINDIR)

clean:
	rm -f $(OBJECTS) $(BINARIES)
	rm -f *.a *.o *.so *~ *.mo core a.out gmon.out 

#
# $Log$
#
