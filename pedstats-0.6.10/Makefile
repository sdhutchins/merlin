#
# PEDSTATS Makefile -- Compiles and installs PEDSTATS
# (c) 2000-2007 Goncalo Abecasis, Janis Wigginton
#

# version information
PEDVERSION = 0.6.10

# default installation directory
INSTALLDIR = /usr/local/bin/

# default C++ compiler
CXX=g++

# default compilation flags are CFLAGS=-O3 -I./libsrc/
# 
# The following special options may also be added to the default
# 
#      Option                        Effect
#      -D_FILE_OFFSET_BITS=64        Enables support for swapfiles larger
#                                    than 2GB on supported systems (tested
#                                    on Linux)
#      -D__USE_LONG_INT              Enables support for markers with up
#                                    to 64 alleles (default is 32). Tested
#                                    on systems where gcc supports the long
#                                    long data type and on Windows.
# 
CFLAGS=-O2 -I./libsrc -I./pedstats -I./janpdf -D__ZLIB_AVAILABLE__

# executable file names and locations
BINDIR = executables
PEDSTATS = $(BINDIR)/pedstats
EXECUTABLES = $(PEDSTATS) 

# PEDSTATS File Set
PEDSTATSSRC = pedstats/Pedstats.cpp pedstats/PedstatsPDF.cpp \
 pedstats/PedstatsAgeCheck.cpp pedstats/PedstatsHWE.cpp pedstats/FilteredPedigree.cpp \
 pedstats/FilteredFamilyStats.cpp pedstats/PedstatsQuality.cpp pedstats/Manners.cpp \
 pedstats/PedigreePairs.cpp pedstats/SimplePairs.cpp
PEDSTATSHDR = $(PEDSTATSSRC:.cpp=.h)
PEDSTATSOBJ = $(PEDSTATSSRC:.cpp=.o)

# Utility Library File Set
LIBFILE = libsrc/lib-goncalo.a
LIBMAIN = libsrc/Error libsrc/IBD libsrc/IntArray libsrc/LongArray libsrc/StringHash libsrc/Hash \
 libsrc/MapFunction libsrc/MathGold libsrc/InputFile \
 libsrc/MathMatrix libsrc/MathStats libsrc/MathVector \
 libsrc/MemoryInfo libsrc/Parameters libsrc/Pedigree libsrc/PedigreeAlleleFreq \
 libsrc/PedigreeDescription libsrc/PedigreeFamily libsrc/PedigreeGlobals \
 libsrc/PedigreePerson libsrc/QuickIndex libsrc/Random libsrc/Sort \
 libsrc/StringArray libsrc/StringBasics libsrc/StringMap libsrc/FortranFormat libsrc/GenotypeLists \
 libsrc/Kinship
LIBPED = libsrc/PedigreeLoader libsrc/PedigreeTwin libsrc/PedigreeTrim
LIBSRC = $(LIBMAIN:=.cpp) $(LIBPED:=.cpp)
LIBHDR = $(LIBMAIN:=.h) libsrc/Constant.h \
 libsrc/MathConstant.h libsrc/PedigreeAlleles.h libsrc/LongInt.h
LIBOBJ = $(LIBSRC:.cpp=.o)
 
# PDF Library File Sets
PDFLIB = janpdf/libpdf.a
PDFFILES = janpdf/PDF janpdf/PDFfont janpdf/PDFinfo janpdf/PDFpage \
 janpdf/PDFpageobject janpdf/PDFchartbasics janpdf/PDFhistogram janpdf/PDFlinechart \
 janpdf/PDFchartaxis janpdf/PDFchartlegend janpdf/PDFchartmarker \
 janpdf/PDFchartline janpdf/PDFchartbar janpdf/PDFchartobject \
 janpdf/PDFgridcell janpdf/PDFgrid
PDFSRC = $(PDFFILES:=.cpp)
PDFHDR = $(PDFFILES:=.h)
PDFOBJ = $(PDFFILES:=.o)

# private parameters
FETCHDIR=$(HOME)/code
DISTRIBDIR=$(HOME)/code/distrib/pedstats-$(PEDVERSION)

# helpful screen listing available options
help : 
	@echo "PEDSTATS Source Distribution"
	@echo " "
	@echo "This Makefile will compile and install pedstats on your system"
	@echo " "
	@echo "Type...           To..."
	@echo "make help         Display this help screen"
	@echo "make all          Compile pedstats"
	@echo "make install      Install binaries in $(INSTALLDIR)"
	@echo "make install INSTALLDIR=directory_for_binaries"
	@echo "                  Install binaries in directory_for_binaries"
	@echo "make clean        Delete temporary files"

# make everything
all : $(EXECUTABLES)

$(EXECUTABLES) : $(BINDIR)

$(BINDIR) :
	mkdir $(BINDIR)

# dependencies for executables
$(PEDSTATS) : $(LIBFILE) $(PDFLIB) $(PEDSTATSOBJ) $(LIBHDR)
	$(CXX) $(CFLAGS) -o $@ $(PEDSTATSOBJ) ${PDFLIB} $(LIBFILE) -lm -lz

$(LIBFILE) : $(LIBOBJ)
	ar -cr $@ $(LIBOBJ)
	ranlib $@

$(PDFLIB) : $(PDFOBJ) 
	ar -cr $@ $(PDFOBJ)
	ranlib $@

$(PEDSTATSOBJ) : $(PEDSTATSHDR) $(LIBHDR) $(PDFHDR)

$(LIBOBJ) : $(LIBHDR)

$(PDFOBJ) : $(PDFHDR)

clean :
	-rm -f */*.a */*.o $(EXECUTABLES) 

install : all $(INSTALLDIR)
	@echo " "
	@echo Installing to directory $(INSTALLDIR)
	@echo To select a different directory, run
	@echo " "
	@echo make install INSTALLDIR=your_preferred_dir
	@echo " "
	@cp $(EXECUTABLES) $(INSTALLDIR)

$(INSTALLDIR) :
	@echo " "
	@echo Creating directory $(INSTALLDIR)
	@echo " "
	@mkdir $(INSTALLDIR)

new-version :
	mkdir -p $(DISTRIBDIR) $(DISTRIBDIR)/pedstats
	mkdir -p $(DISTRIBDIR)/libsrc $(DISTRIBDIR)/janpdf
	cp ChangeLog README $(DISTRIBDIR)
	cp Makefile $(DISTRIBDIR)
	cp -R examples $(DISTRIBDIR)
	
fetch : 
	cd $(FETCHDIR) ; cp $(PEDSTATSSRC) $(PEDSTATSHDR) $(DISTRIBDIR)/pedstats
	cd $(FETCHDIR) ; cp $(LIBSRC) $(LIBHDR) $(DISTRIBDIR)/libsrc
	cd $(FETCHDIR) ; cp $(PDFSRC) $(PDFHDR) $(DISTRIBDIR)/janpdf
	cd $(DISTRIBDIR); csh ../stamp PEDSTATS

.c.o :
	$(CXX) $(CFLAGS) -o $@ -c $*.c

.cpp.X.o : 
	$(CXX) $(CFLAGS) -o $@ -c $*.cpp -DPEDVERSION=\"$(PEDVERSION)\" -D__CHROMOSOME_X__

.cpp.o : 
	$(CXX) $(CFLAGS) -o $@ -c $*.cpp -DPEDVERSION=\"$(PEDVERSION)\"

archive : clean
	mkdir -p pedstats-$(PEDVERSION)
	cp -R README Makefile ChangeLog pedstats-$(PEDVERSION)
	cp -R pedstats libsrc janpdf examples pedstats-$(PEDVERSION)
	tar -cvf pedstats-$(PEDVERSION).tar pedstats-$(PEDVERSION) 
	gzip -f --best pedstats-$(PEDVERSION).tar
	rm -rf pedstats-$(PEDVERSION)

distrib : $(EXECUTABLES)
	mkdir -p pedstats-$(PEDVERSION)
	cp -R README ChangeLog $(EXECUTABLES) examples pedstats-$(PEDVERSION)
	tar -cvf `uname`-pedstats.tar pedstats-$(PEDVERSION)
	gzip -f `uname`-pedstats.tar
	rm -rf pedstats-$(PEDVERSION)

windowszip : $(EXECUTABLES)
	mkdir -p pedstats-$(PEDVERSION)
	cp -R README ChangeLog $(EXECUTABLES) examples pedstats-$(PEDVERSION)
	echo '@ECHO OFF' > pedstats-$(PEDVERSION)/STARTHERE.bat
	echo 'CMD /K "SET PATH=%CD%;%PATH%"' > pedstats-$(PEDVERSION)/STARTHERE.bat
	zip -r Windows-pedstats.zip pedstats-$(PEDVERSION)
	rm -rf pedstats-$(PEDVERSION)

.SUFFIXES : .cpp .c .o $(SUFFIXES)
