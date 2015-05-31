
# to build, create a link to a directory with CCP4 installation that  includes
# distribution libraries and header files
CCP4=./CCP4

ICLPR=${CCP4}/include/clipper
ICCP4=${CCP4}/include
IMMDB=${CCP4}/include/mmdb
LCCP4=${CCP4}/lib
LCLPR=${CCP4}/lib
CXX=g++
# -ltcl8.4
DFLIBS=    -lclipper-contrib -lclipper-ccp4 -lclipper-minimol -lclipper-mmdb -lclipper-phs -lclipper-cif -lclipper-core -lrfftw -lfftw -lmmdb2 -lccp4c    -lm -lc -lgsl
# -pg
#OPTS = -g 
OPTS = -O2   -g
CPPFLAGS= ${OPTS}  ${XCPPFLAGS} -I${ICCP4}


# FIXME execution requires access to CCP4 shared libraries: for now, handle via rpath
DF_LDFLAGS=  -L${LCCP4} -Wl,-rpath,${CCP4}/lib

%:		%.cpp %.o
		${CXX} $(CFLAGS) $(CPPFLAGS) -c $@ $< 

DFGEN_OBJS= dfg.o elemtab.o sphere.o

TARGS = asupad dfg dfsel

all: $(TARGS)

clean: 
	rm -f *.o $(TARGS) 

asupad:		asupad.o
	${CXX} $(CFLAGS) $(CPPFLAGS) $(DF_LDFLAGS)  $< ${DFLIBS} -o $@

dfg:		$(DFGEN_OBJS)
	${CXX} $(CFLAGS) $(CPPFLAGS) $(DF_LDFLAGS)  $(DFGEN_OBJS) ${DFLIBS} -o $@

dfsel:		dfsel.o
	${CXX} $(CFLAGS) $(CPPFLAGS) $(DF_LDFLAGS)  $< ${DFLIBS} -o $@

