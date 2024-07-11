AUTOMAKE_OPTIONS = foreign

ACLOCAL_AMFLAGS = -I m4

lib_LTLIBRARIES = \
  libMbdQA.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I`root-config --incdir` #-isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib\
  -L$(OFFLINE_MAIN)/lib64 \
  -lcentrality_io \
  -lcalotrigger_io

include_HEADERS = \
  MbdQA.h 

ROOTDICTS = \
  MbdQA_Dict.C 

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  MbdQA_Dict_rdict.pcm 

libMbdQA_la_SOURCES = \
  $(ROOTDICTS) \
  MbdQA.cc 

libMbdQA_la_LIBADD = \
-lphool \
-lSubsysReco \
-lg4detectors \
-lmbd_io \
-lfun4allraw \
-lffarawobjects \
-lglobalvertex_io \
-lphg4hit

libMbdQA_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lfun4all \
  -lg4detectors_io \
  -lphg4hit \
  -lcentrality_io\
  `root-config --libs`

# Rule for generating CINT dictionaries from class headers.
%_Dict.C: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^


#just to get the dependency
%_Dict_rdict.pcm: %_Dict.C ;


##########################################
# linking tests
BUILT_SOURCES = \
  testexternals.C

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.C

testexternals_LDADD = \
  libMbdQA.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@


clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm



