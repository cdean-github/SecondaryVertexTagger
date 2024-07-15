AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include \
  -DHomogeneousField

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64 \
  `root-config --libs`

pkginclude_HEADERS = \
  SecondaryVertexTagger.h \
  SecondaryVertexTagger_engine.h \
  SecondaryVertexTagger_nTuple.h

lib_LTLIBRARIES = \
  libsecondaryvertextagger.la

libsecondaryvertextagger_la_SOURCES = \
  SecondaryVertexTagger.cc \
  SecondaryVertexTagger_engine.cc \
  SecondaryVertexTagger_nTuple.cc

libsecondaryvertextagger_la_LIBADD = \
  -lKFParticle \
  -lphool \
  -lphg4hit \
  -lSubsysReco \
  -lffamodules \
  -ltrackbase_historic_io \
  -lglobalvertex \
  -lg4eval

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libsecondaryvertextagger.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
