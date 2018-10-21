AUTOMAKE_OPTIONS = foreign
MAJOR_VERSION=1
MINOR_VERSION=0
VERSION=$(MAJOR_VERSION).$(MINOR_VERSION)
LIBTOOL_VERSION=$(MAJOR_VERSION):$(MINOR_VERSION)

ACLOCAL_AMFLAGS= -I m4
AM_CPPFLAGS = -DGMP -Wall $(INCLUDES)
AM_CXXFLAGS = -DGMP -Wall $(INCLUDES)

# Support both static and dynamic builds
bin_PROGRAMS = corels
# ADD THIS TO PREVIOUS LINE FOR STATIC BUILD: corels_static
include_HEADERS = rule.h queue.hh pmap.hh alloc.hh cache.hh utils.hh
dist_pkgdata_DATA = data/tictactoe.rda data/tictactoe.txt.gz

# Dynamic build
corels_SOURCES = main.cc queue.hh pmap.hh alloc.hh cache.hh utils.hh
corels_LDADD = libcorels.la

lib_LTLIBRARIES = libcorels.la
libcorels_la_SOURCES = alloc.hh cache.cc cache.hh corels.cc pmap.cc pmap.hh queue.hh rule.h rulelib.c utils.cc utils.hh
libcorels_la_LDFLAGS = -version-info $(LIBTOOL_VERSION)
