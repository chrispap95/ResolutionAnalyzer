# name of the library
LIBNAME = ResolutionAnalyzer

#Necessary to use shell built-in commands
SHELL=bash
BOOSTSYS=/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/boost/1.67.0

# Define include paths
USERINCLUDES += -I$(ROOTSYS)/include/
USERINCLUDES += -I$(ROOFITSYS)/include/
USERINCLUDES += -isystem $(BOOSTSYS)/include

# Define libraries to link
USERLIBS += $(shell root-config --glibs) -lRooFit -lGenVector #-lTreePlayer -lTMVA
USERLIBS += -Wl,-rpath
USERLIBS += -L$(BOOSTSYS)/lib -lboost_regex -lboost_program_options -lboost_filesystem

CXXFLAGS = -Wall -W -O2 -std=c++1z
LDFLAGS = -shared -Wall -W


# If possible we'll use the clang compiler, it's faster and gives more helpful error messages
# If it's not available, then fallback to gcc.
CXX=g++
LD=g++
# CLANGPATH := $(shell type -p clang++)
# ifeq ($(CLANGPATH),)
# $(warning clang++ not found, reverting to g++!)
# CXX=g++
# LD=g++
# endif

CXXFLAGS += $(USERINCLUDES)
LIBS += $(USERLIBS)

# A list of directories
BASEDIR = $(shell pwd)
LIBDIR = $(BASEDIR)/lib
EXEDIR = $(BASEDIR)/bin
SRCDIR = $(BASEDIR)/src
OBJDIR = $(BASEDIR)/obj
TESTDIR = $(BASEDIR)/test
DOCDIR= $(BASEDIR)/docs
OBJ_EXT=o
TEST_EXT=cpp

# Build a list of srcs and bins to build
SRCS=$(wildcard $(BASEDIR)/src/*.cc)
EXES=$(wildcard $(BASEDIR)/test/*.cpp)
OBJS=$(subst $(SRCDIR), $(OBJDIR),$(subst cc,$(OBJ_EXT),$(SRCS)))

BINS=$(EXEDIR)/simpleBH

.PHONY: all
all: lib $(BINS)

docs: all
	doxygen Doxyfile

#$(EXEDIR)/%:  $(TESTDIR)/%.cpp $(LIBDIR)/lib$(LIBNAME).so $(wildcard $(BASEDIR)/include/*.h*)
#	$(CXX) -o $@ $(CXXFLAGS) $< $(LIBS) -L$(LIBDIR) -l$(LIBNAME)

$(EXEDIR)/simpleBH:  $(TESTDIR)/simpleBH.cpp $(LIBDIR)/lib$(LIBNAME).so $(wildcard $(BASEDIR)/include/*.h*)
	$(CXX) -o $@ $(CXXFLAGS) $< $(LIBS) -L$(LIBDIR) -l$(LIBNAME)

$(EXEDIR)/simpleBH_barebone:  $(TESTDIR)/simpleBH_barebone.cpp $(LIBDIR)/lib$(LIBNAME).so $(wildcard $(BASEDIR)/include/*.h*)
	$(CXX) -o $@ $(CXXFLAGS) $< $(LIBS) -L$(LIBDIR) -l$(LIBNAME)

$(OBJDIR)/%.$(OBJ_EXT):  $(SRCDIR)/%.cc $(BASEDIR)/include/%.h*
	$(CXX) $(CXXFLAGS) -fPIC -c $<  -o $@

$(LIBDIR)/lib$(LIBNAME).so:  $(OBJS) $(LIBDEPENDS)
	$(LD) $(LDFLAGS) -o $(LIBDIR)/lib$(LIBNAME).so $(OBJS) $(LIBS)

lib: $(LIBDIR)/lib$(LIBNAME).so


#need to put stuff into userlib or it doesn't work....
#dictionary:
#	rootcint -v -f src/dict.cc -c include/MipHit.hh include/LinkDef.h
#	sed "s/include\/Mip/Mip/" src/dict.h > include/dict.h
#	rm src/dict.h



# info:
# 	@echo "LIBS: " $(LIBS)
# 	@echo "CXXFLAGS: " $(CXXFLAGS)
# 	@echo "Source files: " $(SRCS)
# 	@echo "Object files: " $(OBJS)
# 	@echo "Executables:  " $(TARGETS)

clean:
	rm -rf $(OBJS) $(LIBDIR)/lib$(LIBNAME).so $(BINS)

userclean:
	$(MAKE) -C userlib/ clean
