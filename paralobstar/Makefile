# Set compiler
CXX := mpic++

TRGTDIR := bin
TARGET := paralobstar

SRCEXT := cpp
OBJEXT := o
HEADEXT := h
DEPEXT := d

# global include
INCLDIR := include

SRCDIR := src
OBJDIR := obj
DEPDIR := dep

SOURCES := $(wildcard $(SRCDIR)/*.$(SRCEXT))
DEPENDS := $(SOURCES:$(SRCDIR)/%.$(SRCEXT)=$(DEPDIR)/%.$(DEPEXT))
OBJECTS := $(SOURCES:$(SRCDIR)/%.$(SRCEXT)=$(OBJDIR)/%.$(OBJEXT))

# main source file stem name
MAIN := main

# setting general flags
CXXFLAGS := -std=c++11 -I$(INCLDIR) -I/usr/include/hdf5/openmpi
# setting linker flags
LDFLAGS := -std=c++11 -lhdf5 -lboost_filesystem -lboost_mpi -lboost_serialization

# get kernel name to be able to run sed correctly on Darwin (MacOS) or Linux kernels
KERNEL := $(shell uname -s)
ifeq ($(KERNEL), Darwin) 
	SED := sed -i "~"
else
	SED := sed -i
endif

all: $(OBJECTS) $(TRGTDIR)/$(TARGET)

# build all with debug flags
debug: CXXFLAGS += -g -Wall
# show linker invocation when building debug target
debug: LDFLAGS += -v
debug: all

clean: | $(TRGTDIR)
	@rm -r $(TRGTDIR)/$(TRGT) $(OBJECTS)

cleaner: clean | $(TRGTDIR)
	@rm -r $(DEPENDS)

$(TRGTDIR):
	@mkdir -p $@

$(DEPDIR):
	@mkdir -p $@

$(OBJDIR):
	@mkdir -p $@

# generating dependency files for all sources
# sed changes '%.o: ...' to '%.o %.d: ...' in dependency file
$(DEPDIR)/%.$(DEPEXT): $(SRCDIR)/%.$(SRCEXT) | $(DEPDIR)
	@echo "Generating dependency file '$@' ..."
	@$(CXX) -MM $(CXXFLAGS) $< -MF $@
	@$(SED) 's,$(*F).$(OBJEXT),$*.$(OBJEXT) $@,' $@
	@rm -f $@~
	@echo "... done."

# include targets from generated dependency files
include $(DEPENDS)

# build module as shared library and install it on your system
install:
	@echo "Not implemented yet :("

# build main target
# check if target directory 'bin' already exist via prerequisite
$(TRGTDIR)/$(TARGET): $(OBJDIR)/$(MAIN).$(OBJEXT) $(OBJECTS) | $(TRGTDIR)
	@echo "Linking binary '$@'..."
	@$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJECTS) -o $(TRGTDIR)/$(TARGET)
	@echo "... done."

# pattern rule to build object from source file
$(OBJDIR)/%.$(OBJEXT): $(DEPDIR)/%.$(DEPEXT) | $(OBJDIR)
	@echo "Compiling '$@' ..."
	@$(CXX) -c $(CXXFLAGS) $(@:$(OBJDIR)/%.$(OBJEXT)=$(SRCDIR)/%.$(SRCEXT)) -o $@
	@echo "... done."

.PHONY: all debug clean cleaner install
