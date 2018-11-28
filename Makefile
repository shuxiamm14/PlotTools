# --------------------------------------------------------------
# GNUmakefile
# --------------------------------------------------------------

ROOTCFLAGS      := $(shell root-config --cflags) -Iinclude -Iinclude/TtHFitter -Iinclude/fcnc -Iinclude/atlasstyle
ROOTLIBS        := $(shell root-config --libs)
ROOTGLIBS       := $(shell root-config --glibs)

EXTRALIBS +=$(ROOTLIBS) -L./lib
EXTRALIBS +=$(ROOTGLIBS) -lMinuit -lTMVA -lHistFactory -lRooStats -lRooFit -lRooFitCore

FCNCLIB := lib

CPPFLAGS += $(ROOTCFLAGS) -g -I.
CPPFLAGS  += -Wno-long-long -w 

CXX := clang++
MAKESHARED = clang++ -shared -fPIC -dynamiclib -single_module -O2 -mmacosx-version-min=10.10 -m64 -Wl,-install_name,@rpath/$(patsubst $(FCNCLIB)/%,%,$@)
ATLASOBJS 		= $(patsubst src/atlasstyle/Atlas%.C,bin/.Atlas%.o,$(wildcard src/atlasstyle/Atlas*.C))

PLOTOBJS		= $(patsubst src/fcnc/%.cc,bin/.%.o,$(wildcard src/fcnc/*.cc))

all: makebin $(FCNCLIB)/libAtlasStyle.so $(FCNCLIB)/libPlotTool.so 

makebin:
	@echo using compiler: $(CXX)
	@mkdir -p ./bin
	@mkdir -p ./lib
	@echo current directory: $(PWD)

$(FCNCLIB)/libAtlasStyle.so: $(ATLASOBJS)
	@echo Linking $@
	@$(MAKESHARED) $(CPPFLAGS) $(ROOTGLIBS) -o $@ $(ATLASOBJS)

$(FCNCLIB)/libPlotTool.so: $(PLOTOBJS) | $(FCNCLIB)/libAtlasStyle.so 
	@echo Linking $@ with $(PLOTOBJS)
	@$(MAKESHARED) $(CPPFLAGS) $(EXTRALIBS) -lAtlasStyle -o $@ $(PLOTOBJS)

bin/.%.o: src/fcnc/%.cc include/fcnc/%.h
	@echo Compiling $@
	@$(CXX) $(CPPFLAGS) -c $< -o $@

bin/.%.o: src/atlasstyle/%.C include/atlasstyle/%.h
	@echo Compiling $@
	@$(CXX) $(CPPFLAGS) -c $< -o $@

.PHONY:clean

clean:
	@rm bin/.* lib/*
