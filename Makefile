# --------------------------------------------------------------
# GNUmakefile
# --------------------------------------------------------------
INCLUDE := -Iinclude -Iinclude/TtHFitter -Iinclude/fcnc -Iinclude/atlasstyle
ROOTCFLAGS      := $(shell root-config --cflags) $(INCLUDE)
ROOTLIBS := $(shell root-config --libs)
ROOTGLIBS:= $(shell root-config --glibs)

EXTRALIBS +=$(ROOTLIBS) -L./lib
EXTRALIBS +=$(ROOTGLIBS) -lMinuit -lTMVA

FCNCLIB := lib

CPPFLAGS += $(ROOTCFLAGS) -g -I.
CPPFLAGS  += -Wno-long-long -w 

CXX := clang++
MAKESHARED = clang++ -shared -fPIC -dynamiclib -single_module -O2 -mmacosx-version-min=10.10 -m64 -Wl,-install_name,@rpath/$(patsubst $(FCNCLIB)/%,%,$@)
ATLASOBJS 		= $(patsubst src/atlasstyle/Atlas%.C,bin/.Atlas%.o,$(wildcard src/atlasstyle/Atlas*.C))
ATLASSRC		= $(wildcard src/atlasstyle/Atlas*.C)
ATLASINC		= $(wildcard include/atlasstyle/Atlas*.h)

PLOTOBJS		= $(patsubst src/fcnc/%.cc,bin/.%.o,$(wildcard src/fcnc/*.cc))

all: makebin $(FCNCLIB)/libAtlasStyle.so $(FCNCLIB)/libPlotTool.so 

makebin:
	@echo using compiler: $(CXX)
	@mkdir -p ./bin
	@mkdir -p ./lib
	@echo current directory: $(PWD)

bin/.atlasstyle_dict.cc: $(ATLASINC) include/atlasstyle/LinkDef.h
	@rootcint -f $@ -c $^

bin/.plotTool_dict.cc: include/fcnc/HISTFITTER.h include/fcnc/fcnc_include.h include/fcnc/histSaver.h include/fcnc/LinkDef.h       
	@rootcint -f $@ -c $(INCLUDE) $^

bin/.%_dict.o: bin/.%_dict.cc
	@$(CXX) $(CPPFLAGS) -c $< -o $@
	@mv $(patsubst bin/.%_dict.o, bin/.%_dict_rdict.pcm, $@) lib/.

$(FCNCLIB)/libAtlasStyle.so: $(ATLASOBJS) bin/.atlasstyle_dict.o
	@echo Linking $@
	@$(MAKESHARED) $(CPPFLAGS) $(ROOTGLIBS) -o $@ $^

$(FCNCLIB)/libPlotTool.so: $(PLOTOBJS) bin/.plotTool_dict.o | $(FCNCLIB)/libAtlasStyle.so 
	@echo Linking $@ with $^
	@$(MAKESHARED) $(CPPFLAGS) $(EXTRALIBS) -lAtlasStyle -o $@ $^

bin/.%.o: src/fcnc/%.cc include/fcnc/%.h
	@echo Compiling $@
	@$(CXX) $(CPPFLAGS) -c $< -o $@

bin/.%.o: src/atlasstyle/%.C include/atlasstyle/%.h
	@echo Compiling $@
	@$(CXX) $(CPPFLAGS) -c $< -o $@

.PHONY:clean

clean:
	@rm -r bin lib
