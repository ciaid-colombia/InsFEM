#=======================================================================
# User definitions:
#=======================================================================
NAME=InspiraFEM
NPROCS:=$$(($(shell grep -c ^processor /proc/cpuinfo) -1))
#COMPILER=MPINTEL
#COMPILER=GNU
#COMPILER=MPINTEL
#=======================================================================
# System definitions for versions debug(G), release(O) and profile(P)
#=======================================================================
-include InspiraFEM.inc
-include SLEPc.inc 
-include VTK.inc 

# let us put the debug version as default. This works when doing first g and then O
	TARGETNAME=Debug
	VERSION=DEBUG
	DEPS=$(INCLUDES_G) $(INCLUDES_SG) $(INCLUDES_VTK)
	O=Objects_g$(COMPILER)
	FFLAGS=$(FFG) $(VTK_FLAGS) $(SLEPC_FLAGS) $(HDF5_FLAGS) 
	LFLAGS=$(LFG)
	LIBS=$(LIBSG) $(LIBSSG) $(LIBS_VTK)

ifeq "$(V)" "O"
	TARGETNAME=Release
	VERSION=RELEASE
	DEPS=$(INCLUDES_O) $(INCLUDES_SO) $(INCLUDES_VTK)
	O=Objects_$(V)$(COMPILER)
	FFLAGS=$(FFO) $(VTK_FLAGS) $(SLEPC_FLAGS) $(HDF5_FLAGS)
	LFLAGS=$(LFO)
	LIBS=$(LIBSO) $(LIBSSO) $(LIBS_VTK)

endif
ifeq "$(V)" "p"
	TARGETNAME=Profile
	VERSION=PROFILE
	DEPS=$(INCLUDES_P) $(INCLUDES_SP) $(INCLUDES_VTK)
	O=Objects_$(V)$(COMPILER)
	FFLAGS=$(FFP) $(VTK_FLAGS) $(SLEPC_FLAGS) $(HDF5_FLAGS)
	LFLAGS=$(LFP)
	LIBS=$(LIBSP) $(LIBSSP) $(LIBS_VTK)

endif
ifeq "$(V)" "v"
	TARGETNAME=Valgrind
	VERSION=xxx
	O=Objects_$(V)$(COMPILER)
	FFLAGS=$(FFV) $(VTK_FLAGS) $(SLEPC_FLAGS) $(HDF5_FLAGS)
	LFLAGS=$(LFV)
	LIBS=$(LIBSV) $(LIBS_VTK)

endif
ifeq "$(V)" "pg"
	TARGETNAME=Gprof
	VERSION=xxx
	O=Objects_$(V)$(COMPILER)
	FFLAGS=$(FFPG) $(VTK_FLAGS) $(SLEPC_FLAGS) $(HDF5_FLAGS)
	LFLAGS=$(LFPG)
	LIBS=$(LIBSPG) $(LIBS_VTK)

endif
#=======================================================================
# Objects and source directories:
#=======================================================================
-include InspiraFEM.conf.obj
#=======================================================================
# Variables:
#=======================================================================
vpath %f   ${SOURCE_DIRS}
vpath %f90 ${SOURCE_DIRS}
vpath %F90 ${SOURCE_DIRS}
#=======================================================================
# Rules:
#=======================================================================
.PHONY : default all debug release profile ctarget clean cleandebug cleanrelease cleanprofile cleanobj preprocess config test 

default: cleanlog config debug release test

cleanlog:
	 @rm -f *.log

all:     debug release profile

debug:
	$(MAKE) V=g ctarget

release:
	$(MAKE) V=O ctarget

profile:
	$(MAKE) V=p ctarget

ctarget: preprocess
	@echo '---------------------------->>> BEGIN $(TARGETNAME) Version'
	@rm -f $(NAME).$(V)$(COMPILER)
	@mkdir -p $(O)
	@rm -f $(NAME).$(V)$(COMPILER) #delete to avoid using it and thinking it's the new one
ifeq "$(serial)" "yes"	
	$(MAKE) $(NAME).$(V)
else
	$(MAKE) $(NAME).$(V) -j $(NPROCS)
endif	
	@chmod 770 $(NAME).$(V)$(COMPILER)
	@echo '---------------------------->>> END $(TARGETNAME) Version'
	@echo '  '
	@ln -s -f InspiraFEM.sh insFEM

clean: cleandebug cleanrelease cleanprofile cleanpreprocessor
	@rm -f *.log

cleanpreprocessor:
	@rm -rf ../../.SourcesAfterPrep
	@rm -rf ../../.SourcesIntel

cleandebug:
	$(MAKE) "O=Objects_g$(COMPILER)" "V=g" cleanobj

cleanrelease:
	$(MAKE) "O=Objects_O$(COMPILER)" "V=O" cleanobj

cleanprofile:
	$(MAKE) "O=Objects_p$(COMPILER)" "V=p" cleanobj

cleanobj:
	@rm -f $(OFILES)
	@rm -f $(NAME).$(V)$(COMPILER)
	@rm -f $(O)/*.o
	@rm -f $(O)/*.mod

$(NAME)_$V.a: $(OFILES)
	@xiar cr $(NAME)_$(V).a $(O)/*.o $(VTK_O) $(INCLUDES_VTK) $(LIBS)

$(NAME).$(V): $(OFILES)
	$(LINK) $(LFLAGS) -o $(NAME).$(V)$(COMPILER) $(O)/*.o $(VTK_O) $(INCLUDES_VTK) $(LIBS)

#=======================================================================
# Compilation rules
#=======================================================================
$O/%.o: %.f90
	$(FCOMP) $(FFLAGS) $(CODE_COMP_FLAGS) $(DEPS) -o $@ $< 
#@if [ -f $(NAME).opt.log ]; then mv $(NAME).opt.log $@.opt.log; fi

$O/%.o: %.f
	$(FCOMP) $(FFLAGS) $(CODE_COMP_FLAGS) $(DEPS) -o $@ $<
#@if [ -f $(NAME).opt.log ]; then mv $(NAME).opt.log $@.opt.log; fi

#=======================================================================
# Dependencies
#=======================================================================
-include InspiraFEM.conf.dep

#.PHONY: config InspiraFEM.conf.obj

config: preprocess
	@echo 'Building InspiraFEM object list and dependencies'
	@echo $(SRCDIR)
	@./configure InspiraFEM.conf --srcdir=../../$(SRCDIR) --objdir='$$O' --check=InspiraFEM

#=======================================================================
# Preprocessing
#=======================================================================
preprocess:
ifeq ($(PREPROCESS),DONE)
else
	@echo 'Preprocessing'
	PREPROCESS=DONE;
	@cd ../../Scripts/Preprocessor; ./Preprocessor.prep
ifeq ($(COMPILER),MPINTEL)
   ifneq ($(INTEL_NEW),true)
	   @cd ../../Scripts/Preprocessor; ./ToIntel.prep;
   endif
endif
endif

#=======================================================================
# Testing
#=======================================================================
test:
	@cd ../../Tests; ./Tester.prep test

