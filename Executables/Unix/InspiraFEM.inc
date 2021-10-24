#############################################################
#                                                           #
#        SYSTEM DEFINITION DEPENDENT OPTIONS                #
#                                                           #
# These options are used to build a Debug (G), Release (O)  #
# and Profile (P) versions of the code:                     #
#                                                           #
#    * Object and Executables directories (OBJDIR EXEDIR)   #
#    * Fortran compiler and its options (FCOMP,FFG,FFO,FFP) #
#    * C compiler and its options (CCOMP,CFG,CFO,CFP)       #
#    * Liner and its options (LINK,LFG,LFO,LFP)             #
#    * Libraries (LIBG,LIBO,LIBP)                           #
#                                                           #
#############################################################

export
   CODE_COMP_FLAGS =
   ifdef CODE_COMPILER_DIRECTIVES
      space :=
      space +=
      aux = ${CODE_COMPILER_DIRECTIVES}
      CODE_COMP_FLAGS = $(subst $(space), -D, $(aux))
   endif
   PETSC_ARCH := $(PETSC_ARCH_DEBUG)
   PETSC_NEW := $(shell if [ $${PETSC_VERSION_MINOR:-5} -gt 5 ]; then echo true; fi)
   ifneq ($(PETSC_NEW),true)
      include ${PETSC_DIR}/conf/variables
   else
      include ${PETSC_DIR}/lib/petsc/conf/variables
   endif
   INCLUDES_G:= -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include 
   LIBSG := -L$(PETSC_DIR)/$(PETSC_ARCH)/lib ${PETSC_LIB} 
   PETSC_ARCH := $(PETSC_ARCH_RELEASE)
   ifneq ($(PETSC_NEW),true)
      include ${PETSC_DIR}/conf/variables
   else 
      include ${PETSC_DIR}/lib/petsc/conf/variables
   endif
   INCLUDES_O:= -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include 
   #-I$(PETSC_DIR)/include/petsc
   INCLUDES_P:= -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include 
   #-I$(PETSC_DIR)/include/petsc
   LIBSO := ${PETSC_LIB} -L$(PETSC_DIR)/$(PETSC_ARCH)/lib
   LIBSP := -L$(PETSC_DIR)/$(PETSC_ARCH)/lib ${PETSC_LIB} ${PETSC_LIB}

ifneq ($(HDF5),)
   HDF5_FLAGS:= -DHDF5
else
   HDF5_FLAGS:=
endif   

ifeq ($(COMPILER),MPINTEL)
   IPO=
   #OMP = -Dompcode -openmp -openmp-report2
   #OMP = -openmp
   FCOMP = mpif90 -c -module $O -fpp $(OMP)
   INTEL_V=$(shell mpif90 --version | sed '3d' | sed '2d' | sed -r 's/.*(.{8})/\1/')
   #$(info $$INTEL_V is [${INTEL_V}])
   INTEL_NEW := $(shell if [ ${INTEL_V} -gt 20160000 ]; then echo true; fi)
   ifneq ($(INTEL_NEW),true)
      FFG = ${PETSC_FC_INCLUDES} -traceback -g -debug all -check all -nogen-interfaces -ftrapuv -check nooutput_conversion
      # -tcheck -check bounds
      FFO = ${PETSC_FC_INCLUDES} -O3 -no-inline-min-size -no-inline-max-size -no-inline-max-total-size -no-inline-max-per-routine -no-inline-max-per-compile -qopt-report-file=$(NAME).opt.log $(IPO) # xHost
   else
      FFG = -g -traceback -debug all -check all -g3
      # -tcheck -check bounds
      FFO = -O3 -no-inline-min-size -no-inline-max-size -no-inline-max-total-size -no-inline-max-per-routine -no-inline-max-per-compile -qopt-report-file=$(NAME).opt.log $(IPO) # xHost
   endif
   FFP = $(FFO) -g  # for vtune
   ifneq ($(INTEL_NEW),true)
      SRCDIR = .SourcesIntel 
   else
      SRCDIR = Sources
   endif
   LINK= mpif90
endif



#############################################################
# GNU Fortran compiler and linker
ifeq ($(COMPILER),MPIGNU)
   # Activate interprocedural optimization using IPO
   #IPO=-ipo
   IPO=
   FCOMP = mpif90 -c -ffree-line-length-0 -J$(O) -I$(O) -x f95-cpp-input
   #FCOMP += -fopenmp
   FFG  = -fbacktrace  -g -fcheck=bounds -gdwarf-2 -gno-strict-dwarf -g3 -finit-integer=-12 -Wuninitialized -finit-real=snan -finit-logical=false -finit-character=3
#FFG += -fcheck=mem
   #FFG += -fcheck=pointer
   #-fbounds-check  DEPRECATED
   ifneq ($(SKIP_FINIT_DERIVED),)
      FFG += 
   else
      FFG += -finit-derived
   endif
   FFO = -O3 -funroll-loops -ftree-vectorize -msse3 $(IPO) 
   FFP = -pg $(FFO)
   
   LINK=mpif90
   #-openmp
   #LFG = -g
   LFG = $(FFG)
   LFO = $(IPO)
   LFP = $(FFP)
   SRCDIR = Sources
endif

#############################################################
