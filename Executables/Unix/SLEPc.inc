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
ifneq ($(SLEPC_DIR),)
   SLEPC_FLAGS:= -DSLEPC
   PETSC_ARCH := $(PETSC_ARCH_DEBUG)
   PETSC_NEW := $(shell if [ $${PETSC_VERSION_MINOR:-5} -gt 5 ]; then echo true; fi)
   ifneq ($(PETSC_NEW),true)
      include ${SLEPC_DIR}/conf/slepc_variables
   else 
      include ${SLEPC_DIR}/lib/slepc/conf/slepc_variables
   endif
   INCLUDES_SG:= -I$(SLEPC_DIR)/include -I$(SLEPC_DIR)/$(PETSC_ARCH)/include 
   #-I$(SLEPC_DIR)/include/slepc
   LIBSSG := -L$(SLEPC_DIR)/$(PETSC_ARCH)/lib ${SLEPC_LIB}
   PETSC_ARCH := $(PETSC_ARCH_RELEASE)
   ifneq ($(PETSC_NEW),true)
      include ${SLEPC_DIR}/conf/slepc_variables
   else
      include ${SLEPC_DIR}/lib/slepc/conf/slepc_variables
   endif
   INCLUDES_SO:= -I$(SLEPC_DIR)/include -I$(SLEPC_DIR)/$(PETSC_ARCH)/include 
   #-I$(SLEPC_DIR)/include/slepc
   INCLUDES_SP:= -I$(SLEPC_DIR)/include -I$(SLEPC_DIR)/$(PETSC_ARCH)/include 
   #-I$(SLEPC_DIR)/include/slepc
   LIBSSO := ${SLEPC_LIB}
   LIBSSP := -L$(SLEPC_DIR)/$(PETSC_ARCH)/lib ${SLEPC_LIB} ${SLEPC_LIB}
else
   SLEPC_FLAGS:=
   INCLUDES_SG:=
   LIBSSG :=
   INCLUDES_SO:=
   INCLUDES_SP:=
   LIBSSO :=
   LIBSSP :=
endif 
