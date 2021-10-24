#############################################################
#                                                           #
#        SYSTEM DEFINITION DEPENDENT OPTIONS                #
#                                                           #
# These options are used to build VTK post-process          #
#                                                           #
#############################################################
export
ifneq ($(VTK_DIR),)
   VTK_FLAGS:= -DVTK
   VTK_O := $(VTKPOST_DIR)/CMakeFiles/vtkpost.dir/vtkc.cpp.o
   INCLUDES_VTK:= -I$(VTK_DIR)/include/vtk-$(VTK_VERSION) -I$(VTK_DIR)/include
   LIBS_VTK := -rdynamic $(VTK_DIR)/lib/libvtkFiltersGeneric-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkFiltersHyperTree-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkFiltersParallelGeometry-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkFiltersParallelImaging-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkFiltersParallelMPI-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkFiltersPoints-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkFiltersProgrammable-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkFiltersSMP-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkFiltersSelection-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkFiltersTexture-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkIOMPIParallel-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkIOParallelXML-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtklibxml2-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkIOParallel-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkFiltersParallel-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkIOGeometry-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkjsoncpp-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkParallelMPI-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkhdf5_hl-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkhdf5-$(VTK_VERSION).so.1 -lm $(VTK_DIR)/lib/libvtkParallelCore-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkIOLegacy-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkFiltersImaging-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkInfovisLayout-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkFiltersExtraction-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkIOXML-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkIOXMLParser-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkIOCore-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtklz4-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkexpat-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkRenderingFreeType-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkRenderingCore-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkCommonColor-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkFiltersGeometry-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkFiltersSources-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkFiltersGeneral-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkCommonComputationalGeometry-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkFiltersCore-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkCommonExecutionModel-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkCommonDataModel-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkCommonMisc-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkCommonSystem-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtksys-$(VTK_VERSION).so.1 -ldl $(VTK_DIR)/lib/libvtkCommonTransforms-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkCommonMath-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkCommonCore-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkfreetype-$(VTK_VERSION).so.1 $(VTK_DIR)/lib/libvtkzlib-$(VTK_VERSION).so.1 -Wl,-rpath,$(VTK_DIR)/lib 
else
   VTK_FLAGS := 
   VTK_O := 
   INCLUDES_VTK :=
   LIBS_VTK := 
endif 
