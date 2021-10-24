#
#-nsi_ksp_type richardson
#-nsi_ksp_max_it 1
#-nsi_pc_type fieldsplit
#-nsi_pc_fieldsplit_type schur
#-nsi_pc_fieldsplit_schur_factorization_type full
#-nsi_pc_fieldsplit_0_fields 0,1
#-nsi_pc_fieldsplit_1_fields 2
#
#-nsi_fieldsplit_0_ksp_type preonly
#-nsi_fieldsplit_0_pc_type lu
#
#-nsi_pc_fieldsplit_schur_precondition selfp
#-nsi_fieldsplit_1_mat_schur_complement_ainv_type lump
#
#-nsi_fieldsplit_1_ksp_type preonly
#-nsi_fieldsplit_1_pc_type ksp
#-nsi_fieldsplit_1_ksp_ksp_type preonly
#-nsi_fieldsplit_1_ksp_pc_type lu
#
#-nsi_fieldsplit_0_ksp_converged_reason
#-nsi_fieldsplit_1_ksp_converged_reason
#-nsi_fieldsplit_1_ksp_ksp_converged_reason
#-nsi_ksp_converged_reason
##-nsi_ksp_view

-nsi_ksp_type preonly
-nsi_pc_type lu




