
#####################################################
#NAVIERSTOKES MODULE PREFIX: nsc_
#####################################################

########################## ALGEBRAIC FRACTIONAL STEP USING PETSC FIELDSPLIT
#IMPORTANT: USE WITH SPLIT OSS
#-nsc_ksp_type richardson
#-nsc_ksp_max_it 2  #2 ITERATIONS IS THE MINIMUM TO OBTAIN 2ND ORDER FRACTIONAL STEP
#-nsc_pc_type fieldsplit
#-nsc_pc_fieldsplit_type schur
#-nsc_pc_fieldsplit_schur_factorization_type full
#-nsc_pc_fieldsplit_0_fields 0
#-nsc_pc_fieldsplit_1_fields 1,2
#-nsc_pc_fieldsplit_2_fields 3
#
##-nsc_fieldsplit_0_ksp_type preonly
##-nsc_fieldsplit_0_pc_type lu
#
#-nsc_pc_fieldsplit_schur_precondition selfp
#-nsc_fieldsplit_0_mat_schur_complement_ainv_type lump
#-nsc_fieldsplit_3_mat_schur_complement_ainv_type lump
#
#-nsc_fieldsplit_1_ksp_type preonly
#-nsc_fieldsplit_1_pc_type ksp
##-nsc_fieldsplit_1_ksp_ksp_type preonly
##-nsc_fieldsplit_1_ksp_pc_type lu
#
##-nsc_fieldsplit_0_ksp_converged_reason
##-nsc_fieldsplit_1_ksp_converged_reason
##-nsc_fieldsplit_1_ksp_ksp_converged_reason
##-nsc_ksp_converged_reason
##-nsc_ksp_view
############################################################################

#################################### HYPRE
#-nsc_mat_type mpiaij  
#-nsc_ksp_type bcgs
#-nsc_ksp_rtol 1e-8 
#-nsc_pc_type hypre
#-nsc_ksp_diagonal_scale


########################## DIRECT SOLVER USING MUMPS
#-nsc_mat_type mpibaij
#-nsc_ksp_type preonly 
#-nsc_pc_type lu
#-nsc_pc_factor_mat_solver_package mumps

########################## DEFAULT ITERATIVE SOLVER
-nsc_ksp_type bcgs 
-nsc_ksp_rtol 1e-8 
-nsc_pc_type asm 
-nsc_pc_asm_overlap 3 
-nsc_mat_type mpibaij  
#-nsc_ksp_pc_side right 
-nsc_sub_pc_type ilu 
-nsc_sub_pc_factor_mat_ordering rcm 
-nsc_sub_pc_factor_levels 10 

########################## TRILINOS ML SOLVER
#-nsc_mat_type mpiaij  
#-nsc_ksp_type bcgs
#-nsc_ksp_rtol 1e-8 
#-nsc_pc_type ml
#-nsc_pc_mg_cycles 2
#-nsc_pc_mg_smoothup 5 
#-nsc_pc_mg_smoothdown 1
#-nsc_pc_mg_type multiplicative

