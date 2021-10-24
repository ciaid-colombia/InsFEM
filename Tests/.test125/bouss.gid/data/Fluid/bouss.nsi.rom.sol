#ROM_SOLVER_OPTIONS
-nsi_ksp_type bcgs 
-nsi_ksp_rtol 1e-8
-nsi_pc_type asm 
-nsi_pc_asm_overlap 3 
-nsi_mat_type mpibaij  
-nsi_ksp_pc_side right 
-nsi_sub_pc_type ilu 
-nsi_sub_pc_factor_mat_ordering rcm 
-nsi_sub_pc_factor_levels 10 
-nsi_svd_tol 1e-16
-nsi_svd_max_it 1000000
-nsi_svd_type trlanczos
-nsi_svd_nsv 1500
-nsi_bv_type vecs
#END_ROM_SOLVER_OPTIONS
