#ROM_SOLVER_OPTIONS
-tem_ksp_type bcgs 
-tem_ksp_rtol 1e-8
-tem_pc_type asm 
-tem_pc_asm_overlap 3 
-tem_mat_type mpibaij  
-tem_ksp_pc_side right 
-tem_sub_pc_type ilu
-tem_sub_pc_factor_mat_ordering rcm 
-tem_sub_pc_factor_levels 10 
-tem_svd_tol 1e-8
-tem_svd_max_it 1000
-tem_svd_type trlanczos
-tem_svd_nsv 100
-tem_bv_type vecs
#END_ROM_SOLVER_OPTIONS
