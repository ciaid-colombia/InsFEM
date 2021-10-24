#ROM_SOLVER_OPTIONS
-lmn_ksp_type bcgs                                   
-lmn_ksp_rtol 1e-8                                   
-lmn_pc_type asm    
-lmn_pc_asm_overlap 3                                
-lmn_mat_type mpibaij                                
-lmn_ksp_pc_side right                               
-lmn_sub_pc_type ilu                                 
-lmn_sub_pc_factor_mat_ordering rcm                  
-lmn_sub_pc_factor_levels 10                         
-lmn_svd_tol 1e-8
-lmn_svd_max_it 1000
-lmn_svd_type trlanczos                              
-lmn_svd_nsv 100
-lmn_bv_type vecs 
#END_ROM_SOLVER_OPTIONS
