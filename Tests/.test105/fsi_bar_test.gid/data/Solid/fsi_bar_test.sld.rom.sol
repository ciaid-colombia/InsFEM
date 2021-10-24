#ROM_SOLVER_OPTIONS
-sld_ksp_type bcgs                                   
-sld_ksp_rtol 1e-8                                   
-sld_pc_type asm    
-sld_pc_asm_overlap 3                                
-sld_mat_type mpibaij                                
-sld_ksp_pc_side right                               
-sld_sub_pc_type ilu                                 
-sld_sub_pc_factor_mat_ordering rcm                  
-sld_sub_pc_factor_levels 10                         
-sld_svd_tol 1e-8
-sld_svd_max_it 1000
-sld_svd_type trlanczos                              
-sld_svd_nsv 100
-sld_bv_type vecs 
#END_ROM_SOLVER_OPTIONS
