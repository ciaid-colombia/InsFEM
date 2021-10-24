
#####################################################
#LOWMACH MODULE PREFIX: lmn_
#####################################################

########################## DEFAULT ITERATIVE SOLVER
-lmn_ksp_type bcgs 
-lmn_ksp_rtol 1e-8 
-lmn_pc_type asm 
-lmn_pc_asm_overlap 3 
-lmn_mat_type mpibaij  
-lmn_ksp_pc_side right 
-lmn_sub_pc_type ilu 
-lmn_sub_pc_factor_mat_ordering rcm 
-lmn_sub_pc_factor_levels 10 

