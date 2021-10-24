#####################################################
#####################################################


#-plcd_mat_type mpiaij
#-plcd_ksp_type preonly
#-plcd_pc_type lu
#-plcd_pc_factor_mat_solver_package mumps
#-plcd_ksp_diagonal_scale


#-plcd_mat_type mpiaij
#-plcd_ksp_type preonly
#-plcd_pc_type lu
#-plcd_pc_factor_mat_solver_package mumps
#-plcd_ksp_diagonal_scale


-plcd_ksp_type bcgs 
-plcd_ksp_pc_side right 
-plcd_ksp_rtol 1e-8 
-plcd_ksp_max_it 2
#-plcd_pc_type asm 
#-plcd_ksp_diagonal_scale
#-plcd_pc_asm_overlap 3 
#-plcd_sub_pc_type ilu 
#-plcd_sub_pc_factor_mat_ordering rcm 
#-plcd_sub_pc_factor_levels 10 
#-plcd_mat_type mpibaij  

#-plcd_mat_type mpiaij  
#-plcd_ksp_type bcgs
#-plcd_ksp_rtol 1e-8 
#-plcd_pc_type ml
#-plcd_pc_mg_cycles 2
#-plcd_pc_mg_smoothup 5 
#-plcd_pc_mg_smoothdown 1
#-plcd_pc_mg_type multiplicative
#-plcd_ksp_diagonal_scale

