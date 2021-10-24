
-sld_ksp_type bcgs 
-sld_ksp_pc_side right 
-sld_ksp_rtol 1e-8 
-sld_pc_type asm 
-sld_pc_asm_overlap 3 
-sld_sub_pc_type ilu 
-sld_sub_pc_factor_mat_ordering rcm 
-sld_sub_pc_factor_levels 10 
-sld_mat_type mpibaij  

#-sld_mat_type mpiaij  
#-sld_ksp_type bcgs
#-sld_ksp_rtol 1e-8 
#-sld_pc_type ml
#-sld_pc_mg_cycles 2
#-sld_pc_mg_smoothup 5 
#-sld_pc_mg_smoothdown 1
#-sld_pc_mg_type multiplicative

