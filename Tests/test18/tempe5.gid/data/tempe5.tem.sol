#TEMPERATURE MODULE PREFIX: tem_ moco

-tem_ksp_type bcgs 
-tem_ksp_pc_side right 
-tem_ksp_rtol 1e-8 
-tem_pc_type asm 
-tem_pc_asm_overlap 3 
-tem_sub_pc_type ilu 
-tem_sub_pc_factor_mat_ordering rcm 
-tem_sub_pc_factor_levels 10 
-tem_mat_type mpibaij  
