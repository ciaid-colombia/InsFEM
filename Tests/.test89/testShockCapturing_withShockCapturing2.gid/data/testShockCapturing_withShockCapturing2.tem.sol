#TEMPERATURE MODULE PREFIX: tem_

-tem_ksp_type preonly 
#-tem_ksp_pc_side right 
-tem_ksp_rtol 1e-8 
-tem_pc_type lu 
#-tem_pc_asm_overlap 3 
#-tem_sub_pc_type ilu 
#-tem_sub_pc_factor_mat_ordering rcm 
#-tem_sub_pc_factor_levels 10 
#-tem_mat_type mpibaij  

#-tem_mat_type mpiaij  
#-tem_ksp_type bcgs
#-tem_ksp_rtol 1e-8 
#-tem_pc_type ml
#-tem_pc_mg_cycles 2
#-tem_pc_mg_smoothup 5 
#-tem_pc_mg_smoothdown 1
#-tem_pc_mg_type multiplicative

