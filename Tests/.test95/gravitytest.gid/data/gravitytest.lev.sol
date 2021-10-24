#LEVELSET MODULE PREFIX: lev_

-lev_ksp_type bcgs 
#-lev_ksp_pc_side right 
-lev_ksp_rtol 1e-12 
-lev_pc_type asm 
#-lev_pc_asm_overlap 3 
#-lev_sub_pc_type ilu 
#-lev_sub_pc_factor_mat_ordering rcm 
#-lev_sub_pc_factor_levels 10 
-lev_mat_type mpibaij  

#-lev_mat_type mpiaij  
#-lev_ksp_type bcgs
#-lev_ksp_rtol 1e-8 
#-lev_pc_type ml
#-lev_pc_mg_cycles 2
#-lev_pc_mg_smoothup 5 
#-lev_pc_mg_smoothdown 1
#-lev_pc_mg_type multiplicative

