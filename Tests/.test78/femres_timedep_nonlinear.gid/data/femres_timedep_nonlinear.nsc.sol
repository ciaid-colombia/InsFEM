#NSCOMPRESSIBLE MODULE PREFIX: nsc_

-nsc_ksp_type bcgs 
#-nsc_ksp_pc_side right 
-nsc_ksp_rtol 1e-8 
-nsc_pc_type asm 
#-nsc_pc_asm_overlap 3 
#-nsc_sub_pc_type ilu 
#-nsc_sub_pc_factor_mat_ordering rcm 
#-nsc_sub_pc_factor_levels 10 
-nsc_mat_type mpibaij  

#-nsc_mat_type mpiaij  
#-nsc_ksp_type bcgs
#-nsc_ksp_rtol 1e-8 
#-nsc_pc_type ml
#-nsc_pc_mg_cycles 2
#-nsc_pc_mg_smoothup 5 
#-nsc_pc_mg_smoothdown 1
#-nsc_pc_mg_type multiplicative

