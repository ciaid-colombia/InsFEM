#####################################################
#MESH PREFIX: dom_
#####################################################

########################## DIRECT SOLVER USING MUMPS
-dom_ksp_type preonly 
-dom_pc_type asm 
-dom_sub_pc_type lu  
-dom_mat_type mpibaij 
########################## DEFAULT ITERATIVE SOLVER
#-dom_ksp_type bcgs 
#-dom_ksp_rtol 1e-8 
#-dom_pc_type asm 
#-dom_pc_asm_overlap 3 
#-dom_mat_type mpibaij  
#-dom_ksp_pc_side right 
#-dom_sub_pc_type ilu 
#-dom_sub_pc_factor_mat_ordering rcm 
#-dom_sub_pc_factor_levels 10 

########################## TRILINOS ML SOLVER
#-dom_mat_type mpiaij  
#-dom_ksp_type bcgs
#-dom_ksp_rtol 1e-8 
#-dom_pc_type ml
#-dom_pc_mg_cycles 2
#-dom_pc_mg_smoothup 5 
#-dom_pc_mg_smoothdown 1
#-dom_pc_mg_type multiplicative

