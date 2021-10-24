
#####################################################
#NAVIERSTOKES MODULE PREFIX: nsi_
#####################################################

########################## DIRECT SOLVER USING MUMPS
-nsi_mat_type aij
-nsi_ksp_type preonly 
-nsi_pc_type lu
#-nsi_pc_factor_mat_solver_package mumps

########################## DEFAULT ITERATIVE SOLVER
#-nsi_ksp_type bcgs 
#-nsi_ksp_rtol 1e-8 
#-nsi_pc_type asm 
#-nsi_pc_asm_overlap 3 
#-nsi_mat_type mpibaij  
#-nsi_ksp_pc_side right 
#-nsi_sub_pc_type ilu 
#-nsi_sub_pc_factor_mat_ordering rcm 
#-nsi_sub_pc_factor_levels 10 

########################## TRILINOS ML SOLVER
#-nsi_mat_type mpiaij  
#-nsi_ksp_type bcgs
#-nsi_ksp_rtol 1e-8 
#-nsi_pc_type ml
#-nsi_pc_mg_cycles 2
#-nsi_pc_mg_smoothup 5 
#-nsi_pc_mg_smoothdown 1
#-nsi_pc_mg_type multiplicative

######################################################
#FRACTIONAL STEP MODULE, VELOCITY system PREFIX: nsf_
######################################################

-nsf_ksp_type bcgs 
-nsf_mat_type mpibaij  
-nsf_ksp_rtol 1e-8 
-nsf_pc_type asm 
#-nsf_ksp_pc_side right 
#-nsf_pc_asm_overlap 3 
#-nsf_sub_pc_type ilu 
#-nsf_sub_pc_factor_mat_ordering rcm 
#-nsf_sub_pc_factor_levels 10 

######################################################
#FRACTIONAL STEP MODULE, PRESSURE system PREFIX: nsf_
######################################################

-nsp_ksp_type bcgs 
-nsp_ksp_rtol 1e-8 
-nsp_pc_type asm 
-nsp_pc_asm_overlap 3 
-nsp_mat_type mpibaij  
#-nsp_sub_pc_type ilu 
#-nsp_sub_pc_factor_mat_ordering rcm 
#-nsp_sub_pc_factor_levels 10 
#-nsp_ksp_pc_side right 

