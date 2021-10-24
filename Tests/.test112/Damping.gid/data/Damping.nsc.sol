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
