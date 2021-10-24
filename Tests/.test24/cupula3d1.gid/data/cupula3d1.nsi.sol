#NAVIERSTOKES MODULE PREFIX: nsi_

-nsi_ksp_type bcgs 
-nsi_ksp_pc_side right 
-nsi_ksp_rtol 1e-8 
-nsi_pc_type asm 
-nsi_pc_asm_overlap 3 
-nsi_sub_pc_type ilu 
-nsi_sub_pc_factor_mat_ordering rcm 
-nsi_sub_pc_factor_levels 10 
-nsi_mat_type mpibaij  

#FRACTIONAL STEP MODULE, VELOCITY system PREFIX: nsf_

-nsf_ksp_type bcgs 
-nsf_ksp_pc_side right 
-nsf_ksp_rtol 1e-12 
-nsf_pc_type asm 
-nsf_pc_asm_overlap 3 
-nsf_sub_pc_type ilu 
-nsf_sub_pc_factor_mat_ordering rcm 
-nsf_sub_pc_factor_levels 10 
-nsf_mat_type mpibaij  

#FRACTIONAL STEP MODULE, PRESSURE system PREFIX: nsf_

-nsp_ksp_type bcgs 
-nsp_ksp_pc_side right 
-nsp_ksp_rtol 1e-12 
-nsp_pc_type asm 
-nsp_pc_asm_overlap 3 
-nsp_sub_pc_type ilu 
-nsp_sub_pc_factor_mat_ordering rcm 
-nsp_sub_pc_factor_levels 10 
-nsp_mat_type mpibaij  
