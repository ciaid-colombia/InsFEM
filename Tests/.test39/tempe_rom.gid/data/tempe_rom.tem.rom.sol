#TEMPERATURE MODULE PREFIX: tem_

-tem_ksp_type preonly 
#-tem_ksp_pc_side right 
-tem_ksp_rtol 1e-8 
-tem_pc_type lu

-tem_svd_tol 1e-6
-tem_svd_max_it 1000
-tem_svd_type trlanczos
-tem_bv_type mat
