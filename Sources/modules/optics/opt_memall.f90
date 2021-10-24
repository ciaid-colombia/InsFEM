subroutine opt_memall(a)
  !-----------------------------------------------------------------------
  ! DESCRIPTION
  !    This routine allocates memory for the arrays 
  !-----------------------------------------------------------------------
  use typre
  use Mod_Memor
  use Mod_Listen
  use Mod_Mesh
  use Mod_PhysicalProblem
  use Mod_Optics
  implicit none
  
  class(OpticsProblem) :: a
  
  integer(ip) :: ndime,npoin,nelem

  
  !Unknowns
  call a%Mesh%GetNdime(ndime)
  call a%Mesh%GetNpoin(npoin)
  call a%Mesh%GetNelem(nelem)
  
  call a%Memor%alloc(npoin,a%cn2,'cn2','opt_memall')
  call a%Memor%alloc(npoin,a%cn2u53,'cn2u53','opt_memall')
  call a%Memor%alloc(npoin,a%ct2,'ct2','opt_memall')
  
  call a%Memor%alloc(npoin,a%avg_press,'avg_press','opt_memall')
  call a%Memor%alloc(npoin,a%avg_tempe,'avg_tempe','opt_memall')
  call a%Memor%alloc(npoin,a%avg_vnorm,'avg_vnorm','opt_memall')
  call a%Memor%alloc(npoin,a%avg_vdiss,'avg_vdiss','opt_memall')
  call a%Memor%alloc(npoin,a%avg_tdiss,'avg_tdiss','opt_memall')
  call a%Memor%alloc(npoin,2,a%avg_cn2,'avg_cn2','opt_memall')
  call a%Memor%alloc(npoin,2,a%avg_cn2u53,'avg_cn2u53','opt_memall')
  
  !If using Smagorinsky or WALE, we need to average the dissipations over time
  !We will not use the external dissipations, so we use them for exporting to nodes
  if (a%kfl_dissi /= 1) then
     call a%Memor%palloc(npoin,a%NSDissipation,'NSDissipation','opt_memall')
     call a%Memor%palloc(npoin,a%TempeDissipation,'TempeDissipation','opt_memall') 
  endif    
   
   
end subroutine opt_memall

