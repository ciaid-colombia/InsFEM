subroutine opt_turnof(a)
   use typre
   use Mod_Optics
   use def_parame
   implicit none
   class(OpticsProblem) :: a
   
   integer(ip) :: npoin,npoinLocal,nboun,bvnat_coun,iboun,ifunc,istat
   real(rp)    :: cputim1,cputim2
   character(150) :: outstr
   
   interface
      subroutine opt_TurnofComputations(a)
         use Mod_Optics
         implicit none
         class(OpticsProblem) :: a
      end subroutine
   end interface
   
   
   call a%Timer%Total%Tic
   
   !Computations for mean values and rays
   call opt_TurnofComputations(a)
   
   
   call a%Timer%Output%Tic
   !Output results
   call a%Output(one)
   call a%Timer%Output%Toc
   
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNpoinLocal(npoinLocal)
   call a%Mesh%GetNboun(nboun)
   
   !Output
   outstr = adjustl(trim(a%exmod))//'_turnof'

   !MEMORY DEALLOCATION
   call a%Memor%dealloc(npoin,a%ct2,'ct2',outstr)
   call a%Memor%dealloc(npoin,a%cn2,'cn2',outstr)
   call a%Memor%dealloc(npoin,a%cn2u53,'cn2u53',outstr)
   if (allocated(a%beams)) then
      deallocate(a%beams,STAT=istat)
      call a%Memor%deallocObj(istat,'beams',outstr,a%nbeams*6*rp)
   endif
   
   call a%Memor%dealloc(npoin,a%avg_press,'avg_press','opt_memall')
   call a%Memor%dealloc(npoin,a%avg_tempe,'avg_tempe','opt_memall')
   call a%Memor%dealloc(npoin,a%avg_vnorm,'avg_vnorm','opt_memall')
   call a%Memor%dealloc(npoin,a%avg_vdiss,'avg_vdiss','opt_memall')
   call a%Memor%dealloc(npoin,a%avg_tdiss,'avg_tdiss','opt_memall')
   call a%Memor%dealloc(npoin,2,a%avg_cn2,'avg_cn2','opt_memall')
   call a%Memor%dealloc(npoin,2,a%avg_cn2u53,'avg_cn2u53','opt_memall')
   
   !If using Smagorinsky, we need to average the dissipations over time
   !We will not use the external dissipations, so we use them for exporting to nodes
   if (a%kfl_dissi /= 1) then
      call a%Memor%pdealloc(npoin,a%NSDissipation,'NSDissipation','opt_memall')
      call a%Memor%pdealloc(npoin,a%TempeDissipation,'TempeDissipation','opt_memall') 
   endif    
   
   !Close Files
   call a%CloseFiles
   
   call a%Timer%Total%Toc
      
 end subroutine