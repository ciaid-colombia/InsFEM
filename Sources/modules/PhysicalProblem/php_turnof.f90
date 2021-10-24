subroutine php_turnof(a)
   use typre
   use def_parame
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a
   
   integer(ip) :: npoin,npoinLocal,nboun,bvnat_coun,iboun,ifunc
   real(rp)    :: cputim1,cputim2
   character(150) :: outstr
   
   integer(ip) :: ndime
   
   call a%Timer%Total%Tic
   call a%Timer%Turnof%Tic
   
   !Output results
   call a%Timer%Output%Tic
   call a%Output(one)
   call a%Timer%Output%Toc
   
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNpoinLocal(npoinLocal)
   call a%Mesh%GetNboun(nboun)
   
   !Output
   outstr = adjustl(trim(a%exmod))//'_turnof'

   !MEMORY DEALLOCATION
   call a%Memor%dealloc(a%ndofbc,npoin,a%kfl_fixno,'kfl_fixno',outstr)
   call a%Memor%dealloc(nboun,a%kfl_fixbo,'kfl_fixbo',outstr)
   
   bvnat_coun = 0
   do iboun = 1,nboun
      if (associated(a%bvnat(iboun)%a)) then
         bvnat_coun = bvnat_coun + size(a%bvnat(iboun)%a)
         deallocate(a%bvnat(iboun)%a)
      endif
   enddo
   call a%Memor%deallocObj(0,'bvnat%a',outstr,bvnat_coun*rp)
   call a%Memor%dealloc(nboun,a%bvnat,'bvnat',outstr)
   
   if(a%kfl_conbc == 1) then
      call a%Memor%dealloc(a%ndofbc,npoin,1,a%bvess,'bvess',outstr)
   else
      call a%Memor%dealloc(npoin,a%kfl_funno,'kfl_funno',outstr)
      call a%Memor%dealloc(nboun,a%kfl_funbo,'kfl_funbo',outstr)
      do ifunc = 1,10
            if(associated(a%funpa(ifunc)%a)) call a%Memor%pdealloc(a%kfl_funty(ifunc,2),a%funpa(ifunc)%a,'funpa%a',outstr)
      enddo
      call a%Memor%dealloc(10,2 ,a%kfl_funty,'kfl_funty',outstr) 
      call a%Memor%dealloc(10   ,a%funpa,'funpa',outstr)
      call a%Memor%dealloc(a%ndofbc,npoin,2,a%bvess,'bvess',outstr)
   end if

   if(a%nptra>0) then
      call a%TrackingInterpolator%Finalize
      call a%Mesh%GetNdime(ndime)
      !if (a%MPIrank == a%MPIroot) then
         call a%Memor%dealloc(ndime,a%nptra,a%cptra,'cptra','php_reaous')
      !endif
   endif
   
   call a%LinearSystemTurnof

   call a%SpecificTurnof
   
   !Close Files
   call a%CloseFiles
   
   !call a%FilePostpr%closepos(a%MPIrank,a%MPIroot)

   call a%Timer%Turnof%Toc
   call a%Timer%Total%Toc
      
 end subroutine
