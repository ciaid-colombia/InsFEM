subroutine GeneralRestart(a,array,itask)
   !------------------------------------------------------------------------
   !    This routine reads the initial values from the restart file, 
   !    interpolates them if necessary and writes the new restart
   !------------------------------------------------------------------------
   use typre
   use def_parame
   use Mod_iofile
   use Mod_int2str
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a
   integer(ip), intent(in)    :: itask
   
   integer(ip) :: kfl_twost, kfl_twost_pre,kfl_trasg_pre, kfl_repro_pre
   integer(ip) :: ndime,npoin,oldNpoin,ipoin,idofn,ndofn
   real(rp)    :: array(:,:)
   real(rp), allocatable :: auxarray(:,:)
   
   call a%Timer%Total%Tic
   call a%Timer%Restar%Tic
  

   select case (itask)

   case (1)
      
      call a%Mesh%GetNdime(ndime)
      call a%Mesh%GetNpoin(npoin)
      ndofn = size(array,1)

      if (a%kfl_inter /= 0_ip) then
         call a%OldMesh%GetNpoin(oldNpoin)
      else
         oldNpoin = npoin
      endif
         
      call a%Memor%alloc(ndofn,oldNpoin,auxarray,'auxarray','php_GeneralRestart')

      ! Read the previous time scheme used: one integer number
      !read(a%lun_rsta2,*) kfl_twost, kfl_trasg_pre, kfl_repro_pre 

      read(a%lun_rstar) ((auxarray(idofn,ipoin),idofn=1,ndofn),ipoin=1,oldNpoin)

      if (a%kfl_inter/=1_ip) then
         array = auxarray
      else
         call a%Int_Restart%Interpolate(ndofn,auxarray,array)
      endif      

      call a%Memor%dealloc(ndofn,oldNpoin,auxarray,'auxarray','php_GeneralRestart')

   case (2)

      call a%Mesh%GetNpoin(npoin)
      call a%Mesh%GetNdime(ndime)
      ndofn = size(array,1)

      ! Write the previous time scheme used: one integer number
      !write(a%lun_rsta2,*) kfl_twost, kfl_trasg_pre, kfl_repro_pre 

      write(a%lun_rstar) ((array(idofn,ipoin),idofn=1,ndofn),ipoin=1,npoin)

   end select
   
   call a%Timer%Total%Toc
   call a%Timer%Restar%Toc

100 format(5x,'ERROR:   ',a)
101 format(5x,'WARNING: ',a)

end subroutine GeneralRestart
