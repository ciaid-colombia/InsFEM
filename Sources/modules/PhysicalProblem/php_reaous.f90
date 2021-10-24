subroutine php_reaous(a)
   use typre
   use Mod_PhysicalProblem
   implicit none

   class(PhysicalProblem) :: a
   
   integer(ip)   :: ipara,ndime,idime,ipoin,aux1
   character(4)  :: DUMM=CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)
   
   !For Listener
   real(rp), pointer     :: param(:) => NULL()
   character(5), pointer :: words(:) => NULL()
   integer(ip), pointer  :: nnpar => NULL(), nnwor => NULL()
   
   character(150) :: outstr
   
   call a%Listener%SetLunits(a%lun_pdata,a%lun_outpu)
   call a%Listener%getarrs(words,param,nnpar,nnwor)
   
   outstr = adjustl(trim(a%exmod))//'_REAOUS'
   call a%Mesh%GetNdime(ndime)

   ! Initializations.
   a%npp_inits = 0
   a%npp_stepi = 0
   a%pos_tinit = 0.0_rp
   a%pos_times = 0.0_rp
   a%nptra = 0                 ! Number of points to be tracked
   a%kfl_BoundaryConditionsReadStrategy = 0
   
   call a%SpecificReaous(0)
   
   ! Reach the section.
   call a%Listener%listen(outstr)
   do while(words(1)/='OUTPU')
      call a%Listener%listen(outstr)
   end do

   ! Begin to read data.
   do while(words(1)/='ENDOU')
      call a%Listener%listen(outstr)

      !Read strategy for boundary conditions 
      if(words(1)=='READS') then
         if(a%Listener%exists('DEFAU')) then
            a%kfl_BoundaryConditionsReadStrategy = 0
         elseif (a%Listener%exists('MANUF')) then
            a%kfl_BoundaryConditionsReadStrategy = 1
            call a%Listener%listen(outstr)
            a%ManufacturedBoundaryCondition = a%Listener%getint('MANUF',1,'BoundaryConditionType')
         end if

      !Should I print initial BC? by default yes
      elseif(words(1)=='PNTBC') then
         if(a%Listener%exists('NO   ')) then
             a%kfl_printBC = .false.
         endif
      
      !Starting time and step of post-process
      elseif(words(1)=='START') then
         if(a%Listener%exists('STEP ')) then
            a%npp_inits = a%Listener%getint('STEP ',0,'#Initial step to start postprocess')  
         end if
         if(a%Listener%exists('ATTIM')) then
            a%pos_tinit = a%Listener%getrea('ATTIM',0.0_rp,'#Initial step to start postprocess')
         end if

      ! Points to be tracked.
      else if(words(1)=='TRACK') then
         do while(words(1).ne.'ENDTR')
            call a%Listener%listen(outstr)
            !if(a%kfl_timei==1) then
               if(words(1)=='POINT') then
                  a%nptra=int(param(1))
                  call a%Mesh%GetNdime(ndime)
                  call a%Memor%alloc(ndime,a%nptra,a%cptra,'cptra','php_reaous')
                  do ipoin = 1,a%nptra
                     call a%Listener%listen(outstr)
                     ipara=1
                     do idime = 1,ndime
                        a%cptra(idime,ipoin) = (param(ipara))
                        ipara=ipara+1
                     end do
                  end do
               end if
            !end if
         end do
      endif
      
      call a%SpecificReaous(1)
      
   enddo
   
   call a%SpecificReaous(100)  !Final operations
   
end subroutine
   
   
