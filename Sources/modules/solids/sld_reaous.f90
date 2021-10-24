subroutine sld_reaous(a,itask)
   use typre
   use Mod_iofile
   use Mod_Solids
   implicit none

   class(SolidsProblem)  :: a
   integer(ip)           :: itask
   real(rp), pointer     :: param(:) => NULL()
   character(5), pointer :: words(:) => NULL()
   integer(ip), pointer  :: nnpar => NULL(),nnwor => NULL()
   integer(ip)           :: ndime,ipara
   
   
   !Initializations
   
   call a%Listener%getarrs(words,param,nnpar,nnwor)

   if (itask == 0) then  

   elseif (itask == 1) then  

      ! Step or time for post-process
      if(words(1)=='POSTP') then
          if(words(2)=='DISPL') then  
              if(words(3)=='STEPS') then
                  a%npp_stepi(1) = a%Listener%getint('STEPS',1,'#Postprocess step interval for displacements')
                  if(a%Listener%exists('ATTIM')) a%pos_times(1:10,1)=param(4:13)
              end if
          elseif(words(2)=='SIGMA') then  
              a%kfl_printSigma      = .true.
              a%kfl_printNodalSigma = .true.
              if(words(3)=='STEPS') then
                  a%npp_stepi(11) = a%Listener%getint('STEPS',1,'#Postprocess step interval for stresses')
                  if(a%Listener%exists('ATTIM')) a%pos_times(1:10,11)=param(4:13)
              end if
          elseif(words(2)=='GSIGM') then  
              a%kfl_printSigma      = .true.
              a%kfl_printGaussSigma = .true.
              if(words(3)=='STEPS') then
                  a%npp_stepi(20) = a%Listener%getint('STEPS',1,'#Postprocess step interval for stresses')
                  if(a%Listener%exists('ATTIM')) a%pos_times(1:10,20)=param(4:13)
              end if
          elseif(words(2)=='PRESS') then  
              a%kfl_printPress      = .true.
              a%kfl_printNodalPress = .true.
              if(words(3)=='STEPS') then
                  a%npp_stepi(12) = a%Listener%getint('STEPS',1,'#Postprocess step interval for pressure')
                  if(a%Listener%exists('ATTIM')) a%pos_times(1:10,12)=param(4:13)
              end if
          elseif(words(2)=='GPRES') then  
              a%kfl_printPress      = .true.
              a%kfl_printGaussPress = .true.
              if(words(3)=='STEPS') then
                  a%npp_stepi(23) = a%Listener%getint('STEPS',1,'#Postprocess step interval for pressure')
                  if(a%Listener%exists('ATTIM')) a%pos_times(1:10,23)=param(4:13)
              end if
          elseif(words(2)=='STRES') then  
              a%kfl_printStress = .true.
              a%kfl_GaussStress = .true.
              if(words(3)=='STEPS') then
                  a%npp_stepi(2) = a%Listener%getint('STEPS',1,'#Postprocess step interval for stresses')
                  if(a%Listener%exists('ATTIM')) a%pos_times(1:10,2)=param(4:13)
              end if
          elseif(words(2)=='STRAI') then  
              a%kfl_printStrain = .true.
              a%kfl_GaussStrain = .true.
              if(words(3)=='STEPS') then
                  a%npp_stepi(3) = a%Listener%getint('STEPS',1,'#Postprocess step interval for strains')
                  if(a%Listener%exists('ATTIM')) a%pos_times(1:10,3)=param(4:13)
              end if
              if (words(4) == 'MODEL') then
                  a%sld_strain = words(5)
              end if
          elseif(words(2)=='INTTR') then  
              if(words(3)=='STEPS') then
                  a%npp_stepi(4) = a%Listener%getint('STEPS',1,'#Postprocess step interval for internal tractions')
                  if(a%Listener%exists('ATTIM')) a%pos_times(1:10,4)=param(4:13)
              end if
          elseif(words(2)=='FLDTR') then  
              if(words(3)=='STEPS') then
                  a%npp_stepi(5) = a%Listener%getint('STEPS',1,'#Postprocess step interval for fluid traction - FSI')
                  if(a%Listener%exists('ATTIM')) a%pos_times(1:10,5)=param(4:13)
              end if
          elseif(words(2)=='EXTTR') then  
              if(words(3)=='STEPS') then
                  a%npp_stepi(6) = a%Listener%getint('STEPS',1,'#Postprocess step interval for external tractions')
                  if(a%Listener%exists('ATTIM')) a%pos_times(1:10,6)=param(4:13)
              end if
          elseif(words(2)=='VELOC') then  
              if(words(3)=='STEPS') then
                  a%npp_stepi(7) = a%Listener%getint('STEPS',1,'#Postprocess step interval for velocity')
                  if(a%Listener%exists('ATTIM')) a%pos_times(1:10,7)=param(4:13)
              end if
          elseif(words(2)=='PSTRA') then  
              if(words(3)=='STEPS') then
                  a%npp_stepi(8) = a%Listener%getint('STEPS',1,'#Postprocess step interval for ppal strains')
                  if(a%Listener%exists('ATTIM')) a%pos_times(1:10,8)=param(4:13)
              end if
          elseif(words(2)=='ACCEL') then  
              if(words(3)=='STEPS') then
                  a%npp_stepi(9) = a%Listener%getint('STEPS',1,'#Postprocess step interval for velocity')
                  if(a%Listener%exists('ATTIM')) a%pos_times(1:10,9)=param(4:13)
              end if
          elseif(words(2)=='NSTRE') then  
              a%kfl_NodalStress = .true.
              a%kfl_printStress = .true.
              if(words(3)=='STEPS') then
                  a%npp_stepi(21) = a%Listener%getint('STEPS',1,'#Postprocess step interval for stresses')
                  if(a%Listener%exists('ATTIM')) a%pos_times(1:10,21)=param(4:13)
              end if
          elseif(words(2)=='NSTRA') then  
              a%kfl_NodalStrain = .true.
              a%kfl_printStrain = .true.
              if(words(3)=='STEPS') then
                  a%npp_stepi(22) = a%Listener%getint('STEPS',1,'#Postprocess step interval for strains')
                  if(a%Listener%exists('ATTIM')) a%pos_times(1:10,22)=param(4:13)
              end if
              if (words(4) == 'MODEL') then
                  a%sld_strain = words(5)
              end if
          end if

      end if

   elseif (itask == 100) then   

      if (a%nptra/=0) then
         if (a%MPIrank == a%MPIroot) then
            a%fil_trap = adjustl(trim(a%namda))//'.sld.tp'
            call iofile(zero,a%lun_trap,a%fil_trap,'SLD TRACKING OF POINTS')
         endif
      endif 

   endif

   call a%SolidSpecificReaous(itask)

   1 format(a)
   5 format(a,1x,i2)

   101 format(5x,'WARNING: ',a)

end subroutine sld_reaous

