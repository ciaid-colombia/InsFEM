subroutine lev_reaous(a,itask)
   use typre
   use Mod_Listen
   use Mod_LevelSet
   implicit none

   class(LevelSetProblem) :: a
   integer(ip) :: itask
   
   !For Listener
   real(rp), pointer     :: param(:) => NULL()
   character(5), pointer :: words(:) => NULL()
   integer(ip), pointer  :: nnpar => NULL(),nnwor => NULL()
   
   integer(ip) :: igauge,ndime
   
   call a%Listener%getarrs(words,param,nnpar,nnwor)
 
   !Initializations
   if (itask == 0) then
      
      
   !Output section   
   elseif(itask == 1) then
     
      !Step or time for post-process
      if(words(1)=='POSTP') then
         if(words(2)=='LEVEL') then  
            if(words(3)=='STEPS') then
               a%npp_stepi(1) = a%Listener%getint('STEPS',1,'#Postprocess step interval for L')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,1)=param(4:13)
            end if

         end if     

      elseif (words(1)=='HEIGH') then
         a%nHeightGauges = a%Listener%param(1)
         if (a%nHeightGauges > 0) then
            allocate(a%HeightGauges(a%nHeightGauges))
            call a%Memor%allocObj(0,'HeightGauges','lev_reaous',a%nHeightGauges)
            
            call a%Mesh%GetNdime(ndime)
            
            do igauge = 1,a%nHeightGauges
               call a%Listener%Listen('lev_reaous')
               a%HeightGauges(igauge)%Origin = a%Listener%param(1:3)
               a%HeightGauges(igauge)%DirectionVector = a%Listener%param(4:6)
               
               call a%HeightGauges(igauge)%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%Mpirank)
               call a%HeightGauges(igauge)%SetNdime(ndime)
               
               
            enddo
         endif
         
      end if

  !Finalizations   
  elseif (itask == 100) then 
  
  
  
  endif

end subroutine