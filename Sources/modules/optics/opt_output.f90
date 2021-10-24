subroutine opt_output(a,itask)
   use typre
   use Mod_Optics
   use Mod_Postpr
   implicit none
   class(OpticsProblem) :: a
   integer(ip) :: itask
   integer(ip), save       :: dopost(25)
   integer(ip)             :: itime,istat,ifiel,ndime
   real(rp)                :: dummr
   
   select case(itask)

   !End of a time step.
   case(0)
      
   case(1)
      !Cn2
      if(a%npp_stepi(1)/=0) then
         call a%FilePostpr%postpr(a%avg_cn2(:,1),'CN2_MEAN',a%istep,a%ctime,a%Mesh)
         call a%FilePostpr%postpr(a%avg_cn2(:,2),'CN2_FROM_MEAN',a%istep,a%ctime,a%Mesh)
         if (a%kfl_Avg1DCn2 == 1) then
            call a%FilePostpr%postAvg1D(a%Avg1DIdime,a%avg_cn2(:,2),'Averaged_CN2_FROM_MEAN',a%istep,a%ctime,a%Mesh,a%Memor)
            call a%FilePostpr%postAvg1D(a%Avg1DIdime,a%avg_vdiss(:),'OPT_AVG_VDISS',a%istep,a%ctime,a%Mesh,a%Memor)
            call a%FilePostpr%postAvg1D(a%Avg1DIdime,a%avg_tdiss(:),'OPT_AVG_TDISS',a%istep,a%ctime,a%Mesh,a%Memor)
         endif
      end if
      
   end select

  
   !Decide which postprocesses need to be done
   call SetupDoPost(itask,a%istep,size(a%npp_stepi),a%npp_inits,a%npp_stepi,a%pos_alrea,dopost)
   
   !Do the actual postprocess
   !CN2
   if (dopost(1) == 1) then
       call a%FilePostpr%postpr(a%cn2(:),'CN2',a%istep,a%ctime,a%Mesh)
            if (a%kfl_Avg1DCn2 == 1) then
               call a%FilePostpr%postAvg1D(a%Avg1DIdime,a%cn2(:),'Averaged_CN2',a%istep,a%ctime,a%Mesh,a%Memor)
            endif
   endif
   
   !CT2
   if (dopost(2) == 1) then
      call a%FilePostpr%postpr(a%ct2,'CT2',a%istep,a%ctime,a%Mesh)
            if (a%kfl_Avg1DCn2 == 1) then
               call a%FilePostpr%postAvg1D(a%Avg1DIdime,a%ct2(:),'Averaged_CT2',a%istep,a%ctime,a%Mesh,a%Memor)
            endif
   endif

end subroutine
