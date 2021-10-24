subroutine php_updbcs(a)
!-----------------------------------------------------------------------
!    Updates the a%velocity boundary conditions:
!    1. Before a time step begins
!    2. Before a global iteration begins
!    3. Before an inner iteration begins
!-----------------------------------------------------------------------
   use typre
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a
   real(rp), external      :: funcre
   real(rp), pointer       :: exnor(:,:) => NULL()
   integer(ip)             :: ibopo,idofbc,ipoin,npoin
   
   !write(*,*) 'php_updbcs: a%bctime',a%bctime
   
   !Before a time step
   if(a%kfl_conbc==0) then
      if(a%kfl_exacs/=0) then
         !In the case of an exact solution overwrite bc's using exact values
         !in your SpecificUpdbcs
      else
         call a%Mesh%GetNpoin(npoin)
         if (a%kfl_tsche_1st_current=='TR   ') then
            !0.5 t^n+1 + 0.5 t^n
            !write(*,*) 'php_updbcs: a%bctime,a%dtime',a%bctime,a%dtime
            do ipoin = 1,npoin
               call a%Mesh%GetBoundaryPoint(ipoin,ibopo,exnor)
!               if (ibopo==0) cycle ! the ipoin considered is not a boundary point
               do idofbc=1,a%ndofbc
                  if(a%kfl_fixno(idofbc,ipoin)==1) then
                     if(a%kfl_funno(ipoin)>0) then
                        if (a%kfl_funty(1,a%kfl_funno(ipoin)) > 0) then
                           a%bvess(idofbc,ipoin,1)=0.5*a%bvess(idofbc,ipoin,2)&
                              *( &
                                 funcre(a%funpa(a%kfl_funno(ipoin))%a,&
                                 a%kfl_funty(a%kfl_funno(ipoin),2),&
                                 a%kfl_funty(a%kfl_funno(ipoin),1),a%bctime-a%dtime) &
                              +funcre(a%funpa(a%kfl_funno(ipoin))%a,&
                                 a%kfl_funty(a%kfl_funno(ipoin),2),&
                                 a%kfl_funty(a%kfl_funno(ipoin),1),a%bctime) &
                              )
                        endif
                     end if
                  end if
               end do
               !write(*,*) 'php_updbcs: ipoin,a%bvess(:,ipoin,:)',ipoin,a%bvess(:,ipoin,:)
            end do
         else
            !prescribed value at bctime
            do ipoin = 1,npoin
               call a%Mesh%GetBoundaryPoint(ipoin,ibopo,exnor)
!               if (ibopo==0) cycle ! the ipoin considered is not a boundary point
               do idofbc=1,a%ndofbc
                  if(a%kfl_fixno(idofbc,ipoin)==1) then
                     if(a%kfl_funno(ipoin)>0) then
                        if (a%kfl_funty(1,a%kfl_funno(ipoin)) > 0) then
                           a%bvess(idofbc,ipoin,1)=a%bvess(idofbc,ipoin,2)&
                                 *funcre(a%funpa(a%kfl_funno(ipoin))%a,&
                                 a%kfl_funty(a%kfl_funno(ipoin),2),&
                                 a%kfl_funty(a%kfl_funno(ipoin),1),a%ctime)
                           !a%bvess(idofbc,ipoin,1)=a%bvess(idofbc,ipoin,2)&
                           !      *funcre(a%funpa(a%kfl_funno(ipoin))%a,&
                           !      a%kfl_funty(a%kfl_funno(ipoin),2),&
                           !      a%kfl_funty(a%kfl_funno(ipoin),1),a%bctime)
                        endif
                     end if
                  end if
               end do
               !write(*,*) 'php_updbcs: ipoin,a%bvess(:,ipoin,:)',ipoin,a%bvess(:,ipoin,:)
            end do
         end if
      end if
      call a%SpecificUpdbcs
   end if
   
end subroutine php_updbcs
