subroutine lev_endite(a,itask)
   use typre
   use Mod_LevelSet
   use def_parame
   use Mod_int2str
   implicit none
   class(LevelSetProblem) :: a
   integer(ip) :: itask
   
   integer(ip) :: ndime
   
   character(50) :: names
   logical :: isALE
   
   integer(ip) :: ipoin,npoin
   
   interface
      subroutine lev_EnditeElmope(b)
         use typre
         use Mod_LevelSet
         implicit none
         class(LevelSetProblem), target :: b
      end subroutine
      
   end interface    
   
   select case(itask)

   case(1)

      !Update unknowns
      call a%Mesh%GetNdime(ndime)
      call a%Mesh%GetALE(isALE)
      !If it is ALE we do not compute the levelset advection, so we do not need to copy from unkno
      if ((isALE .eqv. .false.) .or. a%kfl_ForceEulerianAdvection == 1) then
         a%level(:,1) = a%unkno(1,:)  
      endif

      call lev_EnditeElmope(a)
      
      !Do not allow the selected nodes to become fluid 1 nodes
      call a%Mesh%GetNpoin(npoin)
      do ipoin = 1,npoin
         if (  (a%kfl_fixno(1,ipoin) == 2) .and. (a%level(ipoin,1) > 0.0_rp) ) then
            !write(*,*) 'correcting it step:' , a%istep,'ipoin: ', ipoin, 'values: ', a%level(ipoin,1), a%level(ipoin,3)
            !pause
            a%level(ipoin,1) = a%level(ipoin,3)
         endif
      enddo
      
!       names = 'level'//trim(adjustl(int2str(a%itera)))
!       call a%FilePostpr%postpr(a%level(:,1),names,a%istep,a%ctime,a%Mesh)
!       
!       names = 'poinstatus'//trim(adjustl(int2str(a%itera)))
!       call a%FilePostpr%postpr(a%CutMesh%Poinstatus,names,a%istep,a%ctime,a%Mesh)
      
      
   case(2)

      !Update u(n,i-1,*) <-- u(n,i,*)
      a%level(:,2) = a%level(:,1)
   
   end select  
   
end subroutine