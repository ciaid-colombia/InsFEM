module Mod_lmn_Linearization
   use typre
   use mod_lmn_BaseElmope
   implicit none
   private
   public SetPointersLinearization

   integer(ip), allocatable :: kfl_IsSet
   
contains
   
   !Set Pointers
   subroutine SetPointersLinearization(itask)
      implicit none
      integer(ip) :: itask
      
      select case (itask)   
      
      case(0)
         allocate(kfl_IsSet)
         call a%Memor%allocObj(0,'kfl_IsSet','InitProcedurePointer',1)
         kfl_IsSet = -1
      
      case(1)
         if (kfl_IsSet == -1) then
            kfl_IsSet = 1
        
            if (a%kfl_linea == 2 .and. a%itera > a%npica .and. a%kfl_linop == 2) then
               call ConcatenateProcedures(ProcHook%InGaussElmats,Jacobian)
               call ConcatenateProcedures(ProcHook%InGaussElmats,Residuals)
               call ConcatenateProcedures(ProcHook%Assembly,AssemblyJacElmats)
               call ConcatenateProcedures(ProcHook%Assembly,ElmatsAssembly)
               call ConcatenateProcedures(ProcHook%ElmatsToZero,ToZeroJacArrays)
               call ConcatenateProcedures(ProcHook%AllocateArrays,AllocateJacArrays)
               call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateJacArrays)
               a%kfl_linop = 2
            else
               call ConcatenateProcedures(ProcHook%Assembly,ElmatsAssembly)
               a%kfl_linop = 1
            end if

         endif 
         
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select   
   end subroutine
   
   subroutine Jacobian
      implicit none
   
      call lmn_elmbuv_NR(e,dvol,acden,grvel,testf_mom,elmuvJ)
      call lmn_elmbtv_NR(e,acden,actex,dvol,grtem(1,:),ticon,elmtvJ)
      call lmn_elmbuw_NR(e,acden,dvol,grtem(1,:),testf_ene,elmuwJ)
      call lmn_elmbuq_NR(e,timom,dvol,acden,grvel,elmuqJ)

   end subroutine
   
   subroutine Residuals
      implicit none

      call lmn_elmbuv_NR_rhs(e,dvol,acden,grvel,gpvel,testf_mom,elrhu)
      call lmn_elmbtv_NR_rhs(e,acden,actex,dvol,grtem(1,:),gpvel,ticon,elrhu)
      call lmn_elmbuw_NR_rhs(e,acden,dvol,grtem(1,:),gpvel,testf_ene,elrht)
      call lmn_elmbuq_NR_rhs(e,timom,dvol,acden,grvel,gpvel,elrhp)

   end subroutine

   subroutine AssemblyJacElmats
      implicit none

      !Assembly J to elmat
      elmat(1:e%ndime,1:e%pnode,1:e%ndime,1:e%pnode) = elmat(1:e%ndime,1:e%pnode,1:e%ndime,1:e%pnode) + elmuvJ(1:e%ndime,1:e%pnode,1:e%ndime,1:e%pnode)
      elmat(e%ndime+1,1:e%pnode,1:e%ndime,1:e%pnode) = elmat(e%ndime+1,1:e%pnode,1:e%ndime,1:e%pnode) + elmuwJ(1,1:e%pnode,1:e%ndime,1:e%pnode)
      elmat(1:e%ndime,1:e%pnode,e%ndime+1,1:e%pnode) = elmat(1:e%ndime,1:e%pnode,e%ndime+1,1:e%pnode) + elmtvJ(1:e%ndime,1:e%pnode,1,1:e%pnode)
      elmat(e%ndime+2,1:e%pnode,1:e%ndime,1:e%pnode) = elmat(e%ndime+2,1:e%pnode,1:e%ndime,1:e%pnode) + elmuqJ(1,1:e%pnode,1:e%ndime,1:e%pnode)

   end subroutine

   subroutine ElmatsAssembly
      implicit none
   
      call ProcHook%PreDirichlet

      !Dirichlet Boundary Conditions
      call lmn_elmdir(a,e,elmat,elrhs)
      call a%LinearSystem%Assembly(e,elmat,elrhs)
   end subroutine

   subroutine AllocateJacArrays
      implicit none
   
      call a%Memor%alloc(e%ndime,e%mnode,e%ndime,e%mnode,elmuvJ,'elmuvJ','lmn_elmope')
      call a%Memor%alloc(1,e%mnode,e%ndime,e%mnode,elmuwJ,'elmuwJ','lmn_elmope')
      call a%Memor%alloc(1,e%mnode,e%ndime,e%mnode,elmuqJ,'elmuqJ','lmn_elmope')
      call a%Memor%alloc(e%ndime,e%mnode,1,e%mnode,elmtvJ,'elmtvJ','lmn_elmope')

   end subroutine

   subroutine DeallocateJacArrays
      implicit none
   
      call a%Memor%dealloc(e%ndime,e%mnode,e%ndime,e%mnode,elmuvJ,'elmuvJ','lmn_elmope')
      call a%Memor%dealloc(1,e%mnode,e%ndime,e%mnode,elmuwJ,'elmuwJ','lmn_elmope')
      call a%Memor%dealloc(1,e%mnode,e%ndime,e%mnode,elmuqJ,'elmuqJ','lmn_elmope')
      call a%Memor%dealloc(e%ndime,e%mnode,1,e%mnode,elmtvJ,'elmtvJ','lmn_elmope')

   end subroutine

   subroutine ToZeroJacArrays

      elmuvJ=0.0_rp
      elmuwJ=0.0_rp
      elmtvJ=0.0_rp
      elmuqJ=0.0_rp
   end subroutine

end module
