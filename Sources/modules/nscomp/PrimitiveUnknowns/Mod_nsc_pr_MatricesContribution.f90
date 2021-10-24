module Mod_nsc_pr_MatricesContribution
   use typre
   use Mod_nsc_pr_BaseElmope
   implicit none
   
   integer(ip), allocatable :: kfl_IsSet

contains

   subroutine SetPointersMatricesContribution(itask)
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
         
            call ConcatenateProcedures(ProcHook_nsc_pr_Initializations,AllocateMatricesContribution)
            call ConcatenateProcedures(ProcHook_nsc_pr_Finalizations,DeallocateMatricesContribution)

         endif
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine
   
   
   !-----------------------------------------------------------------------
   subroutine AllocateMatricesContribution

      call a%Memor%alloc(e%mnode,Edd,'Edd','nsc_pr_elmope')
      call a%Memor%alloc(e%ndime,e%mnode,Edm,'Edm','nsc_pr_elmope')
      call a%Memor%alloc(e%mnode,Ede,'Ede','nsc_pr_elmope')
      call a%Memor%alloc(e%ndime,e%mnode,Emd,'Emd','nsc_pr_elmope')
      call a%Memor%alloc(e%ndime,e%ndime,e%mnode,Emm,'Emm','nsc_pr_elmope')
      call a%Memor%alloc(e%ndime,e%mnode,Eme,'Eme','nsc_pr_elmope')
      call a%Memor%alloc(e%mnode,Eed,'Eed','nsc_pr_elmope')
      call a%Memor%alloc(e%ndime,e%mnode,Eem,'Eem','nsc_pr_elmope')
      call a%Memor%alloc(e%mnode,Eee,'Eee','nsc_pr_elmope')
      call a%Memor%alloc(e%ndime,Erhm,'Erhm','nsc_pr_elmope')

   end subroutine

   subroutine ComputeEulerContribution
      implicit none
      
      Edd = Atdd*LHSdtinv + Add
      Edm =                 Adm
      Ede = Atde*LHSdtinv + Ade
      Emd = Atmd*LHSdtinv + Amd
      Emm = Atmm*LHSdtinv + Amm
      Eme = Atme*LHSdtinv + Ame
      Eed = Ated*LHSdtinv + Aed
      Eem = Atem*LHSdtinv + Aem
      Eee = Atee*LHSdtinv + Aee

   end subroutine

   subroutine ComputeRHSContribution(e,eltemp,eltemv,eltemt,elexp,elexv,elext,Erhd,Erhm,Erhe)
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)       :: eltemp,elexp
      real(rp),    intent(in)       :: eltemv(e%ndime),elexv(e%ndime)
      real(rp),    intent(in)       :: eltemt,elext
      real(rp),    intent(inout)    :: Erhd,Erhe
      real(rp),    intent(inout)    :: Erhm(e%ndime)
  
      real(rp) :: auxme
 
      auxme = dot_product(gpadv,eltemv)
      
      Erhd    = gpden*acbeta*eltemp - gpden*acalpha*eltemt + elexp
      Erhm(:) = gpden*acbeta*gpadv(:)*eltemp + gpden*eltemv(:) -&
                gpden*acalpha*gpadv(:)*eltemt + elexv(:)
      Erhe    = (aux*acbeta-acalpha*gpadt)*eltemp +&
                 gpden*auxme+&
                (gpden*accph-aux*acalpha)*eltemt + elext

   end subroutine

   subroutine DeallocateMatricesContribution

      call a%Memor%dealloc(e%mnode,Edd,'Edd','nsc_pr_elmope')
      call a%Memor%dealloc(e%ndime,e%mnode,Edm,'Edm','nsc_pr_elmope')
      call a%Memor%dealloc(e%mnode,Ede,'Ede','nsc_pr_elmope')
      call a%Memor%dealloc(e%ndime,e%mnode,Emd,'Emd','nsc_pr_elmope')
      call a%Memor%dealloc(e%ndime,e%ndime,e%mnode,Emm,'Emm','nsc_pr_elmope')
      call a%Memor%dealloc(e%ndime,e%mnode,Eme,'Eme','nsc_pr_elmope')
      call a%Memor%dealloc(e%mnode,Eed,'Eed','nsc_pr_elmope')
      call a%Memor%dealloc(e%ndime,e%mnode,Eem,'Eem','nsc_pr_elmope')
      call a%Memor%dealloc(e%mnode,Eee,'Eee','nsc_pr_elmope')
      call a%Memor%dealloc(e%ndime,Erhm,'Erhm','nsc_pr_elmope')

   end subroutine

end module
