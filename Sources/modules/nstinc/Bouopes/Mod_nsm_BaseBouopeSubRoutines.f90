module Mod_nsm_BaseBouopeRoutines

   use Mod_nsm_BaseElmope
   use Mod_nsm_BaseBouope

   implicit none

   contains

   subroutine nsm_bouLoop
     integer(ip)  :: inodb,icount,npoinLocal
     integer(ip), pointer :: lbody => NULL()

     !Initializations
     call a%Mesh%GetNboun(nboun)   
     call a%Mesh%GetNpoinLocal(npoinLocal)     
     
     call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','nsi_forces')      
     call ProcHook_AllocateArrays
     call ProcHook_Initializations
     
     boundaries: do iboun=1,nboun
        !Load Element
        call a%Mesh%BoundaryLoad(iboun,e)
        call a%Mesh%GetLbody(lbody,iboun)
        ibody=lbody
        
        icount = 0
        do inodb = 1,e%pnodb
            if (e%lnodb(inodb) <= npoinLocal) then
                icount = icount +1
            endif
        enddo
        weightfactor = real(icount)/real(e%pnodb)
 
        call ProcHook_ElmatsToZero

        call ProcHook_Gathers

        call e%elmdel
        call e%elmlen
        dsurf0 = 0.0_rp           
        call ProcHook_PreGauss

        do igaub=1,e%pgaub
           e%igaub = igaub

           !Derivatives at the boundary
           call e%elmderb
           !Calculate exterior Normal
           call e%bounor
           
           dsurf=e%weigb(e%igaub)*e%eucta
           dsurf0 = dsurf0 + dsurf

           call ProcHook_Interpolates

           !Physical Parameters
           call a%GetPhysicalParameters(imat,acden,acvis)         
           call ProcHook_PhysicalProp

           call ProcHook_InGauss
           call ProcHook_InGaussElmatsAssembly

        end do 
        call ProcHook_AssemblyEndite
     end do boundaries 

     call ProcHook_PostLoop
     call ProcHook_Finalizations
     
     call ProcHook_DeallocateArrays
     call a%Mesh%ElementDealloc(e,a%Memor,'DefaultRule','nsi_froces')

  end subroutine nsm_bouLoop

end module Mod_nsm_BaseBouopeRoutines
