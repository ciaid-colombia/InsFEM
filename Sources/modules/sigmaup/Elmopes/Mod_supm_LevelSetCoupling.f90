module Mod_supm_LevelSetCoupling
   use typre
   use Mod_supm_BaseElmope
   use Mod_supm_StokesProblem
   use Mod_supm_EnrichElement
   implicit none
   private
   public SetPointersLevelSetCoupling
   integer(ip), allocatable :: kfl_IsSet
   
contains
   
   subroutine SetPointersLevelSetCoupling(itask,task)
      integer(ip) :: itask
      character(6) :: task
      !procedure() :: NULL()
      
      select case (itask)   
      case(0)
         allocate(kfl_IsSet)
         call a%Memor%allocObj(0,'kfl_IsSet','InitProcedurePointer',1)
         kfl_IsSet = -1
      
      case(1)
         if (kfl_IsSet == -1) then
            kfl_IsSet = 1
            if(a%kfl_colev==1)then
               call ConcatenateProcedures(ProcHook_Initializations,AllocLevelsetTF) 
               call ConcatenateProcedures(ProcHook_Finalizations,DeallocLevelSetTF)
               call PrependProcedure(ProcHook_PreGauss,CutelementsTF)
               call PrependProcedure(ProcHook_PhysicalProp,FluidPropertiesTF) 
               
               if (task .eq. 'Elmope') then
                  if(a%kfl_fsurf==1)then
                     call SetPointersStokesProblem(1)
                     call ConcatenateProcedures( ProcPointer%PreAssembly_sup,FreeSurfMatsTF_Elmope)  
                  end if
                  call SetPointersEnrichElement(1)
               else if (task .eq. 'Endite') then
                  if(a%kfl_fsurf==1)then
                     call ConcatenateProcedures( ProcPointer%PreAssembly_sup,FreeSurfMatsTF_Endite)
                  end if 
               end if
               
            end if
         end if
         
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      end select   
      
   end subroutine   
   

   subroutine AllocLevelsetTF
      implicit none   
      integer(ip) :: auxdim,auxdim1
      auxdim=(e%ndime+1)
      auxdim1=(e%ndime-1)*(e%ndime-1)+2
      call a%Memor%alloc(auxdim*auxdim1,weigp,'weigp','nsm_elmope')
      call a%Memor%alloc(e%ndime,auxdim*auxdim1,xloc,'xloc','nsm_elmope') 
   end subroutine
   
   subroutine DeallocLevelSetTF
      implicit none      
      integer(ip) :: auxdim,auxdim1
      auxdim=(e%ndime+1)
      auxdim1=(e%ndime-1)*(e%ndime-1)+2    
      call a%Memor%dealloc(auxdim*auxdim1,weigp,'weigp','nsm_elmope')
      call a%Memor%dealloc(e%ndime,auxdim*auxdim1,xloc,'xloc','nsm_elmope') 
   end subroutine
   
   subroutine CutelementsTF
      implicit none
      integer(ip)  :: elemStatus       
      
      call a%CutMesh%GetElementType(ielem,elemStatus)
      if(elemStatus ==0)then     
         ngaus_total=0
         call a%CutMesh%GetNgaussSide(ielem,ngauss_plus,ngauss_minus)
         call a%CutMesh%GetWeigpCut(ielem,e%ndime,weigp)
         call a%CutMesh%GetXlocCut(ielem,e%ndime,xloc)
         ngaus_total = ngauss_minus+ngauss_plus
         !The rutine give the needed shape functions associated 
         if(a%kfl_fsurf==1) weigp(1:ngauss_minus)=0.0_rp         
         call e%SetParticularGaussPoints(a%Memor,ngaus_total,xloc,weigp(:))
      end if
   
   end subroutine
   
   subroutine FluidPropertiesTF
      implicit none
      integer(ip)  :: elemStatus      
      call a%CutMesh%GetElementType(ielem,elemStatus)
      if(elemStatus==1)then
         imat=1
         acden=a%MatProp(imat)%densi
         acvis=a%MatProp(imat)%visco           
      elseif(elemStatus==-1)then
         imat=2
         acden=a%MatProp(imat)%densi
         acvis=a%MatProp(imat)%visco           
      elseif(elemStatus ==0)then
         ngauss_minus=0
         ngauss_plus=0
         call a%CutMesh%GetNgaussSide(ielem,ngauss_plus,ngauss_minus)
         if(e%igaus<=ngauss_minus)then
            imat=2
            acden=a%MatProp(imat)%densi
            acvis=a%MatProp(imat)%visco               
         elseif(e%igaus>ngauss_minus)then
            imat=1
            acden=a%MatProp(imat)%densi
            acvis=a%MatProp(imat)%visco               
         end if
      end if
   end subroutine
   
   subroutine FreeSurfMatsTF_Elmope
      implicit none 
      integer(ip)  :: elemStatus,poinStatus       
      call a%CutMesh%GetElementType(ielem,elemStatus)
      if(elemStatus==-1)then
         do inode=1,e%pnode
            ipoin=e%lnods(inode)
            call a%CutMesh%GetPointType(ipoin,poinStatus)            
            if(poinStatus==-1)then
               elmat(:,inode,:,1:e%pnode) = 0.0_rp
               elrhs(:,inode)=0.0_rp                   
            end if         
         end do        
      end if
   end subroutine  
   
   subroutine FreeSurfMatsTF_Endite
      implicit none 
      integer(ip)  :: elemStatus       
      integer(ip) :: kfl_nonlinear
      
      call a%CutMesh%GetElementType(ielem,elemStatus)
      !Non-linear elements
      call a%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)      
      if(elemStatus==-1)then
         !OSS 
         if (a%kfl_repro >= 1) then
            ! elrep(:,:) = 0.0_rp
              elres(:,:) = 0.0_rp
         end if
         !SplitOSS
         if (a%MatProp(imat)%lawvi<0) then
            if(a%kfl_repro == 2 .or. a%kfl_repro == 3)then         
               elresdivs(:,:) = 0.0_rp
               elresconv(:,:) = 0.0_rp
               elresgrap(:,:) = 0.0_rp 
               
               if(kfl_nonlinear==1)then
                  elreslapl(:,:) = 0.0_rp
               end if
            end if
         end if 
         !Shock capturing
         if(a%MatProp(imat)%lawvi<0 .and. a%kfl_shock == 4) then  
            elrepGrad(:,:) = 0.0_rp            
         end if         
         !DEVSS
         if (a%MatProp(imat)%lawvi<0 .and. a%kfl_shock == 2) then
            elrepSGrad(:,:) = 0.0_rp 
         end if         
      end if
   end subroutine  

end module
