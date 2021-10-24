subroutine tem_reabcs(a,itask,kflag)
   use Mod_Memor
   use typre
  
   use Mod_Temperature
   implicit none
   class(TemperatureProblem) :: a
   integer(ip) :: itask
   integer(ip), optional :: kflag

   integer(ip) :: npoin
   
   !Initializations
   if (itask == 0) then
      if(a%kfl_sourc==2) then
         call a%Mesh%GetNpoin(npoin)
         call a%Memor%alloc(npoin,a%PointwiseSource,'PointwiseSource','tem_reabcs')
      end if
      
    
   !Header
   elseif(itask == 1) then
      if(a%Listener%exists('WALLD')) then
         a%delta=a%Listener%getrea('WALLD',0.0_rp,'#Distance to the wall')
      end if
      
      
   !Finalization
   elseif(itask == 100) then
      
   
   endif
   
end subroutine


subroutine tem_ReadOnNodes(a)
   use typre
   use Mod_Mesh
   use Mod_MPIObject
   use Mod_Memor
   use Mod_Listen
   use Mod_PhysicalProblem
   use Mod_Temperature
   implicit none
   class(TemperatureProblem) :: a
  
   !On Nodes
   if(a%kfl_conbc==0) then
      if(a%kfl_sourc==2) then
         a%PointwiseSource(a%gipoin) = a%Listener%param(5)
      end if
      
   else
      if(a%kfl_sourc==2) a%PointwiseSource(a%gipoin) = a%Listener%param(4)
   end if
end subroutine



subroutine tem_ReadOnBoundaries(a)
   use typre
   use Mod_Mesh
   use Mod_MPIObject
   use Mod_Memor
   use Mod_Listen
   use Mod_PhysicalProblem
   use Mod_Temperature
   implicit none
   class(TemperatureProblem) :: a
   
   !Robin
   if(a%kfl_fixbo(a%giboun)==3) then
      
      if(a%kfl_conbc==0) a%kfl_funbo(a%giboun) = int(a%Listener%param(a%gipsta+4))
      if(a%kfl_funbo(a%giboun)/=0) then
         allocate(a%bvnat(a%giboun)%a(6))
         
         a%bvnat(a%giboun)%a(1)=a%Listener%param(a%gipsta+1)
         a%bvnat(a%giboun)%a(2)=a%Listener%param(a%gipsta+2)
         a%bvnat(a%giboun)%a(3)=a%Listener%param(a%gipsta+3)
         a%bvnat(a%giboun)%a(4)=a%Listener%param(a%gipsta+1)
         a%bvnat(a%giboun)%a(5)=a%Listener%param(a%gipsta+2)
         a%bvnat(a%giboun)%a(6)=a%Listener%param(a%gipsta+3)
      else
         allocate(a%bvnat(a%giboun)%a(3))
         a%bvnat(a%giboun)%a(1)=a%Listener%param(a%gipsta+1)
         a%bvnat(a%giboun)%a(2)=a%Listener%param(a%gipsta+2)
         a%bvnat(a%giboun)%a(3)=a%Listener%param(a%gipsta+3)
      end if
   !Wall law
   else if(a%kfl_fixbo(a%giboun)==4) then     
      
      if(a%kfl_conbc==0) a%kfl_funbo(a%giboun) = int(a%Listener%param(a%gipsta+2))
      if(a%kfl_funbo(a%giboun)/=0) then
         allocate(a%bvnat(a%giboun)%a(2))
         a%bvnat(a%giboun)%a(1)=a%Listener%param(a%gipsta+1)
         a%bvnat(a%giboun)%a(2)=a%Listener%param(a%gipsta+1)
      else
         allocate(a%bvnat(a%giboun)%a(1))
         a%bvnat(a%giboun)%a(1)=a%Listener%param(a%gipsta+1)              
      end if
   !Dust Transport
   elseif (a%kfl_fixbo(a%giboun) == 5) then
      !Nothing to be done here
   
   !Do nothing boundary condition
   elseif (a%kfl_fixbo(a%giboun) == 6) then
      !Do nothing boundary condition
      !Nothing to be done here
      
   end if
   
end subroutine

subroutine tem_ReadOnElements(a)
   use typre
   use Mod_Listen
   use Mod_Temperature
   implicit none
   class(TemperatureProblem) :: a

   integer(ip) :: ielem,imaterial
   
   ielem = a%gielem
   imaterial = int(a%Listener%param(a%gipsta+1))
   if (imaterial <= 0 .or. imaterial > a%NumberOfMaterials) call runend('Wrong material number')
   a%ElementMaterials(ielem) = imaterial
end subroutine

   
