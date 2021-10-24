subroutine msh_InitializeElementDataStructures(a)
    use Mod_Element
    use Mod_ElementDataStructures
    use Mod_Mesh
    implicit none
    class(FemMesh), target :: a
    class(FiniteElement), pointer :: e => NULL()
    type(ElementDataStructure), pointer :: iED => NULL()
    type(ElementDataStructure), pointer :: EDS(:) => NULL()

    type(BoundaryDataStructure), pointer :: iBD => NULL()
    type(BoundaryDataStructure), pointer :: BDS(:) => NULL()

    integer(ip) :: allocsize

    interface 
        subroutine ElementalInitialize(a,e,EDS,BDS)
            use typre
            use Mod_Element
            use Mod_ElementDataStructures
            use Mod_Mesh
            implicit none
            class(FemMesh), target :: a
            class(FiniteElement), pointer :: e
            type(ElementDataStructure), pointer :: EDS(:)
            type(BoundaryDataStructure), pointer :: BDS(:)
        end subroutine
    end interface

    !Temporarily disable Usedatastructures so that the element info is actually computed
    a%kfl_UseElementDataStructures = .false.
    allocate(a%ElementDataStructure(a%nelem))
    allocate(a%BoundaryDataStructure(a%nboun))
    call a%ElementAlloc(e,a%Memor,'DefaultRule','InitializeElementDataStructures')
    EDS => a%ElementDataStructure 
    BDS => a%BoundaryDataStructure
    call ElementalInitialize(a,e,EDS,BDS)
    call a%ElementDeAlloc(e,a%Memor,'DefaultRule','InitializeElementDataStructures')

    allocate(a%ElementDataStructureClosedRule(a%nelem))
    allocate(a%BoundaryDataStructureClosedRule(a%nboun))
    call a%ElementAlloc(e,a%Memor,'ForceClosedRule','InitializeElementDataStructures')
    EDS => a%ElementDataStructureClosedRule 
    BDS => a%BoundaryDataStructureClosedRule
    call ElementalInitialize(a,e,EDS,BDS)
    call a%ElementDeAlloc(e,a%Memor,'ForceClosedRule','InitializeElementDataStructures')

    !Reactivate use data structures
    a%kfl_UseElementDataStructures = .true.


    allocsize = a%nelem*(a%ndime*a%mnode*a%mgaus+a%ndime*a%mnode+a%ndime*a%ndime+a%ndime+a%ntens*a%mnode*a%mlapl*a%mgaus+a%mgaus)
    allocsize = allocsize + a%nboun*(a%ndime*a%mnode*a%mgaub+a%ndime*a%ndime*a%mgaub+a%mgaub)
    call a%Memor%allocObj(0_ip,'ElementDataStructures','msh_InitializeElementDataStructures',allocsize)

end subroutine


subroutine ElementalInitialize(a,e,EDS,BDS)
   use typre
   use Mod_Mesh
   use Mod_Element
   use Mod_ElementDataStructures
   implicit none

   class(FemMesh), target :: a
   class(FiniteElement), pointer :: e
   type(ElementDataStructure), pointer :: EDS(:)
   type(BoundaryDataStructure), pointer :: BDS(:)
   type(ElementDataStructure), pointer :: iED => NULL()
   type(BoundaryDataStructure), pointer :: iBD => NULL()
   integer(ip) :: ielem,igaus,igaub,iboun,iface

    do ielem = 1,a%nelem
        call a%ElementLoad(ielem,e)

        iED => EDS(ielem)

        !Allocate data structure
        allocate(iED%cartd(e%ndime,e%pnode,e%pgaus), &
            iED%cartg(e%ndime,e%pnode), &
            iED%tragl(e%ndime,e%ndime), &
            iED%hleng(e%ndime), &
            iED%hessi(e%ntens,e%pnode*e%mlapl,e%pgaus), &
            iED%elvec(e%ndime,e%ndime,e%pface,e%pgaub), &
            iED%elnor(e%pface,e%pgaus), &
            iED%detjm(e%pgaus))

        !Cartesian derivatives and Jacobian at center of gravity
        call e%elmdcg

        iED%cartg = e%cartd(:,1:e%pnode)
        iED%tragl = e%xjaci
        iED%detjmg = e%detjm

        !Element length at center of gravity
        call e%elmlen

        iED%hleng = e%hleng

        do iface = 1,e%pface
           e%iface = iface
           call a%FaceLoad(iface,ielem,e)
           do igaub = 1,e%pgaub
              e%igaub = igaub
              call e%elenor
              iED%elvec(:,:,e%iface,e%igaub) = e%baloc
              iED%elnor(e%iface,e%igaub) = e%eucta
           end do
        end do

        !Gauss Point Loop
        gauss_points : do igaus=1,e%pgaus
           e%igaus = igaus

           call e%elmder

           iED%cartd(:,:,e%igaus) = e%cartd(:,1:e%pnode)
           iED%detjm(e%igaus) = e%detjm

           if (e%linea == 1) then
           else
               call e%elmhes
               iED%hessi(:,:,e%igaus) = e%hessi(:,1:e%pnode*e%mlapl)
           endif

        enddo gauss_points
    enddo

    !Boundaries information
    do iboun = 1,a%nboun
       !Load Element
       call a%BoundaryLoad(iboun,e)
       ielem = e%lboel(e%pnodb+1) 

       iBD => BDS(iboun)

       allocate( iBD%cartb(e%ndime,e%pnode,e%pgaub), &
           iBD%baloc(e%ndime,e%ndime,e%pgaub)    ,&
           iBD%eucta(e%pgaub)    )

       !Linear derivatives just in case they are needed
       call e%elmdel

       !Boundary gauss point loop
       do igaub = 1,e%pgaub
          e%igaub = igaub

          !Derivatives at the boundary

          call e%elmderb

          iBD%cartb(:,:,e%igaub) = e%cartb(:,1:e%pnode)

          !Calculate exterior Normal
          call e%bounor

          iBD%baloc(:,:,e%igaub) = e%baloc
          iBD%eucta(e%igaub) = e%eucta
       enddo
    enddo

end subroutine

subroutine msh_FinalizeElementDataStructures(a)
   use Mod_Mesh
   use Mod_Element
   use Mod_ElementDataStructures
   implicit none
   class(FemMesh), target :: a
   class(FiniteElement), pointer :: e => NULL()
   type(ElementDataStructure), pointer :: EDS(:) => NULL()
   type(BoundaryDataStructure), pointer :: BDS(:) => NULL()
   integer(ip) :: allocsize

   interface 
      subroutine ElementalFinalize(a,EDS,BDS)
         use typre
         use Mod_Mesh
         use Mod_Element
         use Mod_ElementDataStructures
         implicit none
         class(FemMesh), target :: a
         type(ElementDataStructure), pointer :: EDS(:)
         type(BoundaryDataStructure), pointer :: BDS(:)
      end subroutine
   end interface

   EDS => a%ElementDataStructure 
   BDS => a%BoundaryDataStructure
   call ElementalFinalize(a,EDS,BDS)
   deallocate(a%ElementDataStructure)
   deallocate(a%BoundaryDataStructure)


   EDS => a%ElementDataStructureClosedRule 
   BDS => a%BoundaryDataStructureClosedRule
   call ElementalFinalize(a,EDS,BDS)
   deallocate(a%ElementDataStructureClosedRule)
   deallocate(a%BoundaryDataStructureClosedRule)

   allocsize = a%nelem*(a%ndime*a%mnode*a%mgaus+a%ndime*a%mnode+a%ndime*a%ndime+a%ndime+a%ntens*a%mnode*a%mlapl*a%mgaus+a%mgaus)
   allocsize = allocsize + a%nboun*(a%ndime*a%mnode*a%mgaub+a%ndime*a%ndime*a%mgaub+a%mgaub)
   call a%Memor%deallocObj(0_ip,'ElementDataStructures','msh_InitializeElementDataStructures',allocsize)

end subroutine

subroutine ElementalFinalize(a,EDS,BDS)
   use typre
   use Mod_Mesh
   use Mod_Element
   use Mod_ElementDataStructures
   implicit none

    class(FemMesh), target :: a
    type(ElementDataStructure), pointer :: iED => NULL()
    type(ElementDataStructure), pointer :: EDS(:)

    type(BoundaryDataStructure), pointer :: iBD => NULL()
    type(BoundaryDataStructure), pointer :: BDS(:)
    integer(ip) :: ielem,igaus,igaub,iboun

    do ielem = 1,a%nelem
        iED => EDS(ielem)

        !DeAllocate data structure
        deallocate(   iED%cartd, &
            iED%cartg, &
            iED%tragl, &
            iED%hleng, &
            iED%hessi, &
            iED%elvec, &
            iED%elnor, &
            iED%detjm   )
    enddo

    !Boundaries information
    do iboun = 1,a%nboun
        iBD => BDS(iboun)
        !Deallocate
        deallocate( iBD%cartb, &
            iBD%baloc    ,&
            iBD%eucta    )
    enddo

end subroutine
