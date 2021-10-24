subroutine nsi_turnon(a)
    ! DESCRIPTION
    !    This routine inits some stuff for the incompressible NS equations
    !-----------------------------------------------------------------------
    use typre
    use Mod_NavierStokes
    implicit none
    class(NavierStokesProblem) :: a
    integer(ip)::ndime,idime,count,localP,nod
    real(rp), allocatable :: localPointRhs(:,:),aux(:,:)
    logical(lg)::pass
    
    interface
      subroutine nsi_InitializeTurbulentInletBoundaryConditions(a)
         import
         implicit none
         class(NavierStokesProblem) :: a
      end subroutine
      
    end interface

    call a%Mesh%GetNdime(ndime)
    call a%Mesh%GetNpoinLocal(localP)

    allocate(localPointRhs(ndime + 1 ,localP))
    allocate(aux(ndime + 1,localP))

    count=0

    !Find non zero entries
    do nod = 1,localP

        pass =.false.
        if(a%PointwiseSource(1,nod)/=0) then 
            pass=.true.
        endif

        if (pass) then
            localPointRhs(:,count+1) = a%PointwiseSource(:,nod)
            count = count +1
        endif

    enddo

    call a%Memor%dealloc(size(a%PointwiseSource,1),size(a%PointwiseSource,2),a%PointwiseSource,'PointwiseSource','nsi_turnon')
    if(count > 0) then
        a%sizeOfPointWise=count
        a%dimOfPointWise = ndime+1
        call a%Memor%alloc(ndime+1,count,a%PointwiseSource,'PointwiseSource','nsi_turnon')

        a%PointwiseSource(:,:)=localPointRhs(:,1:count)
    else
        a%sizeOfPointWise= 0
        a%dimOfPointWise = ndime

    end if

   if (a%kfl_TurbulentInletBoundaryConditions == 1) call nsi_InitializeTurbulentInletBoundaryConditions(a)
        
end subroutine
