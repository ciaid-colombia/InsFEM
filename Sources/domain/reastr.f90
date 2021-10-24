
subroutine reastr(a)
!-----------------------------------------------------------------------
!****f* Domain/reastr
! NAME
!    reastr
! DESCRIPTION
!    This routine reads the numerical strategy needed to build up
!    the mesh data. The variables read are:      
! 
!    LQUAD   Type of integration rule (1 = closed, 0 = open)
!    NGAUS   Number of domain integration points
!    OUMES   Flag to decide whether to dump the mesh data or not 
!    SCALE   Geometrical a%scale factor
! USED BY
!    Domain
!***
!-----------------------------------------------------------------------
  use typre
  use Mod_iofile
  use Mod_Mesh
  
  implicit none
  class(FemMesh) :: a
  
  !For Listener
   real(rp), pointer     :: param(:) => NULL()
   character(5), pointer :: words(:) => NULL()
   integer(ip), pointer  :: nnpar => NULL()
   integer(ip), pointer  :: nnwor => NULL()
  
  integer(ip) :: istat,ielty,mnode,j
  character(150) :: fil_outpu_dom

  call a%Listener%getarrs(words,param,nnpar,nnwor)
  
  call a%Memor%alloc(a%nelty,a%lquad,'LQUAD','reastr')
  call a%Memor%alloc(a%nelty,a%ngaus,'NGAUS','reastr')
  


  a%mitsm      = 1
  a%tolsm      = 1.0_rp
  a%lquad      = 0                                         !Open integration rule
  a%scale(1)   = 1.0_rp                                    !X a%scale factor
  a%scale(2)   = 1.0_rp                                    !Y a%scale factor
  a%scale(3)   = 1.0_rp                                    !Z a%scale factor
  a%kfl_perio  = 0                                         !Periodic bc's
  a%kfl_exnpo  = 0
  a%kfl_TestCommunications = 0                             !Default is do not test communications
  a%kfl_MeshInputType = 0                                  !Default is from disk
  a%kfl_UseElementDataStructures = .false.
  a%kfl_ReadBlock= .false.
  a%numGeoBlocks = 1
  
!
! Reach the section.
!
  call a%Listener%rewind
  call a%Listener%listen('REASTR')
  do while(words(1)/='STRAT')
     call a%Listener%listen('REASTR')
  end do
!
! Read strategy
!      
  do while(words(1)/='ENDST')
     call a%Listener%listen('REASTR')
     if(words(1)=='SMOOT') then                          ! Smoothing iterations
        a%mitsm = a%Listener%getint('ITERA',1,'#Smoothing iterations')
        a%tolsm = a%Listener%getrea('TOLER',1.0_rp,'#Smoothing tolerance')
     else if(words(1)=='INTEG') then                     ! Integration rule
        if(words(3)=='') then
           if(words(2)=='CLOSE') then
              a%lquad=(/(1,ielty=1,a%nelty)/)
           else
              a%lquad=(/(0,ielty=1,a%nelty)/)              
           end if
        else
           do ielty=1,a%nelty
              if(words(ielty+1)=='CLOSE') then
                 a%lquad(ielty)=1
              else
                 a%lquad(ielty)=0
              end if
           end do
        end if
     else if(words(1)=='DOMAI') then                     ! Domain integration pts
        a%ngaus=(/(int(param(ielty)),ielty=1,a%nelty)/)
     else if(words(1)=='SCALE') then                     ! Geometric a%scale factors
        a%scale(1) = a%Listener%getrea('XSCAL',1.0_rp,'#x-factor')
        a%scale(2) = a%Listener%getrea('YSCAL',1.0_rp,'#y-factor')
        a%scale(3) = a%Listener%getrea('ZSCAL',1.0_rp,'#z-factor') 
     else if(words(1)=='PERIO') then
        if(a%Listener%exists('ON   ')) a%kfl_perio = 1
        a%PerDime(1:3) = int(param(2:4))
        if(a%Listener%exists('FORCE')) a%kfl_ForcePeriodicMesh = 1
        a%ForcePeriodicMeshTol = a%Listener%getrea('TOLER',1e-2_rp,'Tolerance for force periodic mesh')
     else if(words(1)=='READS') then
        if(a%Listener%exists('DEFAU')) then
            a%kfl_MeshInputType = 0
        elseif(a%Listener%exists('MANUF')) then
            a%kfl_MeshInputType = 1  
            a%ManufacturedNpoinX = a%Listener%getint('NPOIN',10,'*NUMBER OF POINTS LOCALLY')
            a%ManufacturedType = a%Listener%words(5)
            if(a%Listener%exists('POLAR')) then
               a%kfl_ManuMeshType = 1  
               a%ManufacturedInternalRadius = a%Listener%getrea('INRAD',0.0_rp,'Internal radius')
               a%ManufacturedExternalRadius = a%Listener%getrea('EXRAD',0.0_rp,'External radius')
               a%ManufacturedDivisionsInAngle = a%Listener%getint('NDIVI',10,'*NUMBER OF ANGULAR DIVISIONS LOCALLY')
            elseif(a%Listener%exists('CYLIN')) then
               a%kfl_ManuMeshType = 2  
               a%ManufacturedInternalRadius = a%Listener%getrea('INRAD',0.0_rp,'Internal radius')
               a%ManufacturedExternalRadius = a%Listener%getrea('EXRAD',1.0_rp,'External radius')
               a%ManufacturedAxialLength = a%Listener%getrea('AXIAL',1.0_rp,'Axial length')
               a%ManufacturedDivisionsInAngle = a%Listener%getint('NDIVI',10,'*NUMBER OF ANGULAR DIVISIONS LOCALLY')
               a%ManufacturedNpoinZ = a%Listener%getint('NPOIZ',10,'*NUMBER OF AXIAL POINTS LOCALLY')
            endif
         endif
     else if(words(1)=='EXNOP') then
        if(a%Listener%exists('ON   ')) a%kfl_exnpo = 1   
     else if(words(1)=='TESTC') then
        if(a%Listener%exists('ON   ')) a%kfl_TestCommunications = 1      
     else if(words(1)=='ELEME') then
         if(a%Listener%exists('ON   ')) a%kfl_UseElementDataStructures = .true.  
     else if(words(1)=='GEOBL') then
         if(a%Listener%exists('ON   ')) then

             a%kfl_ReadBlock= .true.  
             a%numGeoBlocks = int(param(2))
             !call a%Memor%alloc(a%numGeoBlocks,a%blockName,'BlockName','reageo')
             !do j=1,a%numGeoBlocks
             !    a%blockName(j) = a%Listener%words(j+2)
             !enddo
         end if
     end if
     
     
  end do
  
  !Check if Gauss points have been assigned
  do ielty=1,a%nelty
     if(a%ngaus(ielty)<=0) then     ! If a%ngaus=0, assign by default:
        a%ngaus(ielty)=a%nnode(ielty) ! a%ngaus=a%nnode
        a%lquad(ielty)=0            ! open rule

!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!* Añadido selección de cuadraturas en rutope para elementos p1, p2, p3, p4
!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        if (a%ndime==2.and.a%nnode(ielty)==3 ) a%ngaus(ielty)=3   ! P1 elements
        if (a%ndime==2.and.a%nnode(ielty)==6 ) a%ngaus(ielty)=6   ! P2 elements
        if (a%ndime==2.and.a%nnode(ielty)==10) a%ngaus(ielty)=7   ! P3 elements
        if (a%ndime==2.and.a%nnode(ielty)==15) a%ngaus(ielty)=13  ! P4 elements
!* * * * * * * *
!* hasta aqui
!* * * * * * * *

     end if     
  end do
  
  !Compute the list of element types
  call a%Memor%alloc(a%mnode,a%ltype,'LTYPE','reageo')
  a%ltype = 0
  do ielty = 1,a%nelty
     a%ltype(a%nnode(ielty)) = ielty
  enddo


  if(a%kfl_mnelt==0) call sortin(a%nelty,a%ngaus)

end subroutine reastr
