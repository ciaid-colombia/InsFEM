module Mod_plcd_Rotators
   use typre
   use Mod_plcd_Rotation
   use Mod_Memor
   implicit none
   private
   public Rotator_Yes, Rotator_Not
   
   type, extends(Rotation) :: Rotator_Yes
   
      integer(ip) :: ndime
      real(rp), allocatable :: angle(:)
      real(rp), allocatable :: OpeA(:,:), OpeR(:,:), OpeRbar(:,:) 

contains

      procedure :: SetNdime => Yes_SetNdime
      procedure :: SetAngle => Yes_SetAngle
      procedure :: GetAngle => Yes_GetAngle
      procedure :: ScatterData => Yes_ScatterData
      procedure :: Initialize => Yes_Initialize
      procedure :: Finalize => Yes_Finalize
      procedure :: Vec_Glo2Loc => Yes_Vec_Glo2Loc
      procedure :: Vec_Loc2Glo => Yes_Vec_Loc2Glo
      procedure :: Mat_Glo2Loc => Yes_Mat_Glo2Loc
      procedure :: Mat_Loc2Glo => Yes_Mat_Loc2Glo
      procedure :: strain_Glo2Loc => Yes_strain_Glo2Loc
      procedure :: strain_Loc2Glo => Yes_strain_Loc2Glo
      procedure :: stress_Glo2Loc => Yes_stress_Glo2Loc
      procedure :: stress_Loc2Glo => Yes_stress_Loc2Glo
      procedure :: CT_Glo2Loc => Yes_CT_Glo2Loc
      procedure :: CT_Loc2Glo => Yes_CT_Loc2Glo
      
   end type

   type, extends(Rotation) :: Rotator_Not

contains

      procedure :: SetNdime => Not_SetNdime
      procedure :: SetAngle => Not_SetAngle
      procedure :: GetAngle => Not_GetAngle
      procedure :: ScatterData => Not_ScatterData
      procedure :: Initialize => Not_Initialize
      procedure :: Finalize => Not_Finalize
      procedure :: Vec_Glo2Loc => Not_Vec_Glo2Loc
      procedure :: Vec_Loc2Glo => Not_Vec_Loc2Glo
      procedure :: Mat_Glo2Loc => Not_Mat_Glo2Loc
      procedure :: Mat_Loc2Glo => Not_Mat_Loc2Glo
      procedure :: strain_Glo2Loc => Not_strain_Glo2Loc
      procedure :: strain_Loc2Glo => Not_strain_Loc2Glo
      procedure :: stress_Glo2Loc => Not_stress_Glo2Loc
      procedure :: stress_Loc2Glo => Not_stress_Loc2Glo
      procedure :: CT_Glo2Loc => Not_CT_Glo2Loc
      procedure :: CT_Loc2Glo => Not_CT_Loc2Glo
      
   end type

contains

!  SetNdime

   subroutine Yes_SetNdime(a,ndime)
      implicit none
      class(Rotator_Yes) :: a
      integer(ip) :: ndime

      a%ndime = ndime

   end subroutine

   subroutine Not_SetNdime(a,ndime)
      implicit none
      class(Rotator_Not) :: a
      integer(ip) :: ndime

   end subroutine
   
!  SetAngle

   subroutine Yes_SetAngle(a,angle,Memor)
      implicit none
      class(Rotator_Yes) :: a
      real(rp) :: angle(:)
      type(MemoryMan) :: Memor
     
      select case(a%ndime)
      case(2)
!         allocate(a%angle(1))
         Call Memor%alloc(1,a%angle,'EulerAngles','Rotator')
         a%angle = angle(1)
         
      case(3)
!         allocate(a%angle(3))
         Call Memor%alloc(3,a%angle,'EulerAngles','Rotator')
         a%angle = angle(1:3)
         
      case default 
         call runend('Rotation not implemented for this ndime data')
      
      end select

   end subroutine

   subroutine Not_SetAngle(a,angle,Memor)
      implicit none
      class(Rotator_Not) :: a
      real(rp) :: angle(:)
      type(MemoryMan) :: Memor
      
   end subroutine

!  GetAngle
   
   subroutine Yes_GetAngle(a,angle)
      implicit none
      class(Rotator_Yes),target :: a
      real(rp), pointer :: angle(:)

      angle => a%angle
   end subroutine

   subroutine Not_GetAngle(a,angle)
      implicit none
      class(Rotator_Not),target :: a
      real(rp), pointer :: angle(:)

   end subroutine
   
! ScatterData
   
   subroutine Yes_ScatterData(a)
      use MPI
      implicit none
      class(Rotator_Yes) :: a
      integer(ip) :: ierr
      call MPI_BCAST(a%angle,size(a%angle),MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      
   end subroutine
   
   subroutine Not_ScatterData(a)
      use MPI
      implicit none
      class(Rotator_Not) :: a
      
   end subroutine
   

!  Initialize

   subroutine Yes_Initialize(a,Memor)
      implicit none
      class(Rotator_Yes) :: a
      type(MemoryMan) :: Memor
      real(rp), parameter :: PI = 3.14159265
      real(rp) :: c1,s1,c2,s2,c3,s3
      real(rp) :: l1,l2,l3,m1,m2,m3,n1,n2,n3

      select case(a%ndime)
      case(2)
!         allocate(a%OpeA(2,2),a%OpeR(3,3),a%OpeRbar(3,3))
         call Memor%alloc(2,2,a%OpeA,'OpeA','Rotator Initialize')
         call Memor%alloc(3,3,a%OpeR,'OpeR','Rotator Initialize')
         call Memor%alloc(3,3,a%OpeRbar,'OpeRbar','Rotator Initialize')
         
         l1 = dcos(a%angle(1)*PI/180.) ; l2 = -dsin(a%angle(1)*PI/180.)
         m1 = dsin(a%angle(1)*PI/180.) ; m2 = dcos(a%angle(1)*PI/180.)
         
         a%OpeA(1,1) = l1 ;  a%OpeA(1,2) = l2
         a%OpeA(2,1) = m1 ;  a%OpeA(2,2) = m2
         
         a%OpeR(1,1) = l1**2  ;  a%OpeR(1,2) = l2**2  ;  a%OpeR(1,3) = 2*l1*l2
         a%OpeR(2,1) = m1**2  ;  a%OpeR(2,2) = m2**2  ;  a%OpeR(2,3) = 2*m1*m2
         a%OpeR(3,1) = l1*m1  ;  a%OpeR(3,2) = l2*m2 ;  a%OpeR(3,3) = l1*m2+m1*l2

         a%OpeRbar(1,1) = l1**2     ;  a%OpeRbar(1,2) = l2**2    ;  a%OpeRbar(1,3) = l1*l2
         a%OpeRbar(2,1) = m1**2     ;  a%OpeRbar(2,2) = m2**2    ;  a%OpeRbar(2,3) = m1*m2
         a%OpeRbar(3,1) = 2*l1*m1   ;  a%OpeRbar(3,2) = 2*l2*m2  ;  a%OpeRbar(3,3) = l1*m2+m1*l2
         
         
      case(3)
!         allocate(a%OpeA(3,3),a%OpeR(6,6),a%OpeRbar(6,6))
         call Memor%alloc(3,3,a%OpeA,'OpeA','Rotator Initialize')
         call Memor%alloc(6,6,a%OpeR,'OpeR','Rotator Initialize')
         call Memor%alloc(6,6,a%OpeRbar,'OpeRbar','Rotator Initialize')
         
         c1 = dcos(a%angle(1)*PI/180.) ; s1 = dsin(a%angle(1)*PI/180.)
         c2 = dcos(a%angle(2)*PI/180.) ; s2 = dsin(a%angle(2)*PI/180.)      
         c3 = dcos(a%angle(3)*PI/180.) ; s3 = dsin(a%angle(3)*PI/180.)
         
         l1 = c1*c3 - s1*c2*s3   ;  l2 = -c1*s3 - s1*c2*c3  ;  l3 = -s1*s2
         m1 = s1*c3 + c1*c2*s3   ;  m2 = -s1*s3 + c1*c2*c3  ;  m3 = -c1*s2
         n1 = s2*s3              ;  n2 = s2*c3              ;  n3 = c2
                  
         a%OpeA(1,1) = l1 ;  a%OpeA(1,2) = l2 ;  a%OpeA(1,3) = l3
         a%OpeA(2,1) = m1 ;  a%OpeA(2,2) = m2 ;  a%OpeA(2,3) = m3
         a%OpeA(3,1) = n1 ;  a%OpeA(3,2) = n2 ;  a%OpeA(3,3) = n3

         
         a%OpeR(1,1) = l1**2    ; a%OpeR(1,2) = l2**2     ;  a%OpeR(1,3) = l3**2
         a%OpeR(1,4) = 2*l1*l2  ; a%OpeR(1,5) = 2*l1*l3   ;  a%OpeR(1,6) = 2*l2*l3
         
         a%OpeR(2,1) = m1**2    ; a%OpeR(2,2) = m2**2     ;  a%OpeR(2,3) = m3**2
         a%OpeR(2,4) = 2*m1*m2  ; a%OpeR(2,5) = 2*m1*m3   ;  a%OpeR(2,6) = 2*m2*m3
         
         a%OpeR(3,1) = n1**2    ; a%OpeR(3,2) = n2**2     ;  a%OpeR(3,3) = n3**2
         a%OpeR(3,4) = 2*n1*n2  ; a%OpeR(3,5) = 2*n1*n3   ;  a%OpeR(3,6) = 2*n2*n3
         
         a%OpeR(4,1) = l1*m1          ; a%OpeR(4,2) = l2*m2           ;  a%OpeR(4,3) = l3*m3
         a%OpeR(4,4) = l1*m2 + l2*m1  ; a%OpeR(4,5) = l1*m3 + l3*m1   ;  a%OpeR(4,6) = l2*m3 + l3*m2
         
         a%OpeR(5,1) = l1*n1          ; a%OpeR(5,2) = l2*n2           ;  a%OpeR(5,3) = l3*n3
         a%OpeR(5,4) = l1*n2 + l2*n1  ; a%OpeR(5,5) = l1*n3 + l3*n1   ;  a%OpeR(5,6) = l2*n3 + l3*n2
                  
         a%OpeR(6,1) = m1*n1          ; a%OpeR(6,2) = m2*n2           ;  a%OpeR(6,3) = m3*n3
         a%OpeR(6,4) = m1*n2 + m2*n1  ; a%OpeR(6,5) = m1*n3 + m3*n1   ;  a%OpeR(6,6) = m2*n3 + m3*n2

         
         a%OpeRbar(1,1) = l1**2  ; a%OpeRbar(1,2) = l2**2   ;  a%OpeRbar(1,3) = l3**2
         a%OpeRbar(1,4) = l1*l2  ; a%OpeRbar(1,5) = l1*l3   ;  a%OpeRbar(1,6) = l2*l3
         
         a%OpeRbar(2,1) = m1**2  ; a%OpeRbar(2,2) = m2**2   ;  a%OpeRbar(2,3) = m3**2
         a%OpeRbar(2,4) = m1*m2  ; a%OpeRbar(2,5) = m1*m3   ;  a%OpeRbar(2,6) = m2*m3
         
         a%OpeRbar(3,1) = n1**2  ; a%OpeRbar(3,2) = n2**2   ;  a%OpeRbar(3,3) = n3**2
         a%OpeRbar(3,4) = n1*n2  ; a%OpeRbar(3,5) = n1*n3   ;  a%OpeRbar(3,6) = n2*n3
         
         a%OpeRbar(4,1) = 2*l1*m1        ; a%OpeRbar(4,2) = 2*l2*m2         ;  a%OpeRbar(4,3) = 2*l3*m3
         a%OpeRbar(4,4) = l1*m2 + l2*m1  ; a%OpeRbar(4,5) = l1*m3 + l3*m1   ;  a%OpeRbar(4,6) = l2*m3 + l3*m2
         
         a%OpeRbar(5,1) = 2*l1*n1        ; a%OpeRbar(5,2) = 2*l2*n2         ;  a%OpeRbar(5,3) = 2*l3*n3
         a%OpeRbar(5,4) = l1*n2 + l2*n1  ; a%OpeRbar(5,5) = l1*n3 + l3*n1   ;  a%OpeRbar(5,6) = l2*n3 + l3*n2
                  
         a%OpeRbar(6,1) = 2*m1*n1        ; a%OpeRbar(6,2) = 2*m2*n2         ;  a%OpeRbar(6,3) = 2*m3*n3
         a%OpeRbar(6,4) = m1*n2 + m2*n1  ; a%OpeRbar(6,5) = m1*n3 + m3*n1   ;  a%OpeRbar(6,6) = m2*n3 + m3*n2         
         
      case default 
         call runend('Rotation not implemented for this ndime data')
      end select
   end subroutine
      
   subroutine Not_Initialize(a,Memor)
      implicit none
      class(Rotator_Not) :: a
      type(MemoryMan) :: Memor
   
   end subroutine
   
! Finalize

   subroutine Yes_Finalize(a,Memor)
      implicit none
      class(Rotator_Yes) :: a
      type(MemoryMan) :: Memor
      
      select case(a%ndime)
      case(2)
!         allocate(a%OpeA(2,2),a%OpeR(3,3),a%OpeRbar(3,3))
         call Memor%dealloc(2,2,a%OpeA,'OpeA','Rotator Initialize')
         call Memor%dealloc(3,3,a%OpeR,'OpeR','Rotator Initialize')
         call Memor%dealloc(3,3,a%OpeRbar,'OpeRbar','Rotator Initialize')
        
      case(3)
!         allocate(a%OpeA(3,3),a%OpeR(6,6),a%OpeRbar(6,6))
         call Memor%dealloc(3,3,a%OpeA,'OpeA','Rotator Initialize')
         call Memor%dealloc(6,6,a%OpeR,'OpeR','Rotator Initialize')
         call Memor%dealloc(6,6,a%OpeRbar,'OpeRbar','Rotator Initialize')
      
      end select
   
   end subroutine
   
   subroutine Not_Finalize(a,Memor)
      implicit none
      class(Rotator_Not) :: a
      type(MemoryMan) :: Memor
   
   end subroutine
 
!  Vec_Glo2Loc

   subroutine Yes_Vec_Glo2Loc(a,VecGloIn,VecLocOut)
      implicit none
      class(Rotator_Yes) :: a
      real(rp) :: VecGloIn(:), VecLocOut(:)
      
      VecLocOut = matmul(transpose(a%OpeA),VecGloIn)
   
   end subroutine
   
   subroutine Not_Vec_Glo2Loc(a,VecGloIn,VecLocOut)
      implicit none
      class(Rotator_Not) :: a
      real(rp) :: VecGloIn(:), VecLocOut(:)

      VecLocOut = VecGloIn
      
   end subroutine   

!  Vec_Loc2Glo

   subroutine Yes_Vec_Loc2Glo(a,VecLocIn,VecGloOut)
      implicit none
      class(Rotator_Yes) :: a
      real(rp) :: VecLocIn(:), VecGloOut(:)
      
      VecGloOut = matmul(a%OpeA,VecLocIn)
   
   end subroutine
   
   subroutine Not_Vec_Loc2Glo(a,VecLocIn,VecGloOut)
      implicit none
      class(Rotator_Not) :: a
      real(rp) :: VecLocIn(:), VecGloOut(:)

      VecGloOut = VecLocIn
      
   end subroutine  
   
!  Mat_Glo2Loc

   subroutine Yes_Mat_Glo2Loc(a,MatGloIn,MatLocOut)
      implicit none
      class(Rotator_Yes) :: a
      real(rp) :: MatGloIn(:,:), MatLocOut(:,:)
      
      MatLocOut = matmul(transpose(a%OpeA),matmul(MatGloIn,a%OpeA))
   
   end subroutine
   
   subroutine Not_Mat_Glo2Loc(a,MatGloIn,MatLocOut)
      implicit none
      class(Rotator_Not) :: a
      real(rp) :: MatGloIn(:,:), MatLocOut(:,:)

      MatLocOut = MatGloIn
      
   end subroutine 

!  Mat_Loc2Glo

   subroutine Yes_Mat_Loc2Glo(a,MatLocIn,MatGloOut)
      implicit none
      class(Rotator_Yes) :: a
      real(rp) :: MatLocIn(:,:), MatGloOut(:,:)
      
      MatGloOut = matmul(a%OpeA,matmul(MatLocIn,transpose(a%OpeA)))
   
   end subroutine
   
   subroutine Not_Mat_Loc2Glo(a,MatLocIn,MatGloOut)
      implicit none
      class(Rotator_Not) :: a
      real(rp) :: MatLocIn(:,:), MatGloOut(:,:)

      MatGloOut = MatLocIn
      
   end subroutine
   
!  strain_Glo2Loc

   subroutine Yes_strain_Glo2Loc(a,strainGloIn,strainLocOut)
      implicit none
      class(Rotator_Yes) :: a
      real(rp) :: strainGloIn(:), strainLocOut(:)
      
      strainLocOut = matmul(transpose(a%OpeR),strainGloIn)
   
   end subroutine
   
   subroutine Not_strain_Glo2Loc(a,strainGloIn,strainLocOut)
      implicit none
      class(Rotator_Not) :: a
      real(rp) :: strainGloIn(:), strainLocOut(:)

      strainLocOut = strainGloIn
      
   end subroutine
   
!  strain_Loc2Glo

   subroutine Yes_strain_Loc2Glo(a,strainLocIn,strainGloOut)
      implicit none
      class(Rotator_Yes) :: a
      real(rp) :: strainLocIn(:), strainGloOut(:)
      
      strainGloOut = matmul(a%OpeRbar,strainLocIn)
   
   end subroutine
   
   subroutine Not_strain_Loc2Glo(a,strainLocIn,strainGloOut)
      implicit none
      class(Rotator_Not) :: a
      real(rp) :: strainLocIn(:), strainGloOut(:)

      strainGloOut = strainLocIn
      
   end subroutine

!  stress_Glo2Loc

   subroutine Yes_stress_Glo2Loc(a,stressGloIn,stressLocOut)
      implicit none
      class(Rotator_Yes) :: a
      real(rp) :: stressGloIn(:), stressLocOut(:)
      
      stressLocOut = matmul(transpose(a%OpeRbar),stressGloIn)
   
   end subroutine
   
   subroutine Not_stress_Glo2Loc(a,stressGloIn,stressLocOut)
      implicit none
      class(Rotator_Not) :: a
      real(rp) :: stressGloIn(:), stressLocOut(:)

      stressLocOut = stressGloIn
      
   end subroutine

!  stress_Loc2Glo

   subroutine Yes_stress_Loc2Glo(a,stressLocIn,stressGloOut)
      implicit none
      class(Rotator_Yes) :: a
      real(rp) :: stressLocIn(:), stressGloOut(:)
      
      stressGloOut = matmul(a%OpeR,stressLocIn)
   
   end subroutine
   
   subroutine Not_stress_Loc2Glo(a,stressLocIn,stressGloOut)
      implicit none
      class(Rotator_Not) :: a
      real(rp) :: stressLocIn(:), stressGloOut(:)

      stressGloOut = stressLocIn
      
   end subroutine

!  CT_Glo2Loc

   subroutine Yes_CT_Glo2Loc(a,CTGloIn,CTLocOut)
      implicit none
      class(Rotator_Yes) :: a
      real(rp) :: CTGloIn(:,:), CTLocOut(:,:)
      
      CTLocOut = matmul(transpose(a%OpeRbar),matmul(CTGloIn,a%OpeRbar))
   
   end subroutine
   
   subroutine Not_CT_Glo2Loc(a,CTGloIn,CTLocOut)
      implicit none
      class(Rotator_Not) :: a
      real(rp) :: CTGloIn(:,:), CTLocOut(:,:)

      CTLocOut = CTGloIn
      
   end subroutine

!  CT_Loc2Glo

   subroutine Yes_CT_Loc2Glo(a,CTLocIn,CTGloOut)
      implicit none
      class(Rotator_Yes) :: a
      real(rp) :: CTLocIn(:,:), CTGloOut(:,:)
      
      CTGloOut = matmul(a%OpeR,matmul(CTLocIn,transpose(a%OpeR)))
   
   end subroutine
   
   subroutine Not_CT_Loc2Glo(a,CTLocIn,CTGloOut)
      implicit none
      class(Rotator_Not) :: a
      real(rp) :: CTLocIn(:,:), CTGloOut(:,:)

      CTGloOut = CTLocIn
      
   end subroutine
   
end module
