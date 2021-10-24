module Mod_plcd_Rotation
   use typre
   use Mod_MPIObject
   use Mod_Memor
   implicit none
   private
   public Rotation

   type, extends(MPIObjectLight), abstract :: Rotation

contains
      procedure(SetNdime), deferred             :: SetNdime
      procedure(SetAngle), deferred             :: SetAngle
      procedure(GetAngle), deferred             :: GetAngle
      procedure(ScatterData), deferred          :: ScatterData
      procedure(Initialize), deferred           :: Initialize
      procedure(Finalize), deferred             :: Finalize
      procedure(Vec_Glo2Loc), deferred          :: Vec_Glo2Loc
      procedure(Vec_Loc2Glo), deferred          :: Vec_Loc2Glo
      procedure(Mat_Glo2Loc), deferred          :: Mat_Glo2Loc
      procedure(Mat_Loc2Glo), deferred          :: Mat_Loc2Glo
      procedure(strain_Glo2Loc), deferred       :: strain_Glo2Loc
      procedure(strain_Loc2Glo), deferred       :: strain_Loc2Glo
      procedure(stress_Glo2Loc), deferred       :: stress_Glo2Loc
      procedure(stress_Loc2Glo), deferred       :: stress_Loc2Glo      
      procedure(CT_Glo2Loc), deferred           :: CT_Glo2Loc
      procedure(CT_Loc2Glo), deferred           :: CT_Loc2Glo


   end type

   abstract interface 
      Subroutine SetNdime(a,ndime)
        import
        implicit none
        class(Rotation) :: a
        integer(ip) :: ndime
      end subroutine
      
      Subroutine SetAngle(a,angle,Memor)
        import
        implicit none
        class(Rotation) :: a
        real(rp) :: angle(:)
        type(MemoryMan) :: Memor
      end subroutine
      
      subroutine GetAngle(a,angle)
         import
         class(Rotation), target :: a
         real(rp), pointer :: angle(:)
      end subroutine
      
      subroutine ScatterData(a)
         import
         implicit none
         class(Rotation) :: a
      end subroutine
     
      Subroutine Initialize(a,Memor)
        import
        implicit none
        class(Rotation) :: a
        type(MemoryMan) :: Memor
      end subroutine
      
      Subroutine Finalize(a,Memor)
        import
        implicit none
        class(Rotation) :: a
        type(MemoryMan) :: Memor
      end subroutine      
    
      
      subroutine Vec_Glo2Loc(a,VecGloIn,VecLocOut)
         import
         class(Rotation) :: a
         real(rp) :: VecGloIn(:), VecLocOut(:)
         
      end subroutine
   
      subroutine Vec_Loc2Glo(a,VecLocIn,VecGloOut)
         import
         class(Rotation) :: a
         real(rp) :: VecLocIn(:), VecGloOut(:)
         
      end subroutine
      
      subroutine Mat_Glo2Loc(a,MatGloIn,MatLocOut)
         import
         class(Rotation) :: a
         real(rp) :: MatGloIn(:,:), MatLocOut(:,:)
         
      end subroutine
   
      subroutine Mat_Loc2Glo(a,MatLocIn,MatGloOut)
         import
         class(Rotation) :: a
         real(rp) :: MatLocIn(:,:), MatGloOut(:,:)
         
      end subroutine

      subroutine strain_Glo2Loc(a,strainGloIn,strainLocOut)
         import
         class(Rotation) :: a
         real(rp) :: strainGloIn(:), strainLocOut(:)
         
      end subroutine
   
      subroutine strain_Loc2Glo(a,strainLocIn,strainGloOut)
         import
         class(Rotation) :: a
         real(rp) :: strainLocIn(:), strainGloOut(:)
         
      end subroutine
      
      subroutine stress_Glo2Loc(a,stressGloIn,stressLocOut)
         import
         class(Rotation) :: a
         real(rp) :: stressGloIn(:), stressLocOut(:)
         
      end subroutine
   
      subroutine stress_Loc2Glo(a,stressLocIn,stressGloOut)
         import
         class(Rotation) :: a
         real(rp) :: stressLocIn(:), stressGloOut(:)
         
      end subroutine
      
      subroutine CT_Glo2Loc(a,CTGloIn,CTLocOut)
         import
         class(Rotation) :: a
         real(rp) :: CTGloIn(:,:), CTLocOut(:,:)
         
      end subroutine
   
      subroutine CT_Loc2Glo(a,CTLocIn,CTGloOut)
         import
         class(Rotation) :: a
         real(rp) :: CTLocIn(:,:), CTGloOut(:,:)
         
      end subroutine
   
   end interface
   
end module
