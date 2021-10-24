module  Mod_DistributedContainer
   use Mod_DCHashCharSize
   use typre
   implicit none
   private
   public :: DistributedContainer, DistributedContainer_const

   type DistributedContainer
     !class(DistributedContainer), pointer :: next => null()
     
     character(DCHashCharSize) :: key
     character(DCHashCharSize) :: nick
     logical :: isSetKey = .false.
     logical :: isSetNick= .false.
contains

      procedure :: SetKey
      procedure :: GetKey
      procedure :: SetNick
      procedure :: GetNick
   
   end type DistributedContainer

   interface DistributedContainer_const
      procedure constructor   
   end interface

contains

   function constructor()
     class(DistributedContainer), pointer :: constructor
     
     allocate(constructor)
   end function constructor

   subroutine SetNick(this,nick)
      class(DistributedContainer) :: this
      character(*) :: nick
      
      this%nick = nick
      this%isSetNick = .true.
   end
   
   subroutine GetNick(this,nick)
      class(DistributedContainer) :: this
      character(*) :: nick
      
      nick = this%nick
   end

   subroutine SetKey(this,key)
      class(DistributedContainer) :: this
      character(*) :: key
      
      this%key = key
      this%isSetKey = .true.
   end
   
   subroutine GetKey(this,key)
      class(DistributedContainer) :: this
      character(*) :: key
      
      key = this%key
   end
   
  
   

end module Mod_DistributedContainer
