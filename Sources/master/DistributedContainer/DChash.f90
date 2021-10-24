MODULE Mod_DChash
   use Mod_DCHashCharSize
   use Mod_DistributedContainer
   

  IMPLICIT NONE ! Use strong typing
  INTEGER, PARAMETER :: tbl_size = 50



  TYPE sllist
     TYPE(sllist), POINTER :: child => NULL()
     character(DCHashCharSize) :: key
     logical :: associatedKey = .false. 
     class(DistributedContainer), pointer :: val => NULL()

   CONTAINS
     PROCEDURE :: put  => put_sll
     PROCEDURE :: get  => get_sll
     PROCEDURE :: free => free_sll
  END TYPE sllist

  TYPE DChash
     TYPE(sllist), DIMENSION(:), ALLOCATABLE :: vec
     INTEGER                                 :: vec_len = 0
     LOGICAL                                 :: is_init = .FALSE.
   CONTAINS
     PROCEDURE :: init => init_DChash
     PROCEDURE :: put  => put_DChash
     PROCEDURE :: get  => get_DChash
     PROCEDURE :: free => free_DChash
     
     
  END TYPE DChash

  PUBLIC :: DChash, DCHashCharSize

  CONTAINS

  RECURSIVE SUBROUTINE put_sll(list,key,val)
    CLASS(sllist),    INTENT(inout) :: list
    CHARACTER(*) :: key
    class(DistributedContainer), target :: val 
    INTEGER                         :: keylen
    
    integer :: i

    keylen = LEN(key)
    !vallen = LEN(val)
    IF (list%associatedKey) THEN
       IF (list%key /= key) THEN
          IF ( .NOT. ASSOCIATED(list%child) ) ALLOCATE(list%child)
          CALL put_sll(list%child,key,val)
       END IF
    ELSE
        IF (.NOT. list%associatedKey) then
            list%associatedKey = .true.
        END IF
       list%key = key
       list%val => val
    END IF
  END SUBROUTINE put_sll


  RECURSIVE SUBROUTINE get_sll(list,key,val)
    !CLASS(sllist),                 INTENT(in)    :: list
    CLASS(sllist)                     :: list
    CHARACTER(*),              INTENT(in)    :: key
    class(DistributedContainer), pointer :: val 
    integer :: i
    
    IF (list%associatedKey) THEN
        IF (list%key == key) THEN
            val => list%val
        ELSE IF(ASSOCIATED(list%child)) THEN ! keep going
           CALL get_sll(list%child,key,val)
        ELSE
           val => NULL()!exit status
           RETURN
        ENDIF
    ELSE ! At the end of the list, no key found
        IF (ASSOCIATED(list%val)) THEN
            list%val => NULL()!exit status
        END IF
        val => NULL()!exit status
        RETURN
    END IF
  END SUBROUTINE get_sll


  RECURSIVE SUBROUTINE free_sll(list)
    CLASS(sllist), INTENT(inout) :: list
    IF (ASSOCIATED(list%child)) THEN
       CALL free_sll(list%child)
       DEALLOCATE(list%child)
    END IF
    list%child => NULL()
    IF (ASSOCIATED(list%val)) list%val => NULL()
  END SUBROUTINE free_sll

  SUBROUTINE init_DChash(tbl,tbl_len)
    CLASS(DChash),   INTENT(inout) :: tbl
    INTEGER,     OPTIONAL, INTENT(in)    :: tbl_len

    IF (ALLOCATED(tbl%vec)) DEALLOCATE(tbl%vec)
    IF (PRESENT(tbl_len)) THEN
       ALLOCATE(tbl%vec(0:tbl_len-1))
       tbl%vec_len = tbl_len
    ELSE
       ALLOCATE(tbl%vec(0:tbl_size-1))
       tbl%vec_len = tbl_size
    END IF
    tbl%is_init = .TRUE.
  END SUBROUTINE init_DChash

  ! The first part of the hashing procedure using the string
  ! collating sequence
  FUNCTION sum_string(str) RESULT(sig)
    CHARACTER, INTENT(in)   :: str(:)
    INTEGER                        :: sig
    INTEGER :: i

   
    sig = SUM(ICHAR(str))
  END FUNCTION sum_string


  SUBROUTINE put_DChash(tbl,key,val)
    CLASS(DChash), INTENT(inout) :: tbl
    CHARACTER(*),   INTENT(in)    :: key
    class(DistributedContainer), target :: val 
    INTEGER                            :: hash

    character :: akey(len(key))
    
    akey = Copy_s2a(key)

    hash = MOD(sum_string(akey),tbl%vec_len)
    CALL tbl%vec(hash)%put(key=key,val=val)
  END SUBROUTINE put_DChash

  PURE FUNCTION Copy_s2a(s)  RESULT (a)   ! copy s(1:Clen(s)) to char array
   CHARACTER(*),INTENT(IN) :: s
   CHARACTER :: a(LEN(s))
   INTEGER :: i
   DO i = 1,LEN(s)
      a(i) = s(i:i)
   END DO
  END FUNCTION Copy_s2a


  SUBROUTINE get_DChash(tbl,key,val)
    CLASS(DChash),           INTENT(in)    :: tbl
    CHARACTER(*),        INTENT(in)    :: key
    class(DistributedContainer), pointer :: val 
    !CHARACTER(len=:), ALLOCATABLE, INTENT(out)   :: val
    INTEGER                                      :: hash
    
    character :: akey(len(key))
    
    akey = Copy_s2a(key)
    
    hash = MOD(sum_string(akey),tbl%vec_len)
    CALL tbl%vec(hash)%get(key=key,val=val)
    
    
  END SUBROUTINE get_DChash
  
  


  SUBROUTINE free_DChash(tbl)
    CLASS(DChash), INTENT(inout) :: tbl    
    INTEGER     :: i, low, high

    low  = LBOUND(tbl%vec,dim=1)
    high = UBOUND(tbl%vec,dim=1) 
    IF (ALLOCATED(tbl%vec)) THEN
       DO i=low,high
          CALL tbl%vec(i)%free()
       END DO
       DEALLOCATE(tbl%vec)
    END IF
    tbl%is_init = .FALSE.
  END SUBROUTINE free_DChash

END MODULE Mod_DChash
